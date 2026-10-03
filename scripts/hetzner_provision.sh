#!/usr/bin/env bash
# Run this on YOUR machine (needs the `hcloud` CLI installed + authenticated:
# `hcloud context create roots` and paste an API token from the Hetzner
# Cloud console -> Security -> API Tokens).
#
# FRESH-START VERSION: intentionally does NOT upload the 99,901 ligands
# already sitting in step4/ligands locally. The remote box starts at
# 0/328869 and ligand_prep.py reprocesses every compound from scratch.
# (Your local 99,901 files are left untouched on disk -- delete them
# yourself later if you want the space back.)
#
# It will:
#   1. Spin up a CCX63 (48 vCPU / 188 GB, dedicated AMD EPYC) server
#   2. Install Docker on it
#   3. Upload your step3/ input data, Dockerfile/build context, and the
#      two scripts (NOT any existing step4/ligands -- see above)
#   4. Build the roots-ligand-prep image on the box
#   5. Kick off run_with_notify.sh from zero, auto-restarting on crash
#
# Fill in the paths below first. This connects as root@<server-ip>, so
# everything lands under /root/roots_pipeline1 on the remote box, even
# though the local copy lives under /home/genome-alchemist/roots_pipeline1.

set -euo pipefail

# ---- EDIT THESE ----
SERVER_NAME="step4-ligand-prep"
SSH_KEY_NAME="step4-key"                          # `hcloud ssh-key list` to see names
LOCATION="nbg1"                                   # nbg1 / fsn1 / hel1 (Germany/Finland, cheapest) or ash/hil (US)
LOCAL_STEP3_DIR="$HOME/roots_pipeline1/step3"       # must contain library_final.csv
LOCAL_DOCKER_BUILD_DIR="$HOME/roots_pipeline1"      # directory containing your Dockerfile
LOCAL_LIGAND_PREP_PY="$HOME/roots_pipeline1/ligand_prep.py"
LOCAL_RUN_SCRIPT="$HOME/roots_pipeline1/run_with_notify.sh"
# ---------------------

echo "== Creating $SERVER_NAME (ccx63, $LOCATION) =="
hcloud server create \
  --name "$SERVER_NAME" \
  --type ccx63 \
  --image ubuntu-24.04 \
  --location "$LOCATION" \
  --ssh-key "$SSH_KEY_NAME"

SERVER_IP=$(hcloud server ip "$SERVER_NAME")
echo "Server IP: $SERVER_IP"

echo "== Waiting for SSH =="
until ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=5 "root@${SERVER_IP}" true 2>/dev/null; do
  sleep 5
done
echo "SSH is up."

echo "== Installing Docker on the box =="
ssh "root@${SERVER_IP}" 'curl -fsSL https://get.docker.com | sh'

echo "== Creating clean remote directory layout (fresh start -- no leftover ligands) =="
ssh "root@${SERVER_IP}" 'rm -rf ~/roots_pipeline1/step4 && mkdir -p ~/roots_pipeline1/step3 ~/roots_pipeline1/step4/ligands ~/roots_pipeline1/step4/logs'

echo "== Uploading step3 library (this is the input, needed before anything can run) =="
rsync -az --info=progress2 "${LOCAL_STEP3_DIR}/" "root@${SERVER_IP}:~/roots_pipeline1/step3/"

echo "== Skipping ligand upload on purpose -- reprocessing all 328,869 compounds from scratch =="

echo "== Uploading Docker build context, and the two scripts =="
rsync -az --exclude '.git' "${LOCAL_DOCKER_BUILD_DIR}/Dockerfile" "root@${SERVER_IP}:~/roots_pipeline1/Dockerfile"
scp "$LOCAL_LIGAND_PREP_PY" "root@${SERVER_IP}:~/roots_pipeline1/ligand_prep.py"
scp "$LOCAL_RUN_SCRIPT" "root@${SERVER_IP}:~/roots_pipeline1/run_with_notify.sh"
ssh "root@${SERVER_IP}" 'chmod +x ~/roots_pipeline1/run_with_notify.sh'

echo "== Building roots-ligand-prep image on the box =="
ssh "root@${SERVER_IP}" 'cd ~/roots_pipeline1 && docker build -t roots-ligand-prep .'

echo "== Launching run_with_notify.sh (fresh start, survives you disconnecting) =="
ssh "root@${SERVER_IP}" 'cd ~/roots_pipeline1 && setsid ./run_with_notify.sh < /dev/null > /tmp/launch.log 2>&1 &'
sleep 2
ssh "root@${SERVER_IP}" 'cat /tmp/launch.log'

cat <<EOF

== Done ==
Server:  $SERVER_NAME ($SERVER_IP)
SSH in with:      ssh root@${SERVER_IP}
Check status:     ssh root@${SERVER_IP} 'cat ~/roots_pipeline1/step4/logs/run_status.txt'
Tail the log:     ssh root@${SERVER_IP} 'tail -f ~/roots_pipeline1/step4/logs/docker_run.log'
Count finished:   ssh root@${SERVER_IP} 'ls ~/roots_pipeline1/step4/ligands | wc -l'

Don't forget: delete the server when it's done (hcloud server delete $SERVER_NAME)
-- CCX63 runs about \$1.56/hr (~\$37/day) whether it's working or idle.
EOF
