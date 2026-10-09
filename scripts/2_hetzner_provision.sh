#!/usr/bin/env bash
# 2_hetzner_provision.sh -- Create a Hetzner CCX23 instance, upload data, install deps.
#
# Prerequisites (on your local machine):
#   export HCLOUD_TOKEN=your_read_write_token   # needs WRITE access for server creation
#   brew install hcloud  /  sudo apt install hcloud-cli
#   ssh key registered in Hetzner Cloud Console (Project > Security > SSH Keys)
#
# This creates a CCX23 (4 dedicated AMD vCPU, 16 GB RAM, ~EUR 0.0641/hr = ~$0.07/hr).
# Estimated total runtime: ~8h for steps 2+3 -> ~$0.56 compute + negligible egress.
set -euo pipefail

# ---- configuration ----
SERVER_NAME="${SERVER_NAME:-roots-analysis}"
SERVER_TYPE="${SERVER_TYPE:-ccx23}"
IMAGE="${IMAGE:-ubuntu-24.04}"
LOCATION="${LOCATION:-nbg1}"        # Nuremberg; alternatives: fsn1 (Falkenstein), hel1 (Helsinki)
SSH_KEY_NAME="${SSH_KEY_NAME:-}"     # Name of your SSH key in Hetzner Console; leave blank to auto-detect

HARVEST="${HARVEST:-$HOME/roots_pipeline1/GPUai_harvest}"
TOPO_BASE="${TOPO_BASE:-$HOME/roots_pipeline1}"

# Topology paths (CHARMM-GUI output directories)
declare -A TOPO_DIRS
TOPO_DIRS[Sys1]="$TOPO_BASE/Sys1+Sys2_results/Sys1/charmm-gui-9082326193/gromacs"
TOPO_DIRS[Sys2]="$TOPO_BASE/Sys1+Sys2_results/Sys2/charmm-gui-9086874719/gromacs"
TOPO_DIRS[Sys3]="$TOPO_BASE/Sys3+Sys4_results/Sys3/charmm-gui-9087526694/gromacs"
TOPO_DIRS[Sys4]="$TOPO_BASE/Sys3+Sys4_results/Sys4/charmm-gui-9087521260/gromacs"

need_bin() { command -v "$1" >/dev/null || { echo "'$1' not found. $2" >&2; exit 1; }; }
need_bin hcloud "Install: brew install hcloud OR snap install hcloud"
need_bin rsync  "Install rsync."
need_bin ssh    "Install openssh-client."

# ---- detect SSH key ----
if [[ -z "$SSH_KEY_NAME" ]]; then
  SSH_KEY_NAME=$(hcloud ssh-key list -o noheader -o columns=name | head -1)
  [[ -n "$SSH_KEY_NAME" ]] || { echo "No SSH keys found in Hetzner Console. Add one first." >&2; exit 1; }
  echo "Auto-detected SSH key: $SSH_KEY_NAME"
fi

# ---- create server ----
echo "Creating $SERVER_TYPE server '$SERVER_NAME' in $LOCATION..."
hcloud server create \
  --name "$SERVER_NAME" \
  --type "$SERVER_TYPE" \
  --image "$IMAGE" \
  --location "$LOCATION" \
  --ssh-key "$SSH_KEY_NAME"

IP=$(hcloud server ip "$SERVER_NAME")
echo "Server IP: $IP"
echo "Waiting for SSH to become available..."
for i in $(seq 1 30); do
  ssh -o StrictHostKeyChecking=no -o ConnectTimeout=5 root@"$IP" true 2>/dev/null && break
  sleep 5
done

# ---- install dependencies on server ----
echo "Installing analysis stack on $SERVER_NAME..."
ssh root@"$IP" bash <<'REMOTE_SETUP'
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

apt-get update -qq
apt-get install -y -qq python3-pip python3-venv gromacs libgromacs-dev \
  build-essential cmake pkg-config libfftw3-dev liblapack-dev

# Create analysis virtualenv
python3 -m venv /opt/roots_env
source /opt/roots_env/bin/activate

pip install --upgrade pip
pip install \
  gmx_MMPBSA==1.6.3 \
  MDAnalysis==2.8.0 \
  plip==2.4.1 \
  matplotlib==3.9.3 \
  seaborn==0.13.2 \
  scipy==1.14.1 \
  scikit-learn==1.6.0 \
  py3Dmol==2.4.2 \
  pandas==2.2.3 \
  openbabel-wheel \
  parmed

echo "Dependencies installed successfully."
gmx --version 2>&1 | head -1
python3 -c "import MDAnalysis; print(f'MDAnalysis {MDAnalysis.__version__}')"
REMOTE_SETUP

# ---- upload data ----
echo "Creating workspace on server..."
ssh root@"$IP" "mkdir -p /root/roots/{harvest,topo}/{Sys1,Sys2,Sys3,Sys4} /root/roots/results"

echo "Uploading harvest data (trajectories, topologies, checkpoints)..."
for sys in Sys1 Sys2 Sys3 Sys4; do
  echo "  $sys harvest..."
  rsync -az --progress "$HARVEST/$sys/" root@"$IP":/root/roots/harvest/"$sys"/

  echo "  $sys topology..."
  rsync -az --progress "${TOPO_DIRS[$sys]}/" root@"$IP":/root/roots/topo/"$sys"/
done

echo
echo "============================================"
echo "Hetzner server ready: $SERVER_NAME @ $IP"
echo "Server type: $SERVER_TYPE (~\$0.07/hr)"
echo "SSH: ssh root@$IP"
echo "Workspace: /root/roots/"
echo "============================================"
echo
echo "Next: scp the analysis scripts to the server and run them:"
echo "  scp 3_run_mmpbsa.sh 4_run_analysis.py root@$IP:/root/roots/"
echo "  ssh root@$IP 'cd /root/roots && bash 3_run_mmpbsa.sh'"
echo "  ssh root@$IP 'cd /root/roots && source /opt/roots_env/bin/activate && python3 4_run_analysis.py'"
echo
echo "When done, download results and delete the server:"
echo "  rsync -az root@$IP:/root/roots/results/ ./roots_results/"
echo "  hcloud server delete $SERVER_NAME"
