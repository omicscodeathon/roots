#!/usr/bin/env bash
# 5_download_and_cleanup.sh -- Download results from Hetzner and delete the server.
#
# Run LOCALLY after steps 3 and 4 have completed on the server.
set -euo pipefail

SERVER_NAME="${SERVER_NAME:-roots-analysis}"
LOCAL_RESULTS="${LOCAL_RESULTS:-$HOME/roots_pipeline1/results}"

need_bin() { command -v "$1" >/dev/null || { echo "'$1' not found. $2" >&2; exit 1; }; }
need_bin hcloud "Install: brew install hcloud"
need_bin rsync  "Install rsync."

IP=$(hcloud server ip "$SERVER_NAME" 2>/dev/null || true)
if [[ -z "$IP" ]]; then
  echo "Server '$SERVER_NAME' not found. Already deleted?" >&2
  exit 1
fi

echo "Downloading results from $SERVER_NAME ($IP)..."
mkdir -p "$LOCAL_RESULTS"

rsync -avz --progress root@"$IP":/root/roots/results/ "$LOCAL_RESULTS/"

echo
echo "Results downloaded to: $LOCAL_RESULTS/"
echo
echo "Contents:"
find "$LOCAL_RESULTS" -type f | head -40
echo

# ---- Confirm before deleting ----
read -rp "Delete server '$SERVER_NAME' ($IP)? This stops billing. [y/N] " confirm
if [[ "$confirm" =~ ^[Yy] ]]; then
  hcloud server delete "$SERVER_NAME"
  echo "Server deleted. Billing stopped."
else
  echo "Server NOT deleted. It's still running and billing at ~\$0.07/hr."
  echo "To delete later: hcloud server delete $SERVER_NAME"
fi

echo
echo "=== ROOTS analysis pipeline complete ==="
echo
echo "Your results:"
echo "  $LOCAL_RESULTS/mmpbsa/     -- MM-PBSA binding free energies (GB & PB)"
echo "  $LOCAL_RESULTS/analysis/   -- trajectory analysis plots & data"
echo
echo "Key files for your manuscript:"
echo "  mmpbsa/Sys*/FINAL_RESULTS_MMPBSA_*.dat  -- binding DeltaG values"
echo "  mmpbsa/Sys*/FINAL_DECOMP_MMPBSA_*.csv   -- per-residue decomposition"
echo "  analysis/comparative/summary_statistics.csv  -- all-system summary table"
echo "  analysis/comparative/comparative_rmsd.png    -- Figure: RMSD overlay"
echo "  analysis/comparative/comparative_rmsf.png    -- Figure: RMSF overlay"
echo "  analysis/Sys*/pca.png                        -- Figure: PCA per system"
echo "  analysis/Sys*/hbonds_occupancy.png            -- Figure: H-bond persistence"
