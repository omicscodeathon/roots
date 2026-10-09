#!/usr/bin/env bash
# 1_verify_trajectories.sh -- Run LOCALLY before deploying to Hetzner.
# Installs GROMACS (if missing), runs `gmx check` on every harvest .xtc,
# and prints the true frame count + simulation time for each system.
# This takes ~30 seconds and tells you whether you have 50ns or ~34ns.
set -euo pipefail

HARVEST="${HARVEST:-$HOME/roots_pipeline1/GPUai_harvest}"

echo "=== ROOTS trajectory verification ==="
echo "Harvest dir: $HARVEST"
echo

# ---- ensure gmx is available ----
if ! command -v gmx &>/dev/null; then
  echo "gmx not found -- installing via conda..."
  # Retry-friendly: transient HTTP errors are common on conda-forge CDN
  for attempt in 1 2 3; do
    if conda install -y -c conda-forge gromacs 2>&1; then
      break
    fi
    echo "  conda install attempt $attempt failed, retrying in 10s..."
    sleep 10
  done
  if ! command -v gmx &>/dev/null; then
    echo "ERROR: could not install gromacs after 3 attempts." >&2
    echo "Try: conda install -y -c conda-forge gromacs   (or mamba)" >&2
    exit 1
  fi
fi
echo "Using: $(gmx --version 2>&1 | head -1)"
echo

# ---- check each system ----
for sys in Sys1 Sys2 Sys3 Sys4; do
  xtc=$(find "$HARVEST/$sys" -name '*.xtc' -type f | head -1)
  if [[ -z "$xtc" ]]; then
    echo "[$sys] WARNING: no .xtc found in $HARVEST/$sys"
    continue
  fi
  echo "[$sys] Checking: $xtc ($(du -h "$xtc" | cut -f1))"
  # gmx check prints frame count and last time to stderr
  check_out=$(gmx check -f "$xtc" 2>&1 || true)
  last_time=$(echo "$check_out" | grep -i "last frame" | tail -1 || true)
  n_frames=$(echo "$check_out" | grep -i "frame" | grep -i "step" | tail -1 || true)
  # Fallback: look for "Step" line with the largest time value
  if [[ -z "$last_time" ]]; then
    last_time=$(echo "$check_out" | grep -oP 'time\s+\K[0-9.]+' | sort -n | tail -1 || true)
    [[ -n "$last_time" ]] && last_time="Last time: ${last_time} ps"
  fi
  echo "  $last_time"
  echo "  $n_frames"
  echo
done

echo "=== Verification complete ==="
echo "If all systems show ~50000 ps (50 ns): full trajectories, proceed as planned."
echo "If any show ~33940 ps (34 ns): trajectories stopped early. Still publishable --"
echo "  use frames from the equilibrated portion (typically last 20+ ns) for MM-PBSA."
