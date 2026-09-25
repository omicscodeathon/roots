#!/usr/bin/env bash
#
# ROOTS -- Step 6: High-Throughput Docking Worker (Layer 2)
#
# Runs on a single fleet instance. Docks one chunk of the compound library
# against the validated ALDH3B2 receptor using QuickVina2, resumable and
# checkpointed so a mid-run interruption loses at most the in-progress
# compound, never completed work.
#
# QuickVina2 vs. standard Vina: this step deliberately uses QuickVina2, a
# CPU-only, Vina-compatible drop-in with a published average 2.3-fold
# speedup (max ~20-fold) over standard Vina at equal exhaustiveness, and a
# 0.967 correlation to Vina's own binding-energy predictions on the top
# pose (Vina-compatible scoring, not a different, incomparable one). This
# combination -- real speedup with near-identical scoring -- is why
# QuickVina2 is used for this large-scale step while Step 5's validation
# redock deliberately used standard Vina at higher exhaustiveness.
#
# Exhaustiveness is deliberately lower here (8) than Step 5's validation
# redock (32): this is the standard, published practice for large-library
# virtual screening, where the goal is reliable relative ranking across
# ~10^5-10^6 compounds rather than the single highest-precision pose. Any
# compound that becomes a real hit candidate should be re-docked at high
# exhaustiveness individually before being reported as a validated result.
#
# Usage (on each fleet instance, after deployment):
#   bash dock_worker.sh <CHUNK_ID>
# where <CHUNK_ID> is a two-digit string (e.g. "00") matching
# /root/chunk_<CHUNK_ID>.txt, a newline-separated list of ligand filenames
# for this instance's share of the library.

set -e

CHUNK_ID=$1
if [ -z "$CHUNK_ID" ]; then
    echo "Usage: bash dock_worker.sh <CHUNK_ID>"
    exit 1
fi

RECEPTOR="/root/step5_docking/receptors/ALDH3B2_holo_receptor.pdbqt"
LIGAND_DIR="/root/ligands_final"
OUTPUT_DIR="/root/docking_results"
LOG_DIR="/root/docking_logs"
CHUNK_FILE="/root/chunk_${CHUNK_ID}.txt"

# Box parameters -- identical to Step 5's validation run. Never change these
# between validation and production without re-validating: the whole point
# of Step 5 is confirming the docking protocol against THIS box.
CENTER_X=-11.727
CENTER_Y=4.248
CENTER_Z=6.293
SIZE_X=18.7
SIZE_Y=13.8
SIZE_Z=20.3
EXHAUSTIVENESS=8
SEED=42

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

if [ ! -f "$CHUNK_FILE" ]; then
    echo "Chunk file not found: $CHUNK_FILE"
    exit 1
fi

TOTAL=$(wc -l < "$CHUNK_FILE")
DONE=0
FAILED=0

echo "Starting chunk ${CHUNK_ID}: ${TOTAL} compounds"
echo "Receptor: ${RECEPTOR}"
echo "Box: center=(${CENTER_X},${CENTER_Y},${CENTER_Z}) size=(${SIZE_X},${SIZE_Y},${SIZE_Z})"
echo "Exhaustiveness: ${EXHAUSTIVENESS}, seed: ${SEED}"

while read -r compound_filename; do
  compound_id="${compound_filename%.pdbqt}"
  ligand_path="${LIGAND_DIR}/${compound_filename}"
  out_path="${OUTPUT_DIR}/${compound_id}_docked.pdbqt"

  # Resume support: skip anything already docked from a prior, interrupted run
  if [ -f "$out_path" ]; then
    DONE=$((DONE+1))
    continue
  fi

  if /usr/local/bin/qvina2.1 \
      --receptor "$RECEPTOR" \
      --ligand "$ligand_path" \
      --center_x "$CENTER_X" --center_y "$CENTER_Y" --center_z "$CENTER_Z" \
      --size_x "$SIZE_X" --size_y "$SIZE_Y" --size_z "$SIZE_Z" \
      --exhaustiveness "$EXHAUSTIVENESS" \
      --seed "$SEED" \
      --out "$out_path" \
      > "${LOG_DIR}/${compound_id}.log" 2>&1; then
    DONE=$((DONE+1))
  else
    FAILED=$((FAILED+1))
  fi

  # Checkpoint progress every 500 compounds, not just at the end -- a real
  # progress record on disk, not a claim that could be lost with the process.
  if [ $(( (DONE+FAILED) % 500 )) -eq 0 ]; then
    echo "[$(date +%H:%M:%S)] Progress: $((DONE+FAILED))/$TOTAL (done=$DONE, failed=$FAILED)" >> /root/progress.log
  fi
done < "$CHUNK_FILE"

echo "[$(date +%H:%M:%S)] COMPLETE: $((DONE+FAILED))/$TOTAL (done=$DONE, failed=$FAILED)" >> /root/progress.log
touch /root/DOCKING_COMPLETE

echo ""
echo "Chunk ${CHUNK_ID} complete: ${DONE} docked, ${FAILED} failed"
