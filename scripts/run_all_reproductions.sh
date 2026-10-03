#!/usr/bin/env bash
# Batch reproduction test: runs REPRODUCE_sys1-4.sh in an isolated scratch copy,
# never touching the real working directory, and reports pass/fail against the
# confirmed reference values from tonight's actual run.
set -uo pipefail

SOURCE_DIR="$HOME/roots_pipeline1/step7_admet/charmmgui_ligand_inputs"
SCRATCH_DIR="$HOME/roots_pipeline1/step7_admet/REPRO_TEST_$(date +%Y%m%d_%H%M%S)"

echo "=== Setting up isolated scratch copy ==="
echo "Source (untouched): $SOURCE_DIR"
echo "Scratch (test runs here): $SCRATCH_DIR"
mkdir -p "$SCRATCH_DIR"

# Copy only the REAL INPUT files each script needs -- the receptor, the four docked
# PDBQTs, and the fix script -- not tonight's already-generated outputs, so this is
# a genuine from-scratch reproduction, not a no-op re-copy of existing results.
cp "$SOURCE_DIR/ALDH3B2_DOCKING_receptor.pdb" "$SCRATCH_DIR/" || { echo "MISSING: ALDH3B2_DOCKING_receptor.pdb"; exit 1; }
cp "$HOME/fix_ligand_for_charmmgui.py" "$SCRATCH_DIR/../fix_ligand_for_charmmgui.py" 2>/dev/null \
  || cp "$SOURCE_DIR/../fix_ligand_for_charmmgui.py" "$SCRATCH_DIR/../" 2>/dev/null \
  || echo "WARNING: fix_ligand_for_charmmgui.py not found at an expected path -- scripts will fail at Step 1"

declare -A PDBQT_FILES=(
  [sys1]="CNP0220660.1_docked.pdbqt"
  [sys2]="CNP0358193.1_docked.pdbqt"
  [sys3]="CNP0384182.2_docked.pdbqt"
  [sys4]="CNP0088713.1_docked.pdbqt"
)
for sys in sys1 sys2 sys3 sys4; do
  f="${PDBQT_FILES[$sys]}"
  if [ -f "$SOURCE_DIR/$f" ]; then
    cp "$SOURCE_DIR/$f" "$SCRATCH_DIR/"
  else
    echo "MISSING: $f for $sys -- this system's test will fail at Step 1"
  fi
done

declare -A EXPECTED_DIST=( [sys1]="1.97" [sys2]="1.88" [sys3]="2.18" [sys4]="2.24" )
declare -A EXPECTED_LIG=( [sys1]="37" [sys2]="32" [sys3]="33" [sys4]="35" )

echo
echo "=== Running all 4 reproduction scripts in the scratch copy ==="
cd "$SCRATCH_DIR"
RESULTS=()
for n in 1 2 3 4; do
  sys="sys${n}"
  script="$SOURCE_DIR/REPRODUCE_${sys}.sh"
  echo
  echo "--- $sys ---"
  if [ ! -f "$script" ]; then
    echo "MISSING SCRIPT: $script"
    RESULTS+=("$sys: FAIL (script not found)")
    continue
  fi
  LOG="$SCRATCH_DIR/${sys}_run.log"
  if bash "$script" > "$LOG" 2>&1; then
    DIST=$(grep -oP 'min_dist=\K[0-9.]+' "$LOG" | tail -1)
    LIG=$(grep -oP 'ligand=\K[0-9]+' "$LOG" | tail -1)
    if [ "$DIST" = "${EXPECTED_DIST[$sys]}" ] && [ "$LIG" = "${EXPECTED_LIG[$sys]}" ]; then
      RESULTS+=("$sys: PASS (ligand=$LIG, min_dist=$DIST A, matches confirmed reference)")
    else
      RESULTS+=("$sys: RAN BUT MISMATCH (got ligand=$LIG min_dist=$DIST, expected ligand=${EXPECTED_LIG[$sys]} min_dist=${EXPECTED_DIST[$sys]}) -- see $LOG")
    fi
  else
    RESULTS+=("$sys: FAIL (script exited with error) -- see $LOG")
  fi
done

echo
echo "=== BATCH REPRODUCTION SUMMARY ==="
for r in "${RESULTS[@]}"; do
  echo "$r"
done
echo
echo "Scratch directory preserved for inspection: $SCRATCH_DIR"
echo "(delete manually once reviewed -- this script never touches your real working files)"
