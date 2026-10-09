#!/usr/bin/env bash
# 3_run_mmpbsa.sh -- MM-PBSA binding free energy calculation for ROOTS pipeline
#
# Runs ON the Hetzner CCX23 instance (after 2_hetzner_provision.sh).
# Uses gmx_MMPBSA (Valdes-Tresanco et al., JCTC 2021) with CHARMM36 topologies.
#
# For each of 4 systems:
#   1. Builds stripped trajectory (protein + ligand only, no solvent/ions)
#   2. Creates GROMACS index groups for receptor and ligand
#   3. Runs gmx_MMPBSA with:
#      - GB (igb=8, GBneck2) for fast, reasonably accurate results
#      - PB (linearized Poisson-Boltzmann) for rigorous comparison
#      - Per-residue energy decomposition (identifies key binding residues)
#      - Entropy via interaction entropy (IE) method
#   4. Saves all outputs to /root/roots/results/mmpbsa/
#
# Expected runtime: ~4-6h on CCX23 (4 dedicated vCPU) for all 4 systems.
set -euo pipefail

source /opt/roots_env/bin/activate

BASE="/root/roots"
RESULTS="$BASE/results/mmpbsa"
mkdir -p "$RESULTS"

# ---- configuration ----
# Frame selection: use the equilibrated portion of the trajectory.
# START_NS: skip initial equilibration; END_NS: 0 = use all remaining frames.
# DT_PS: extract one frame every DT_PS picoseconds (100 ps = 1 frame/0.1ns).
# With ~30ns of usable trajectory at dt=100ps, this gives ~300 frames -- adequate
# for MM-PBSA convergence (literature: 100-500 frames typical).
START_NS="${START_NS:-10}"    # discard first 10 ns (equilibration)
END_NS="${END_NS:-0}"         # 0 = to end of trajectory
DT_PS="${DT_PS:-100}"         # one frame every 100 ps

LIGAND_RES="UNL"              # ligand residue name from CHARMM-GUI
NPROC=$(nproc)

echo "=== ROOTS MM-PBSA analysis ==="
echo "Start: $(date)"
echo "CPUs: $NPROC"
echo "Frame selection: start=${START_NS}ns, end=${END_NS}ns (0=all), dt=${DT_PS}ps"
echo

for sys in Sys1 Sys2 Sys3 Sys4; do
  echo "====== Processing $sys ======"
  SYSDIR="$BASE/harvest/$sys"
  TOPODIR="$BASE/topo/$sys"
  OUTDIR="$RESULTS/$sys"
  mkdir -p "$OUTDIR"

  # Find input files
  TPR=$(find "$SYSDIR" -name '*.tpr' | head -1)
  XTC=$(find "$SYSDIR" -name '*.xtc' | sort | tail -1)  # latest part
  NDX="$TOPODIR/index.ndx"
  TOP="$TOPODIR/topol.top"

  if [[ -z "$TPR" || -z "$XTC" ]]; then
    echo "  WARNING: missing .tpr or .xtc in $SYSDIR -- skipping."
    continue
  fi
  echo "  TPR: $TPR"
  echo "  XTC: $XTC"
  echo "  TOP: $TOP"

  # ---- Step 1: Determine trajectory bounds ----
  echo "  Checking trajectory length..."
  CHECK_OUT=$(gmx check -f "$XTC" 2>&1 || true)
  echo "$CHECK_OUT" > "$OUTDIR/gmx_check.log"

  # Extract last time from gmx check output
  LAST_TIME=$(echo "$CHECK_OUT" | grep -oP 'Last frame\s+\d+\s+time\s+\K[0-9.]+' || \
              echo "$CHECK_OUT" | grep -oP 'time\s+\K[0-9.]+' | sort -n | tail -1 || \
              echo "unknown")
  echo "  Trajectory last time: ${LAST_TIME} ps"

  # Convert start/end to ps for gmx trjconv
  START_PS=$((START_NS * 1000))
  if [[ "$END_NS" -eq 0 ]]; then
    END_FLAG=""
  else
    END_PS=$((END_NS * 1000))
    END_FLAG="-e $END_PS"
  fi

  # ---- Step 2: Make index groups ----
  # Create custom index with Receptor (Protein) and Ligand (UNL) groups
  echo "  Creating index groups..."
  cat > "$OUTDIR/make_ndx.input" <<NDXEOF
q
NDXEOF

  # Use gmx make_ndx to identify group numbers, then build receptor+ligand complex
  # First, get available groups
  GROUPS_OUT=$(echo "q" | gmx make_ndx -f "$TPR" -n "$NDX" -o "$OUTDIR/temp.ndx" 2>&1 || true)
  echo "$GROUPS_OUT" > "$OUTDIR/groups.log"

  # Build receptor-ligand index for gmx_MMPBSA
  # gmx_MMPBSA needs: complex group (protein+ligand), receptor group, ligand group
  cat > "$OUTDIR/ndx_complex.input" <<NDXEOF
"Protein" | r $LIGAND_RES
name 0 Complex
"Protein"
name 1 Receptor
r $LIGAND_RES
name 2 Ligand
q
NDXEOF

  # Note: gmx select syntax may vary; using make_ndx for robustness
  echo '1 | r UNL' > "$OUTDIR/mk1.txt"
  echo 'q' >> "$OUTDIR/mk1.txt"

  # Create the complex index
  gmx make_ndx -f "$TPR" -o "$OUTDIR/complex.ndx" <<MKNDX 2>&1 | tee "$OUTDIR/make_ndx.log"
"Protein" | r ${LIGAND_RES}
q
MKNDX

  # ---- Step 3: Extract complex-only trajectory (strip solvent/ions) ----
  echo "  Extracting complex trajectory (protein + ligand)..."
  # Identify the group number for Protein_UNL (the combined group make_ndx created)
  # It's typically the last group created. We'll use the group name.
  # For gmx trjconv, pipe the group selection
  echo "Protein_${LIGAND_RES}" | gmx trjconv \
    -f "$XTC" \
    -s "$TPR" \
    -n "$OUTDIR/complex.ndx" \
    -o "$OUTDIR/complex_traj.xtc" \
    -b "$START_PS" \
    $END_FLAG \
    -dt "$DT_PS" \
    -pbc mol \
    -center <<TRJEOF 2>&1 | tee "$OUTDIR/trjconv.log"
Protein
Protein_${LIGAND_RES}
TRJEOF

  # Also extract a reference structure (first frame)
  echo "Protein_${LIGAND_RES}" | gmx trjconv \
    -f "$XTC" \
    -s "$TPR" \
    -n "$OUTDIR/complex.ndx" \
    -o "$OUTDIR/complex_ref.gro" \
    -b "$START_PS" \
    -e "$START_PS" \
    -pbc mol \
    -center <<TRJEOF2 2>&1 | tee -a "$OUTDIR/trjconv.log"
Protein
Protein_${LIGAND_RES}
TRJEOF2

  # ---- Step 4: Create gmx_MMPBSA input files ----
  # GB calculation (fast, ~1h per system)
  cat > "$OUTDIR/mmpbsa_gb.in" <<'GBIN'
&general
  sys_name="ROOTS",
  startframe=1,
  endframe=9999999,
  interval=1,
  forcefields="charmm",
  PBRadii=3,
  temperature=310,
/
&gb
  igb=8, saltcon=0.150,
  ifqnt=0,
/
&decomp
  idecomp=1, dec_verbose=1,
  print_res="within 6"
/
GBIN

  # PB calculation (rigorous, ~3-4h per system)
  cat > "$OUTDIR/mmpbsa_pb.in" <<'PBIN'
&general
  sys_name="ROOTS",
  startframe=1,
  endframe=9999999,
  interval=1,
  forcefields="charmm",
  PBRadii=3,
  temperature=310,
/
&pb
  istrng=0.150, fillratio=4.0,
  inp=1, radiopt=0,
  indi=1.0, exdi=80.0,
/
&decomp
  idecomp=1, dec_verbose=1,
  print_res="within 6"
/
PBIN

  # ---- Step 5: Run gmx_MMPBSA ----
  echo "  Running gmx_MMPBSA (GB model, igb=8)..."
  echo "  This typically takes 1-2h per system on 4 CPUs..."

  cd "$OUTDIR"

  # gmx_MMPBSA command
  # -i: input file, -cs: complex structure, -ci: complex index,
  # -cg: complex group receptor group ligand group (0-indexed from the complex.ndx we made),
  # -ct: complex trajectory, -cp: complex topology
  # -nogui: don't launch the GUI analyzer
  gmx_MMPBSA -O \
    -i mmpbsa_gb.in \
    -cs "$TPR" \
    -ci complex.ndx \
    -cg "Protein" "${LIGAND_RES}" \
    -ct complex_traj.xtc \
    -cp "$TOP" \
    -o "$OUTDIR/FINAL_RESULTS_MMPBSA_GB.dat" \
    -eo "$OUTDIR/FINAL_RESULTS_MMPBSA_GB.csv" \
    -do "$OUTDIR/FINAL_DECOMP_MMPBSA_GB.dat" \
    -deo "$OUTDIR/FINAL_DECOMP_MMPBSA_GB.csv" \
    -nogui \
    2>&1 | tee "$OUTDIR/mmpbsa_gb_run.log"

  echo "  GB calculation complete for $sys."

  echo "  Running gmx_MMPBSA (PB model)..."
  gmx_MMPBSA -O \
    -i mmpbsa_pb.in \
    -cs "$TPR" \
    -ci complex.ndx \
    -cg "Protein" "${LIGAND_RES}" \
    -ct complex_traj.xtc \
    -cp "$TOP" \
    -o "$OUTDIR/FINAL_RESULTS_MMPBSA_PB.dat" \
    -eo "$OUTDIR/FINAL_RESULTS_MMPBSA_PB.csv" \
    -do "$OUTDIR/FINAL_DECOMP_MMPBSA_PB.dat" \
    -deo "$OUTDIR/FINAL_DECOMP_MMPBSA_PB.csv" \
    -nogui \
    2>&1 | tee "$OUTDIR/mmpbsa_pb_run.log"

  echo "  PB calculation complete for $sys."
  echo
done

echo "=== MM-PBSA analysis complete for all systems ==="
echo "End: $(date)"
echo
echo "Results in: $RESULTS/"
echo "Key output files per system:"
echo "  FINAL_RESULTS_MMPBSA_GB.dat  -- GB binding free energy summary"
echo "  FINAL_RESULTS_MMPBSA_PB.dat  -- PB binding free energy summary"
echo "  FINAL_DECOMP_MMPBSA_*.csv    -- per-residue decomposition (for binding hotspot analysis)"
echo "  *_run.log                    -- full run logs"
echo
echo "Next: run 4_run_analysis.py for trajectory analysis and visualization."
