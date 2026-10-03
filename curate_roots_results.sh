#!/usr/bin/env bash
# curate_roots_results.sh -- run this LOCALLY on your Arch machine, not on
# any pod. It copies only small, final, publishable artifacts from
# ~/roots_pipeline1 into ~/roots (your actual git repo), leaving every
# heavy, raw, or duplicate file behind. Nothing here touches a GPU pod or
# costs any compute money -- this is pure local file organization.
#
# EXCLUDED BY DESIGN (confirmed from the real structure dump):
#   - venv/                              (328k-file Python virtualenv, never belongs in git)
#   - prodigy_results/, docking results/ (raw Hetzner fleet tar.gz bundles, 5 nodes each)
#   - step4/ligands/, step4/chunks/      (280k+ individual docked PDBQT files)
#   - ALDH3B2_af3_raw/                   (raw AlphaFold3 job output: MSAs, templates, 5 model variants)
#   - *.backup_*, *.before_*, *.pre_*    (superseded script versions)
#   - __pycache__/                       (compiled bytecode)
#   - OldSys2Dump/, roots-md-staging*/, Sys_outputs/, REPRO_TEST_*/  (duplicate/superseded staging copies)
#   - charmmgui_ligand_inputs/*.tgz, *.pdbqt (large structure/archive files -- scripts that GENERATE them are kept)
set -euo pipefail

SRC="$HOME/roots_pipeline1"
DST="$HOME/roots"

[[ -d "$SRC" ]] || { echo "SRC not found: $SRC -- adjust the path at the top of this script." >&2; exit 1; }
[[ -d "$DST" ]] || { echo "DST not found: $DST -- cd to where you cloned the repo first, or adjust DST." >&2; exit 1; }

echo "=== Creating target directory structure under $DST ==="
mkdir -p "$DST/data/step1" "$DST/data/step2" "$DST/data/step3" \
         "$DST/data/step5" "$DST/data/step6" "$DST/data/step7" \
         "$DST/scripts" "$DST/docs"

echo
echo "=== Step 1: differential expression / target shortlist (corrected versions) ==="
# The CORRECTED set is the real final one (tool-identity and comparator fixes
# applied) -- copying those, not the earlier uncorrected step1_* files.
cp -v "$SRC/ROOTS_Step1B_corrected_outputs/final_shortlist_CORRECTED.csv" "$DST/data/step1/" 2>/dev/null || true
cp -v "$SRC/ROOTS_Step1B_corrected_outputs/step1_target_shortlist_CORRECTED.csv" "$DST/data/step1/" 2>/dev/null || true
cp -v "$SRC/ROOTS_Step1B_corrected_outputs/step1b_provenance_CORRECTED.json" "$DST/data/step1/" 2>/dev/null || true
cp -v "$SRC/ROOTS_Step1B_corrected_outputs/hpa_ovary_ihc_results.csv" "$DST/data/step1/" 2>/dev/null || true
cp -v "$SRC/ROOTS_Step1B_corrected_outputs/geo_validated.csv" "$DST/data/step1/" 2>/dev/null || true

echo
echo "=== Step 2: validated receptor structures (small CSV/metadata only, not raw AF3 dumps) ==="
cp -v "$SRC/ROOTS_relevant_files/Step2/step2_validated_structures.csv" "$DST/data/step2/" 2>/dev/null || true
cp -v "$SRC/ROOTS_relevant_files/Step2/logs/step2_provenance.json" "$DST/data/step2/" 2>/dev/null || true

echo
echo "=== Step 3: final compound library ==="
cp -v "$SRC/ROOTS_Step3_deliverables/library_final.csv" "$DST/data/step3/" 2>/dev/null || true

echo
echo "=== Step 5: PRODIGY-LIG cross-validation (merged scores + the script, not the 5 raw fleet bundles) ==="
cp -v "$SRC/prodigy_scores_merged.csv" "$DST/data/step5/" 2>/dev/null || true
cp -v "$SRC/prodigy_cross_validation.py" "$DST/scripts/" 2>/dev/null || true
cp -v "$SRC/merge_prodigy.py" "$DST/scripts/" 2>/dev/null || true

echo
echo "=== Step 6: docking summary (extracted scores + failure log, not 280k raw PDBQTs) ==="
# vina_scores_extracted.csv and failures.csv are the real summary artifacts;
# the per-ligand PDBQT files and chunk tarballs are intentionally excluded.
find "$SRC" -maxdepth 2 -iname 'vina_scores_extracted.csv' -exec cp -v {} "$DST/data/step6/" \; 2>/dev/null || true
cp -v "$SRC/failures.csv" "$DST/data/step6/" 2>/dev/null || true
cp -v "$SRC/step6_dock_worker.sh" "$DST/scripts/" 2>/dev/null || true 2>/dev/null || true

echo
echo "=== Step 7: ADMET triage + hit definition (small CSVs, scripts, NOT the REPRO_TEST dump) ==="
cp -v "$SRC/step7_admet/admet_triage_full.csv" "$DST/data/step7/" 2>/dev/null || true
cp -v "$SRC/step7_admet/tier2_shortlist.csv" "$DST/data/step7/" 2>/dev/null || true
cp -v "$SRC/step7_admet/hit_smiles.csv" "$DST/data/step7/" 2>/dev/null || true
cp -v "$SRC/step7_admet/swissadme.csv" "$DST/data/step7/" 2>/dev/null || true
cp -v "$SRC/step7_admet/admet_triage.py" "$DST/scripts/" 2>/dev/null || true
cp -v "$SRC/step7_admet/extract_hit_smiles.py" "$DST/scripts/" 2>/dev/null || true
cp -v "$SRC/step7_admet/pre_admet_diagnostics.py" "$DST/scripts/" 2>/dev/null || true

echo
echo "=== MD system-prep scripts (the REPRODUCE_*.sh + their helpers -- small, text, fully reproducible) ==="
for N in 1 2 3 4; do
  cp -v "$SRC/step7_admet/charmmgui_ligand_inputs/REPRODUCE_sys${N}.sh" "$DST/scripts/" 2>/dev/null || true
done
cp -v "$SRC/step7_admet/charmmgui_ligand_inputs/reorder_sdf.py" "$DST/scripts/" 2>/dev/null || true
cp -v "$SRC/step7_admet/fix_ligand_for_charmmgui.py" "$DST/scripts/" 2>/dev/null || true
cp -v "$SRC/step7_admet/charmmgui_ligand_inputs/run_all_reproductions.sh" "$DST/scripts/" 2>/dev/null || true

echo
echo "=== Deployment / infra scripts (small, text -- real, reusable for a reviewer or re-run) ==="
cp -v "$SRC/hetzner_provision.sh" "$DST/scripts/" 2>/dev/null || true
cp -v "$SRC/prodigy_hetzner_fleet.py" "$DST/scripts/" 2>/dev/null || true
# Deliberately NOT copying the 3 .backup_*/.before_*/.pre_* variants of this file.

echo
echo "=== EXPLICITLY EXCLUDED (left behind in $SRC, not copied) ==="
cat <<'EXCLUDED'
  venv/                                  (Python virtualenv)
  prodigy_results/roots-prodigy-0{0..4}/ (raw per-node PRODIGY tar.gz bundles)
  step4/ligands/, step4/chunks/          (280k+ individual docked PDBQT files)
  ALDH3B2_af3_raw/                       (raw AlphaFold3 job: MSAs, templates, 5 model .cif variants)
  step7_admet/REPRO_TEST_20261001_190150/ (superseded reproducibility test dump)
  OldSys2Dump/, roots-md-staging/, roots-md-staging2/, Sys_outputs/ (duplicate staging copies)
  __pycache__/, *.backup_*, *.before_*, *.pre_*  (build artifacts, superseded script versions)
  All raw CHARMM-GUI .tgz archives and merged_*.pdb structure files
  (these are LARGE and REGENERATABLE by the REPRODUCE_sys*.sh scripts already copied above)
EXCLUDED

echo
echo "=== Done. Review what landed in $DST before committing: ==="
find "$DST/data" "$DST/scripts" -type f -newer /tmp/curate_roots_results.sh 2>/dev/null | sort
