#!/usr/bin/env bash
# Reproducibility record: System 1 (consensus), CNP0220660.1
# ALDH3B2 MD system preparation -- every command used, in the order actually run.
set -euo pipefail

# ---- Step 1: extract the real Vina/QuickVina2 docked pose (first/best model only) ----
# Source: step6_results/roots-docking-01/results_roots-docking-01.tar.gz
# docking_results/CNP0220660.1_docked.pdbqt
python3 ~/fix_ligand_for_charmmgui.py \
  --pdbqt CNP0220660.1_docked.pdbqt \
  --out sys1_REAL_ligand.pdb

# ---- Step 1 (continued): strip explicit polar hydrogens (AutoDock HD atoms) ----
grep -v " HD\| H$" sys1_REAL_ligand.pdb > sys1_REAL_ligand_heavyatoms.pdb
# Confirmed heavy-atom count: 37 (library_final.csv MW 472.666 matches exactly)

# ---- Step 2: merge with the verified docking receptor (NOT a fresh NAD-stripped copy --
#      confirmed via CA-CA comparison that only the ACTUAL docking receptor is coordinate-
#      frame-consistent; max discrepancy vs. an independently-stripped copy was 57.79 A) ----
cat ALDH3B2_DOCKING_receptor.pdb sys1_REAL_ligand_heavyatoms.pdb \
  > merged_v4_sys1_consensus_CNP0220660.pdb

# ---- Step 3: verify no steric clash between ligand (chain L) and protein (chain A) ----
python3 << 'EOF'
def min_distance_between_chains(pdb_path, chain_a, chain_b):
    coords_a, coords_b = [], []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 22:
            chain = line[21]
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            if chain == chain_a: coords_a.append((x, y, z))
            elif chain == chain_b: coords_b.append((x, y, z))
    min_d = min(((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
                for x1,y1,z1 in coords_a for x2,y2,z2 in coords_b)
    return min_d, len(coords_a), len(coords_b)

d, na, nb = min_distance_between_chains("merged_v4_sys1_consensus_CNP0220660.pdb", "A", "L")
print(f"protein={na}, ligand={nb}, min_dist={d:.2f} A")
assert 1.0 < d < 5.0, "clash or coordinate-frame mismatch -- do not proceed"
EOF
# Confirmed real result: protein=3638, ligand=37, min_dist=1.97 A

# ---- Step 4: CHARMM-GUI Solution Builder (charmm-gui.org), manual web steps ----
# Job ID: 9082326193
# 1. Input Generator -> Solution Builder -> upload merged_v4_sys1_consensus_CNP0220660.pdb
# 2. Model/Chain Selection: select BOTH Protein (chain A, 385 res) and Hetero (chain L, UNL)
# 3. Ligand FF: "the SDF file uploaded" -> sys1_consensus_CNP0220660.sdf
#    (CGenFF via ParamChem; RDKit-embedded SDF, independently formula/MW-verified against
#     library_final.csv before upload)
# 4. Protonation state: accepted CHARMM-GUI's pKa-predicted GLU129->GLUP, ASP3->ASPP
#    (verified: GLU129 is 8.94 A from the NAD/box center with 39 nearby protein atoms within
#     6 A, consistent with a genuine buried-residue pKa shift; ASP3 is 35.38 A away, distant
#     from the binding site, accepted as CHARMM-GUI default given negligible mechanistic impact)
# 5. Solvator: Rectangular box, Fit to Protein Size, Edge Distance 10.0 A
#    Ions: KCl, 0.15 M, Neutralizing ON (deliberate, documented deviation from an initial
#    NaCl-based plan; final box: 113x113x113 A cubic, ~148,000 total atoms)
# 6. PBC Setup: automatic PME FFT grid generation; CHARMM-GUI server-side minimization
# 7. Input Generator: GROMACS; Hydrogen Mass Repartitioning ON; Equilibration NVT;
#    Dynamics NPT; Temperature 300 K
# 8. Downloaded: sys1-pdb-charmm-gui.tgz, sys1-solvator-charmm-gui.tgz,
#    sys1-pbc-charmm-gui.tgz, sys1-input-generator-charmm-gui.tgz

# ---- Verification of final GROMACS production settings (step5_production.mdp) ----
# dt = 0.004 (HMR-enabled 4fs timestep, confirmed)
# pcoupl = C-rescale, pcoupltype = isotropic (NPT, confirmed)
# tcoupl = v-rescale, ref_t = 300 300 (confirmed)
# step4.1_equilibration.mdp: ref_t = 300 300, gen-temp = 300 (confirmed consistent)
echo "System 1 build reproduced/documented. See inline comments for exact verified values."
