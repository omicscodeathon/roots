#!/usr/bin/env bash
# Reproducibility record: System 2 (consensus), CNP0358193.1
# Real SMILES: CC1(C)C(=O)C=C[C@]2(C)[C@@H]3CC[C@H]4[C@H](C5=COC=C5)OC(=O)[C@H]5O[C@]54[C@]3(C)C(=O)C[C@]12O
# Formula: C25H28O7, MW: 440.492, heavy atoms: 32
set -euo pipefail

# ---- Step 1: extract real docked pose ----
# Source: step6_results/roots-docking-03/results_roots-docking-03.tar.gz
python3 ~/fix_ligand_for_charmmgui.py \
  --pdbqt CNP0358193.1_docked.pdbqt \
  --out sys2_REAL_ligand.pdb
grep -v " HD\| H$" sys2_REAL_ligand.pdb > sys2_REAL_ligand_heavyatoms.pdb
# Verified: 32 heavy atoms, matching library_final.csv MW 440.492 exactly

# ---- Step 1c: fix the AutoDock-contaminated element column (col 77-78 held "A"/"OA"/"HD"
#      atom-type codes instead of real element symbols -- rejected by RDKit's strict PDB
#      parser with "Post-condition Violation: Element 'A' not found") ----
python3 << 'EOF'
import re
def guess_element(atom_name):
    name = re.sub(r'^\d+', '', atom_name.strip())
    m = re.match(r'[A-Za-z]+', name)
    letters = m.group(0) if m else "C"
    two_letter = {"CL","BR","NA","MG","CA","ZN","FE"}
    return letters[:2].capitalize() if letters[:2].upper() in two_letter else letters[0].upper()

fixed = []
for line in open("sys2_REAL_ligand_heavyatoms.pdb"):
    if line.startswith(("ATOM","HETATM")) and len(line) >= 54:
        element = guess_element(line[12:16])
        line = line[:76].rstrip("\n").ljust(76) + f"{element:>2}" + "\n"
    fixed.append(line)
open("sys2_REAL_ligand_heavyatoms_fixed.pdb", "w").writelines(fixed)
EOF

# ---- Step 1d: bond-order correction via template (CRITICAL step for this compound) ----
# Root cause: the independently RDKit-embedded SDF (sys2_consensus_CNP0358193.sdf) was
# structurally INCOMPLETE for this furan-fused pentacyclic scaffold -- confirmed via
# formula check: C21H28*4O7 (4 unidentified "*" dummy atoms, only 31 bonds/5 rings vs.
# the real 37 bonds/6 rings). Separately, naive PDB->RDKit bond perception (distance-based,
# no explicit bond orders) on the real docked pose over-saturated double bonds: C25H40O7
# instead of the correct C25H28O7 (12 extra H, consistent with 6 misassigned double bonds).
# Fix: use the VERIFIED library SMILES as a bond-order template applied to the REAL docked
# 3D coordinates -- correct chemistry + correct geometry, simultaneously.
python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

real_smi = "CC1(C)C(=O)C=C[C@]2(C)[C@@H]3CC[C@H]4[C@H](C5=COC=C5)OC(=O)[C@H]5O[C@]54[C@]3(C)C(=O)C[C@]12O"
template = Chem.MolFromSmiles(real_smi)
pdb_mol = Chem.MolFromPDBFile("sys2_REAL_ligand_heavyatoms_fixed.pdb", removeHs=True, sanitize=False)
fixed_mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
assert rdMolDescriptors.CalcMolFormula(fixed_mol) == rdMolDescriptors.CalcMolFormula(template)
Chem.MolToMolFile(fixed_mol, "sys2_REAL_ligand_corrected.sdf")

mi = Chem.AtomPDBResidueInfo()
mi.SetResidueName("UNL"); mi.SetResidueNumber(1); mi.SetChainId("L"); mi.SetIsHeteroAtom(True)
for atom in fixed_mol.GetAtoms():
    atom.SetMonomerInfo(mi)
Chem.MolToPDBFile(fixed_mol, "sys2_FINAL_ligand.pdb")
print("Formula-verified: C25H28O7, atom order IDENTICAL between PDB and SDF by construction")
EOF

# ---- Step 2: merge ----
cat ALDH3B2_DOCKING_receptor.pdb sys2_FINAL_ligand.pdb > merged_FINAL_sys2_consensus_CNP0358193.pdb

# ---- Step 3: clash verification ----
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

d, na, nb = min_distance_between_chains("merged_FINAL_sys2_consensus_CNP0358193.pdb", "A", "L")
print(f"protein={na}, ligand={nb}, min_dist={d:.2f} A")
assert 1.0 < d < 5.0
EOF
# Confirmed real result: protein=3638, ligand=32, min_dist=1.88 A

# ---- Step 4: CHARMM-GUI Solution Builder ----
# Job ID: 9086874719
# SDF used: sys2_REAL_ligand_corrected.sdf (bond-order-corrected, atom-order-matched)
# Result: NO atom-order mismatch table (unlike the first attempt with the broken SDF)
# Same protonation (GLU129/ASP3), solvation (KCl 0.15M), PBC, and Input Generator
# settings as System 1 (GROMACS, HMR, NPT, 300K) for cross-system uniformity.
# Verified production .mdp: dt=0.004, pcoupl=C-rescale, tcoupl=v-rescale, ref_t=300 300
# Verified equilibration .mdp: ref_t=300 300, gen-temp=300 (consistent with production)
echo "System 2 build reproduced/documented."
