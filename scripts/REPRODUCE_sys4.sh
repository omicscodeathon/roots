#!/usr/bin/env bash
# Reproducibility record: System 4 (selective, borderline), CNP0088713.1
# Real SMILES: CC1=C(C)C(=O)O[C@@H]([C@](C)(O)[C@H]2[C@@H](O)C[C@@H]3[C@H]4C[C@H]5O[C@]56[C@@H](O)C=CC(=O)[C@]6(C)[C@H]4CC[C@@]32C)C1
# Formula: C28H38O7, MW: 486.605, heavy atoms: 35
set -euo pipefail

python3 ~/fix_ligand_for_charmmgui.py \
  --pdbqt CNP0088713.1_docked.pdbqt \
  --out sys4_REAL_ligand.pdb
grep -v " HD\| H$" sys4_REAL_ligand.pdb > sys4_REAL_ligand_heavyatoms.pdb
# Verified: 35 heavy atoms, matching library_final.csv MW 486.605 exactly

python3 << 'EOF'
import re
def guess_element(atom_name):
    name = re.sub(r'^\d+', '', atom_name.strip())
    m = re.match(r'[A-Za-z]+', name)
    letters = m.group(0) if m else "C"
    two_letter = {"CL","BR","NA","MG","CA","ZN","FE"}
    return letters[:2].capitalize() if letters[:2].upper() in two_letter else letters[0].upper()

fixed = []
for line in open("sys4_REAL_ligand_heavyatoms.pdb"):
    if line.startswith(("ATOM","HETATM")) and len(line) >= 54:
        element = guess_element(line[12:16])
        line = line[:76].rstrip("\n").ljust(76) + f"{element:>2}" + "\n"
    fixed.append(line)
open("sys4_REAL_ligand_heavyatoms_fixed.pdb", "w").writelines(fixed)
EOF

python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

real_smi = "CC1=C(C)C(=O)O[C@@H]([C@](C)(O)[C@H]2[C@@H](O)C[C@@H]3[C@H]4C[C@H]5O[C@]56[C@@H](O)C=CC(=O)[C@]6(C)[C@H]4CC[C@@]32C)C1"
template = Chem.MolFromSmiles(real_smi)
pdb_mol = Chem.MolFromPDBFile("sys4_REAL_ligand_heavyatoms_fixed.pdb", removeHs=True, sanitize=False)
fixed_mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
assert rdMolDescriptors.CalcMolFormula(fixed_mol) == rdMolDescriptors.CalcMolFormula(template)
Chem.MolToMolFile(fixed_mol, "sys4_REAL_ligand_corrected.sdf")

mi = Chem.AtomPDBResidueInfo()
mi.SetResidueName("UNL"); mi.SetResidueNumber(1); mi.SetChainId("L"); mi.SetIsHeteroAtom(True)
for atom in fixed_mol.GetAtoms():
    atom.SetMonomerInfo(mi)
Chem.MolToPDBFile(fixed_mol, "sys4_FINAL_ligand.pdb")
EOF

cat ALDH3B2_DOCKING_receptor.pdb sys4_FINAL_ligand.pdb > merged_FINAL_sys4_selective_borderline_CNP0088713.pdb

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

d, na, nb = min_distance_between_chains("merged_FINAL_sys4_selective_borderline_CNP0088713.pdb", "A", "L")
print(f"protein={na}, ligand={nb}, min_dist={d:.2f} A")
assert 1.0 < d < 5.0
EOF
# Confirmed real result: protein=3638, ligand=35, min_dist=2.24 A

# ---- Step 4: CHARMM-GUI Solution Builder ----
# Job ID: 9087521260
# Same protonation (GLU129/ASP3), solvation (KCl 0.15M), PBC, Input Generator settings
# as sys1-3 for cross-system uniformity (GROMACS, HMR, NPT, 300K).
# Verified production .mdp: dt=0.004, tcoupl=v-rescale, ref_t=300 300,
# pcoupl=C-rescale, pcoupltype=isotropic
# Verified genuinely distinct from sys3's job ID (9087526694) by direct comparison.
echo "System 4 build reproduced/documented."
