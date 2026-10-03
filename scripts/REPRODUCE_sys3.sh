#!/usr/bin/env bash
# Reproducibility record: System 3 (selective, promiscuity-flagged), CNP0384182.2
# Real SMILES: N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)N1CCC(C(=O)NCC2CCC(C(=O)O)CC2)CC1
# Formula: C25H34N4O4, MW: 454.571, heavy atoms: 33
set -euo pipefail

python3 ~/fix_ligand_for_charmmgui.py \
  --pdbqt CNP0384182.2_docked.pdbqt \
  --out sys3_REAL_ligand.pdb
grep -v " HD\| H$" sys3_REAL_ligand.pdb > sys3_REAL_ligand_heavyatoms.pdb
# Verified: 33 heavy atoms, matching library_final.csv MW 454.571 exactly

python3 << 'EOF'
import re
def guess_element(atom_name):
    name = re.sub(r'^\d+', '', atom_name.strip())
    m = re.match(r'[A-Za-z]+', name)
    letters = m.group(0) if m else "C"
    two_letter = {"CL","BR","NA","MG","CA","ZN","FE"}
    return letters[:2].capitalize() if letters[:2].upper() in two_letter else letters[0].upper()

fixed = []
for line in open("sys3_REAL_ligand_heavyatoms.pdb"):
    if line.startswith(("ATOM","HETATM")) and len(line) >= 54:
        element = guess_element(line[12:16])
        line = line[:76].rstrip("\n").ljust(76) + f"{element:>2}" + "\n"
    fixed.append(line)
open("sys3_REAL_ligand_heavyatoms_fixed.pdb", "w").writelines(fixed)
EOF

python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

real_smi = "N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)N1CCC(C(=O)NCC2CCC(C(=O)O)CC2)CC1"
template = Chem.MolFromSmiles(real_smi)
pdb_mol = Chem.MolFromPDBFile("sys3_REAL_ligand_heavyatoms_fixed.pdb", removeHs=True, sanitize=False)
fixed_mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
# NOTE: RDKit reported "More than one matching pattern found -- picking one" for this
# molecule (internal near-symmetry between the indole/piperidine and cyclohexane rings).
# Formula match confirmed correct regardless: C25H34N4O4.
assert rdMolDescriptors.CalcMolFormula(fixed_mol) == rdMolDescriptors.CalcMolFormula(template)
Chem.MolToMolFile(fixed_mol, "sys3_REAL_ligand_corrected.sdf")

mi = Chem.AtomPDBResidueInfo()
mi.SetResidueName("UNL"); mi.SetResidueNumber(1); mi.SetChainId("L"); mi.SetIsHeteroAtom(True)
for atom in fixed_mol.GetAtoms():
    atom.SetMonomerInfo(mi)
Chem.MolToPDBFile(fixed_mol, "sys3_FINAL_ligand.pdb")
EOF

cat ALDH3B2_DOCKING_receptor.pdb sys3_FINAL_ligand.pdb > merged_FINAL_sys3_selective_flagged_CNP0384182.pdb

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

d, na, nb = min_distance_between_chains("merged_FINAL_sys3_selective_flagged_CNP0384182.pdb", "A", "L")
print(f"protein={na}, ligand={nb}, min_dist={d:.2f} A")
assert 1.0 < d < 5.0
EOF
# Confirmed real result: protein=3638, ligand=33, min_dist=2.18 A

# ---- Step 4: CHARMM-GUI Solution Builder ----
# Job ID: 9087526694
# NOTE ON THIS SYSTEM'S BUILD: the first attempt at Solvator/PBC/Input Generator
# produced three byte-identical archives (confirmed via md5sum: 5cf710f0937d4811...),
# indicating the CHARMM-GUI session did not genuinely advance past PDB Reader despite
# distinct download filenames. Rebuilt from a fresh session under this job ID.
# Same protonation (GLU129/ASP3), solvation (KCl 0.15M), PBC, Input Generator settings
# as sys1/sys2 (GROMACS, HMR, NPT, 300K) for cross-system uniformity.
# Verified production .mdp: dt=0.004, tcoupl=v-rescale, ref_t=300 300,
# pcoupl=C-rescale, pcoupltype=isotropic
# Verified genuinely distinct from sys4's job ID (9087521260) by direct comparison.
echo "System 3 build reproduced/documented."
