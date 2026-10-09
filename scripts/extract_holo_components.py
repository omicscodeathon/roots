#!/usr/bin/env python3
"""Extract NAD+ (chain B) and protein-only (chain A) from the ALDH3B2 holo CIF,
using manual line parsing since this OpenFold3 CIF lacks the occupancy field
Biopython's MMCIFParser expects."""

import sys

CIF_PATH = "ALDH3B2_openfold3_holo_NAD_seed_2746317213_sample_5_model.cif"
PROTEIN_OUT = "step5_docking/receptors/ALDH3B2_protein_only.pdb"
NAD_OUT = "step5_docking/ligands_ref/NAD_reference.pdb"

with open(CIF_PATH) as f:
    lines = f.readlines()

protein_lines = []
nad_lines = []
atom_serial = 1

for line in lines:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        continue
    parts = line.split()
    if len(parts) < 18:
        continue

    atom_name = parts[2]      # label_atom_id
    comp_id = parts[4]        # label_comp_id
    chain_id = parts[5]       # label_asym_id
    x, y, z = float(parts[15]), float(parts[16]), float(parts[17])
    element = parts[1]        # type_symbol
    auth_seq = parts[8]       # auth_seq_id

    pdb_line = (
        f"ATOM  {atom_serial:>5} {atom_name:<4} {comp_id:<3} {chain_id}"
        f"{auth_seq:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{1.00:>6.2f}{0.00:>6.2f}          {element:>2}\n"
    )

    if comp_id == "NAD":
        nad_lines.append(pdb_line)
    else:
        protein_lines.append(pdb_line)

    atom_serial += 1

with open(PROTEIN_OUT, "w") as f:
    f.writelines(protein_lines)
    f.write("END\n")

with open(NAD_OUT, "w") as f:
    f.writelines(nad_lines)
    f.write("END\n")

print(f"Protein atoms written: {len(protein_lines)} -> {PROTEIN_OUT}")
print(f"NAD+ atoms written: {len(nad_lines)} -> {NAD_OUT}")

if len(nad_lines) != 44:
    print(f"WARNING: expected 44 NAD+ atoms, got {len(nad_lines)} — check extraction.")
    sys.exit(1)
