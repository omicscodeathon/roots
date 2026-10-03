#!/usr/bin/env python3
"""
Reorder an SDF's atoms to match a reference PDB's atom order exactly, using RDKit's graph-based
substructure matching (GetSubstructMatch), not manual/visual inspection. This guarantees a correct
atom-to-atom correspondence even for fused/bridged ring systems where visual matching is error-prone,
exactly the failure mode that produced CHARMM-GUI's XX1-XX4 placeholder atoms.

Usage: python3 reorder_sdf.py --pdb reference.pdb --sdf input.sdf --out reordered.sdf
"""
import argparse
import sys
from rdkit import Chem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="the REAL docked-pose PDB (heavy atoms only) -- the reference order")
    ap.add_argument("--sdf", required=True, help="the SDF to reorder so it matches the PDB's atom order")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pdb_mol = Chem.MolFromPDBFile(args.pdb, removeHs=True, sanitize=True)
    if pdb_mol is None:
        sys.exit(f"ERROR: RDKit could not parse {args.pdb}. Check it's a valid heavy-atom-only PDB.")

    sdf_mol = next(Chem.SDMolSupplier(args.sdf, removeHs=True))
    if sdf_mol is None:
        sys.exit(f"ERROR: RDKit could not parse {args.sdf}.")

    n_pdb = pdb_mol.GetNumAtoms()
    n_sdf = sdf_mol.GetNumAtoms()
    print(f"PDB atoms: {n_pdb}")
    print(f"SDF atoms: {n_sdf}")
    if n_pdb != n_sdf:
        sys.exit(f"ERROR: atom count mismatch ({n_pdb} vs {n_sdf}) -- this is NOT an ordering issue, "
                 f"these are structurally different molecules. Do not proceed with reordering; "
                 f"re-check which SMILES/structure actually generated each file.")

    # Match by connectivity/topology (not coordinates, not names) -- this is what correctly handles
    # fused/bridged rings where a simple traversal-order heuristic fails.
    match = sdf_mol.GetSubstructMatch(pdb_mol)
    if not match:
        # try the reverse direction in case of template/query asymmetry
        match = pdb_mol.GetSubstructMatch(sdf_mol)
        if match:
            # match gives PDB-atom-index -> SDF-atom-index order; invert it
            inverse = [0] * len(match)
            for pdb_idx, sdf_idx in enumerate(match):
                inverse[sdf_idx] = pdb_idx
            match = tuple(inverse)

    if not match:
        sys.exit("ERROR: RDKit could not find a substructure match between the PDB and SDF at all.\n"
                 "This means they are likely different molecules (different connectivity/bonding),\n"
                 "not just differently ordered. Re-verify both files represent the same compound\n"
                 "before proceeding -- do NOT force a reorder on a genuine structural mismatch.")

    print(f"Substructure match found: {len(match)} atoms matched")
    if len(match) != n_pdb:
        sys.exit(f"ERROR: partial match only ({len(match)} of {n_pdb} atoms) -- the match is incomplete.\n"
                 "This usually means a real connectivity difference (e.g. different bond order,\n"
                 "different protonation/tautomer) rather than just atom-order. Investigate before\n"
                 "trusting this reorder.")

    # 'match' is: for each atom in pdb_mol (reference order), which atom index in sdf_mol corresponds.
    # RenumberAtoms needs the inverse: for each NEW position (following PDB order), which OLD sdf index.
    reordered = Chem.RenumberAtoms(sdf_mol, list(match))

    writer = Chem.SDWriter(args.out)
    writer.write(reordered)
    writer.close()
    print(f"Written: {args.out}")

    # Verify: element-by-element, does the reordered SDF now match the PDB's element sequence?
    pdb_elements = [a.GetSymbol() for a in pdb_mol.GetAtoms()]
    sdf_elements = [a.GetSymbol() for a in reordered.GetAtoms()]
    mismatches = sum(1 for p, s in zip(pdb_elements, sdf_elements) if p != s)
    print(f"\nElement-order verification: {len(pdb_elements) - mismatches}/{len(pdb_elements)} positions match")
    if mismatches == 0:
        print("PERFECT MATCH -- the reordered SDF's atom sequence now exactly mirrors the PDB's.")
    else:
        print(f"WARNING: {mismatches} position(s) still differ in element identity -- "
              f"do not trust this output, investigate further before using it.")


if __name__ == "__main__":
    main()
