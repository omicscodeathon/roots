#!/usr/bin/env python3
"""
Convert a REAL docked PDBQT pose into a CHARMM-GUI-ready PDB file: correct chain ID (L),
correct residue name (UNL), real docked coordinates preserved exactly (NOT re-embedded).

This fixes the root cause confirmed tonight: prep_charmmgui_inputs_v2.py's RDKit embedding
was only ever meant to produce a CGenFF-format-testable ligand, its coordinates (near-origin,
arbitrary from SMILES embedding) were never meant to be merged with the real protein structure.
The actual binding pose lives in the Vina/QuickVina2 docked PDBQT files
(step6_results/roots-docking-*/results_*.tar.gz), and THOSE coordinates -- not RDKit's fresh
embedding -- are what must be merged with the protein for a physically meaningful MD system.

Confirmed real symptom of using the wrong file: a 0.44 A minimum distance between ligand and
protein atoms (effectively coincident atoms), producing a single-point energy of ~9.8 billion
in CHARMM-GUI's PDB Reader step -- this script's output should instead produce close but
non-overlapping contacts (confirmed earlier: 1.4-2.9 A ligand-centroid-to-NAD-centroid distance,
consistent with a real, non-clashing bound pose once NAD itself is excluded).

Usage:
    python3 fix_ligand_for_charmmgui.py --pdbqt CNP0220660.1_docked.pdbqt --out sys1_real_ligand.pdb
"""
import argparse
import sys
from pathlib import Path

CHAIN_ID = "L"
RESNAME = "UNL"
RESNUM = 1


def convert(pdbqt_path, out_path, chain_id=CHAIN_ID, resname=RESNAME, resnum=RESNUM):
    """Read only the FIRST model's ATOM/HETATM lines (Vina's best-ranked pose, confirmed
    convention from earlier tonight's work), rewrite as HETATM with correct chain/residue
    fields, preserving the real x/y/z coordinates exactly -- no re-embedding, no re-centering."""
    lines_out = []
    in_first_model = True
    n_atoms = 0

    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54:
                print(f"WARNING: skipping malformed line (too short): {line!r}", file=sys.stderr)
                continue

            # Confirmed real column layout from an actual docked PDBQT line (verified by direct
            # inspection): columns 0-16 = record type/serial/atom name (preserved as-is),
            # columns 17-19 = residue name (3 chars), column 20 = blank, column 21 = chain ID
            # (blank in the source -- this is what we fill in), columns 22-25 = residue seq num
            # (right-justified, 4 chars), column 26 = blank, column 27 onward = coordinates.
            # Rebuilding explicitly by fixed column position rather than strings-plus-slicing,
            # since the earlier attempt's off-by-one concatenation (producing "UNLL" with no
            # separating space) showed that approach is too easy to get subtly wrong.
            prefix = line[:17]                      # record type + serial + atom name
            new_resname = f"{resname:<3}"            # columns 17-19
            new_line = (
                prefix
                + new_resname                        # 17-19
                + " "                                 # 20
                + chain_id[:1]                        # 21
                + f"{resnum:>4}"                       # 22-25
                + " "                                  # 26
                + line[27:]                            # coordinates onward, REAL pose preserved
            )
            new_line = "HETATM" + new_line[6:]
            lines_out.append(new_line)
            n_atoms += 1

    if n_atoms == 0:
        sys.exit(f"ERROR: no ATOM/HETATM lines found in the first model of {pdbqt_path}. "
                 f"Check the file is a real docked PDBQT, not empty or malformed.")

    with open(out_path, "w") as f:
        f.writelines(lines_out)
        f.write("END\n")

    return n_atoms


def verify(out_path, expected_chain=CHAIN_ID):
    """Confirm every written line has the correct chain ID and that coordinates were preserved
    (not reset to near-zero, which was the original bug's symptom)."""
    coords = []
    bad_chain = 0
    for line in open(out_path):
        if line.startswith(("ATOM", "HETATM")):
            if len(line) < 22 or line[21] != expected_chain:
                bad_chain += 1
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                pass
    if not coords:
        return bad_chain, None
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return bad_chain, (cx, cy, cz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbqt", required=True, help="the REAL docked pose PDBQT file")
    ap.add_argument("--out", required=True, help="output PDB path, CGenFF/CHARMM-GUI-ready")
    ap.add_argument("--chain-id", default=CHAIN_ID)
    ap.add_argument("--resname", default=RESNAME)
    args = ap.parse_args()

    pdbqt_path = Path(args.pdbqt)
    if not pdbqt_path.exists():
        sys.exit(f"File not found: {pdbqt_path}")

    n_atoms = convert(pdbqt_path, args.out, chain_id=args.chain_id, resname=args.resname)
    print(f"Converted {n_atoms} atoms from {pdbqt_path} -> {args.out}")

    bad_chain, centroid = verify(args.out, expected_chain=args.chain_id)
    if bad_chain:
        print(f"WARNING: {bad_chain} lines have an incorrect or missing chain ID")
    else:
        print(f"Chain ID verified correct on all {n_atoms} atoms")

    if centroid:
        cx, cy, cz = centroid
        print(f"Ligand centroid in output file: ({cx:.2f}, {cy:.2f}, {cz:.2f})")
        # Sanity check: this should be near the known docking box center, NOT near the origin
        box_center = (-11.727, 4.248, 6.293)
        d = sum((a - b) ** 2 for a, b in zip(centroid, box_center)) ** 0.5
        print(f"Distance from known docking box center {box_center}: {d:.2f} A")
        if d > 10:
            print("WARNING: centroid is far from the expected docking box region -- "
                  "this does NOT look like a real docked pose. Do not proceed to merge "
                  "until this is understood.")
        else:
            print("Centroid is consistent with the real docking box region -- "
                  "this looks like a genuine docked pose, safe to proceed to merge.")


if __name__ == "__main__":
    main()
