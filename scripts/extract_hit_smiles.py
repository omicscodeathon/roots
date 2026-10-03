#!/usr/bin/env python3
import argparse, csv, sys

def main():
    p = argparse.ArgumentParser(description="Attach SMILES to hit compound_ids from a manifest.")
    p.add_argument("--hits", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--manifest-id-col", required=True)
    p.add_argument("--manifest-smiles-col", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    manifest_map = {}
    with open(args.manifest, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if args.manifest_id_col not in reader.fieldnames:
            sys.exit(f"ERROR: '{args.manifest_id_col}' not in manifest columns: {reader.fieldnames}")
        if args.manifest_smiles_col not in reader.fieldnames:
            sys.exit(f"ERROR: '{args.manifest_smiles_col}' not in manifest columns: {reader.fieldnames}")
        for row in reader:
            manifest_map[row[args.manifest_id_col]] = row[args.manifest_smiles_col]

    with open(args.hits, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "compound_id" not in reader.fieldnames:
            sys.exit(f"ERROR: 'compound_id' not in hits columns: {reader.fieldnames}")
        fieldnames = reader.fieldnames + ["smiles"]
        rows = list(reader)

    matched, unmatched = [], []
    for row in rows:
        smi = manifest_map.get(row["compound_id"])
        if smi is None:
            unmatched.append(row["compound_id"])
        else:
            row["smiles"] = smi
            matched.append(row)

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(matched)

    unmatched_path = args.out.rsplit(".", 1)[0] + "_unmatched.txt"
    with open(unmatched_path, "w", encoding="utf-8") as f:
        for cid in unmatched:
            f.write(cid + "\n")

    print(f"Matched: {len(matched)}  Unmatched: {len(unmatched)}")
    print(f"Wrote: {args.out}")
    print(f"Wrote: {unmatched_path}")

if __name__ == "__main__":
    main()
