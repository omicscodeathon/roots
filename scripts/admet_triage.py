#!/usr/bin/env python3
import argparse, csv, sys

def main():
    p = argparse.ArgumentParser(description="Tier-1 ADMET gate + rank triage.")
    p.add_argument("--admet", required=True)
    p.add_argument("--herg-col", required=True)
    p.add_argument("--ames-col", required=True)
    p.add_argument("--dili-col", required=True)
    p.add_argument("--carcinogen-col", default=None)
    p.add_argument("--gate-threshold", type=float, default=0.5)
    p.add_argument("--shortlist-size", type=int, default=100)
    p.add_argument("--out-full", required=True)
    p.add_argument("--out-shortlist", required=True)
    args = p.parse_args()

    with open(args.admet, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        required = [args.herg_col, args.ames_col, args.dili_col]
        missing = [c for c in required if c not in cols]
        if missing:
            sys.exit(f"ERROR: columns {missing} not found. Actual columns: {cols}")

        has_carc = args.carcinogen_col is not None and args.carcinogen_col in cols
        if args.carcinogen_col is not None and not has_carc:
            print(f"NOTE: carcinogen column '{args.carcinogen_col}' not found in {cols} — "
                  f"ranking proceeds WITHOUT it.", file=sys.stderr)

        rows = list(reader)

    gate_cols = [args.herg_col, args.ames_col, args.dili_col]
    gated_out, survived = [], []

    for row in rows:
        try:
            vals = [float(row[c]) for c in gate_cols]
        except ValueError:
            sys.exit(f"ERROR: non-numeric value in gate columns for row {row}")

        if any(v >= args.gate_threshold for v in vals):
            row["gate_status"] = "GATED_OUT"
            gated_out.append(row)
        else:
            score_cols = list(gate_cols) + ([args.carcinogen_col] if has_carc else [])
            row["gate_status"] = "SURVIVED"
            row["triage_score"] = f"{sum(1.0 - float(row[c]) for c in score_cols):.6f}"
            survived.append(row)

    survived.sort(key=lambda r: float(r["triage_score"]), reverse=True)
    shortlist = survived[: args.shortlist_size]

    full_fieldnames = cols + ["gate_status", "triage_score"]
    with open(args.out_full, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=full_fieldnames, restval="")
        writer.writeheader()
        writer.writerows(gated_out + survived)

    with open(args.out_shortlist, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=full_fieldnames, restval="")
        writer.writeheader()
        writer.writerows(shortlist)

    print(f"Gated out: {len(gated_out)}")
    print(f"Survived gate: {len(survived)}")
    print(f"Shortlist size: {len(shortlist)}")

if __name__ == "__main__":
    main()
