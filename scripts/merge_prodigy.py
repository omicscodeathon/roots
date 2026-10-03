#!/usr/bin/env python3
"""Merge and verify the 5 pulled PRODIGY-LIG fleet results into one final scores file.

Checks for duplicate compound_ids across instances (would indicate an overlap/reconciliation
bug from the docking-stage rebalance), reports per-instance and combined stats, and writes a
single merged CSV ready for downstream ADMET/hit-reconciliation steps.

Usage: python3 merge_prodigy.py [--results-dir ~/roots_pipeline1/prodigy_results] [--out prodigy_scores_merged.csv]
"""
import argparse
import csv
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=str(Path.home() / "roots_pipeline1" / "prodigy_results"))
    ap.add_argument("--out", default=str(Path.home() / "roots_pipeline1" / "prodigy_scores_merged.csv"))
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    instance_dirs = sorted(results_dir.glob("roots-prodigy-*"))
    if not instance_dirs:
        raise SystemExit(f"No roots-prodigy-* directories found in {results_dir}")

    all_rows = []
    seen_ids = Counter()
    per_instance_stats = []

    for d in instance_dirs:
        csv_path = d / "prodigy_scores.csv"
        verified = (d / "VERIFIED.json").exists()
        if not csv_path.exists():
            print(f"WARNING: {d.name} has no prodigy_scores.csv, skipping")
            continue
        with open(csv_path, newline="") as f:
            rows = list(csv.DictReader(f))
        ok_rows = [r for r in rows if r["status"] == "ok"]
        for r in rows:
            seen_ids[r["compound_id"]] += 1
            r["source_instance"] = d.name
            all_rows.append(r)
        per_instance_stats.append({
            "instance": d.name, "verified": verified, "total": len(rows),
            "ok": len(ok_rows), "failed": len(rows) - len(ok_rows),
        })

    print(f"{'instance':<20}{'verified':>10}{'total':>10}{'ok':>10}{'failed':>10}")
    grand_total = grand_ok = grand_failed = 0
    for s in per_instance_stats:
        print(f"{s['instance']:<20}{str(s['verified']):>10}{s['total']:>10,}{s['ok']:>10,}{s['failed']:>10,}")
        grand_total += s["total"]
        grand_ok += s["ok"]
        grand_failed += s["failed"]
    print(f"{'TOTAL':<20}{'':>10}{grand_total:>10,}{grand_ok:>10,}{grand_failed:>10,}")

    dupes = {cid: n for cid, n in seen_ids.items() if n > 1}
    print(f"\nDuplicate compound_ids across instances: {len(dupes)}")
    if dupes:
        print("  This means the same compound was scored by more than one instance -- check")
        print("  whether the docking-stage rebalance's chunk assignment overlapped. Sample:")
        for cid in list(dupes)[:10]:
            print(f"    {cid}: appears {dupes[cid]} times")

    ok_vals = sorted(float(r["prodigy_dg"]) for r in all_rows if r["status"] == "ok")
    n = len(ok_vals)
    if n:
        print(f"\nCombined PRODIGY-LIG dG stats (n={n:,}):")
        print(f"  min={ok_vals[0]:.2f}  max={ok_vals[-1]:.2f}  "
              f"median={ok_vals[n // 2]:.2f}  mean={sum(ok_vals) / n:.2f}")

    fail_reasons = Counter(r["status"] for r in all_rows if r["status"] != "ok")
    if fail_reasons:
        print(f"\nFailure reasons across the fleet:")
        for reason, count in fail_reasons.most_common():
            print(f"  {count:>6}  {reason}")

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["compound_id", "prodigy_dg", "status", "source_instance"])
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nMerged file written: {args.out} ({len(all_rows):,} rows)")


if __name__ == "__main__":
    main()
