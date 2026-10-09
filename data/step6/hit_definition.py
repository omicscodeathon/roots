#!/usr/bin/env python3
"""
ROOTS -- Hit Definition and Reconciliation

Defines "hit" status from the Vina/QuickVina2 score distribution (percentile-based, computed
from the real, measured distribution of this screen, not an assumed fixed cutoff, per the
earlier confirmed finding that PRODIGY-LIG's ~1.1 kcal/mol range across 280,507 compounds gives
it no discriminative power here). PRODIGY-LIG is used as a post-hoc sanity/promiscuity check on
selected hits, not as a second independent gate.

Tiers:
  - "hit"            top HIT_PERCENTILE (default 1%) of the real Vina score distribution
  - "high_confidence" top HIGH_CONF_PERCENTILE (default 0.1%) of the same distribution,
                      a subset of "hit"

Promiscuity sanity check: for each hit, computes where its PRODIGY-LIG score falls within
PRODIGY-LIG's own (narrow) observed range. A hit whose Vina score is very favorable but whose
PRODIGY-LIG score sits at the UNfavorable end of PRODIGY-LIG's own narrow range is flagged as
"promiscuity_flag" -- consistent with the beta-elemene concern raised earlier (a greasy,
non-specific binder scoring well in Vina via raw contact count, without genuine complementarity).
This is a flag for manual review, not an automatic exclusion.

Usage:
    python3 hit_definition.py --vina vina_scores.csv --prodigy prodigy_scores_merged.csv \
        --out hit_reconciliation.csv --hit-percentile 1.0 --high-conf-percentile 0.1
"""
import argparse
import csv
import statistics
import sys
from pathlib import Path


def load_vina(path):
    scores = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        id_col = "compound_id" if "compound_id" in reader.fieldnames else reader.fieldnames[0]
        score_col = next((c for c in reader.fieldnames if "vina" in c.lower() or "score" in c.lower()
                          or "affinity" in c.lower()), None)
        if score_col is None:
            sys.exit(f"Could not find a score column in {path}. Columns found: {reader.fieldnames}")
        for row in reader:
            try:
                scores[row[id_col]] = float(row[score_col])
            except (ValueError, KeyError):
                continue
    return scores


def load_prodigy(path):
    scores = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("status") == "ok":
                try:
                    scores[row["compound_id"]] = float(row["prodigy_dg"])
                except ValueError:
                    continue
    return scores


def percentile_threshold(values, percentile):
    """The Vina score at which the top `percentile`% of compounds (by favorability, i.e. most
    negative) begins. E.g. percentile=1.0 -> the score at the 1st percentile of the sorted
    (ascending, most-favorable-first) distribution."""
    sorted_vals = sorted(values)  # ascending: most negative (most favorable) first
    n = len(sorted_vals)
    idx = max(0, int(n * percentile / 100) - 1)
    return sorted_vals[idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vina", required=True)
    ap.add_argument("--prodigy", required=True)
    ap.add_argument("--out", default="hit_reconciliation.csv")
    ap.add_argument("--hit-percentile", type=float, default=1.0,
                    help="top N%% of the real Vina distribution = 'hit' (default 1.0)")
    ap.add_argument("--high-conf-percentile", type=float, default=0.1,
                    help="top N%% = 'high_confidence' subset of hit (default 0.1)")
    args = ap.parse_args()

    for label, p in (("--vina", args.vina), ("--prodigy", args.prodigy)):
        if not Path(p).exists():
            sys.exit(f"{label} file not found: {p}")

    print(f"Loading Vina scores from {args.vina}...")
    vina = load_vina(args.vina)
    print(f"  {len(vina):,} Vina scores loaded")

    print(f"Loading PRODIGY-LIG scores from {args.prodigy}...")
    prodigy = load_prodigy(args.prodigy)
    print(f"  {len(prodigy):,} PRODIGY-LIG (ok) scores loaded")

    vina_vals = list(vina.values())
    n_total = len(vina_vals)
    vmin, vmax = min(vina_vals), max(vina_vals)
    vmedian = statistics.median(vina_vals)
    print(f"\nReal Vina distribution (n={n_total:,}): min={vmin:.2f} max={vmax:.2f} median={vmedian:.2f}")

    hit_threshold = percentile_threshold(vina_vals, args.hit_percentile)
    high_conf_threshold = percentile_threshold(vina_vals, args.high_conf_percentile)
    print(f"\nHit threshold (top {args.hit_percentile}% of real distribution): Vina <= {hit_threshold:.2f}")
    print(f"High-confidence threshold (top {args.high_conf_percentile}%): Vina <= {high_conf_threshold:.2f}")

    # PRODIGY-LIG's own observed range, for the promiscuity sanity check
    prodigy_vals = list(prodigy.values())
    if prodigy_vals:
        p_min, p_max = min(prodigy_vals), max(prodigy_vals)
        p_mid = (p_min + p_max) / 2
        print(f"PRODIGY-LIG observed range: {p_min:.2f} to {p_max:.2f} (midpoint {p_mid:.2f})")
    else:
        p_min = p_max = p_mid = None
        print("WARNING: no PRODIGY-LIG scores loaded; promiscuity flag will be unavailable for all hits.")

    rows = []
    n_hit = n_high_conf = n_flagged = n_no_prodigy = 0
    for cid, v_score in vina.items():
        is_hit = v_score <= hit_threshold
        if not is_hit:
            continue
        is_high_conf = v_score <= high_conf_threshold
        p_score = prodigy.get(cid)
        promiscuity_flag = ""
        if p_score is not None and p_mid is not None:
            # Flagged if PRODIGY-LIG places this compound in the LESS favorable half of its own
            # narrow range, despite a strongly favorable Vina score -- a real, disclosed proxy
            # for "favorable Vina score not corroborated even weakly by an independent estimator",
            # not a claim that PRODIGY-LIG's absolute value is itself meaningful here.
            promiscuity_flag = "flag" if p_score > p_mid else ""
            if promiscuity_flag:
                n_flagged += 1
        else:
            n_no_prodigy += 1
        rows.append({
            "compound_id": cid, "vina_score": v_score, "prodigy_dg": p_score if p_score is not None else "",
            "is_hit": True, "is_high_confidence": is_high_conf, "promiscuity_flag": promiscuity_flag,
        })
        n_hit += 1
        if is_high_conf:
            n_high_conf += 1

    rows.sort(key=lambda r: r["vina_score"])  # most favorable first

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["compound_id", "vina_score", "prodigy_dg",
                                                "is_hit", "is_high_confidence", "promiscuity_flag"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n=== Hit reconciliation summary ===")
    print(f"Total compounds screened: {n_total:,}")
    print(f"Hits (top {args.hit_percentile}%): {n_hit:,}")
    print(f"  of which high-confidence (top {args.high_conf_percentile}%): {n_high_conf:,}")
    print(f"  of which promiscuity-flagged: {n_flagged:,} ({100*n_flagged/n_hit:.1f}% of hits)" if n_hit else "")
    print(f"  hits missing a PRODIGY-LIG score: {n_no_prodigy:,}")
    print(f"\nOutput: {args.out}")
    print(f"\nNOTE: promiscuity_flag is a disclosed sanity-check heuristic, not an exclusion.")
    print(f"Given PRODIGY-LIG's own narrow observed range on this screen ({p_min:.2f} to {p_max:.2f} "
          f"if loaded), this flag should be reported and interpreted cautiously, not treated as a")
    print(f"confident independent confirmation either way.")


if __name__ == "__main__":
    main()
