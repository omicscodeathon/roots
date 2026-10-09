# RedDustReclaimer

A small, reproducible computational study of whether an engineered perchlorate
detoxification pathway could benefit *Bacillus subtilis* under the
oxygen-limited conditions relevant to Martian regolith.

Perchlorate salts are abundant in Martian soil. The bacterial enzymes that
break them down — perchlorate reductase (PcrAB) and chlorite dismutase (Cld) —
are notable because the second step releases molecular oxygen. This repository
integrates that pathway into the published genome-scale metabolic model of
*B. subtilis* (iYO844) and asks a focused question: does the pathway change
predicted growth, and if so, under what conditions?

**Result in one line:** the pathway is metabolically silent when oxygen is
plentiful, but under oxygen limitation the oxygen released by chlorite
dismutase supplements respiration and raises predicted growth by 25–200% —
and below the model's aerobic feasibility threshold it enables growth that is
otherwise impossible. The benefit is *oxygen self-supply*, not detoxification
per se, which is precisely why it matters for Mars' oxygen-poor atmosphere.

---

## Contents

```
reddust/
  perchlorate_fba.py   pathway integration + growth / pFBA / NGAM analysis for iYO844
  codon_usage.py       codon-usage scoring and synonymous recoding
scripts/
  01_codon_analysis.py            -> results/codon_usage.csv
  02_perchlorate_o2_analysis.py   -> results/perchlorate_o2_benefit.csv
  03_sensitivity_checks.py        -> results/sensitivity_checks.csv
models/iYO844.xml      base genome-scale model (Oh et al., 2007)
data/                  candidate coding sequences and enzyme references
tests/                 unit tests
paper/manuscript.md    the write-up
```

## Reproducing the results

```bash
# environment (conda/micromamba recommended)
micromamba create -y -f environment.yml
micromamba activate reddust
# or: pip install -r requirements.txt

python scripts/01_codon_analysis.py
python scripts/02_perchlorate_o2_analysis.py
python scripts/03_sensitivity_checks.py
pytest -q
```

Each script writes a CSV to `results/` and prints its key numbers. The pathway
reactions are checked for elemental and charge balance at build time, so an
unbalanced edit cannot pass silently. The 9-test suite verifies model
dimensions, reaction balance, and the neutral-when-aerobic /
beneficial-when-limited behaviour.

---

## The pathway

Two reactions are added to iYO844, both elementally and charge balanced
(verified programmatically at construction):

```
PCR (perchlorate reductase, pcrAB):   ClO4- + 2 FADH2  ->  ClO2- + 2 H2O + 2 FAD
CLD (chlorite dismutase, cld):        ClO2-            ->  Cl- + O2
```

The decisive feature is that chlorite dismutase releases molecular oxygen into
the host's cytosolic O₂ pool, where it can serve as a terminal electron
acceptor.

**Base model.** iYO844 — 844 genes, 1250 reactions, 990 metabolites, objective
`BIOMASS_BS_10`. Non-growth maintenance (`ATPM`) fixed at 9.0 mmol gDW⁻¹ h⁻¹.
Baseline growth at the default glucose bound (−1.7): **0.11797 h⁻¹**.

---

## Headline result

Predicted specific growth rate (h⁻¹) as the external oxygen uptake cap is
tightened, with and without perchlorate available (perchlorate uptake cap = 20;
PCR and CLD carry equal flux):

| External O₂ cap (mmol gDW⁻¹ h⁻¹) | Growth, no perchlorate | Growth, perchlorate available | PCR = CLD flux | Change |
|---|---|---|---|---|
| unlimited | 0.11797 | 0.11797 | 0.0000 | none |
| 6 | 0.11797 | 0.11797 | 0.0000 | none |
| 4 | 0.07993 | 0.10589 | 0.9699 | +32% |
| 2 | 0.02524 | 0.03162 | 0.4910 | +25% |
| 1.5 | 0.00854 | 0.01305 | 0.3713 | +53% |
| 1.3 | 0.00186 | 0.00562 | 0.3234 | +202% |
| 1.25 | 0.00019 | 0.00376 | 0.3114 | +1879% |
| 1.2 | infeasible | 0.00190 | 0.2995 | **growth rescued** |

Feasibility threshold without perchlorate: feasible at an O₂ cap of 1.25
(growth 0.0002); infeasible at 1.2444 and below. Flux is reported only for
solver-optimal solutions.

### Mechanism

At an external-O₂ cap of 1.3, the cell draws its full external oxygen allowance
and chlorite dismutase supplies more oxygen internally:

| Quantity (mmol gDW⁻¹ h⁻¹) | Value |
|---|---|
| External O₂ uptake | 1.3000 |
| CLD-generated O₂ | 0.3234 |
| Total O₂ to cell | 1.6234 |
| Increase over external alone | **+25%** |

Equal PCR and CLD flux across every condition confirms flux runs through the
complete pathway to the oxygen-releasing step, not through perchlorate
reduction alone.

---

## Robustness checks

**Parsimonious FBA (pFBA).** To confirm the result is not an artifact of one
arbitrary flux distribution among many optima, each scenario was re-solved with
pFBA (minimum total flux at the growth optimum). Growth rates matched plain FBA
to within 1e-6 at **all 16 conditions** (8 oxygen caps × 2 perchlorate
conditions), including agreement on infeasibility at the rescue point.

**NGAM sensitivity.** Non-growth maintenance (`ATPM`) was swept ±20% at the
rescue scenario (O₂ cap = 1.2). The strict rescue (infeasible without
perchlorate, feasible with) appears at the default maintenance and shifts
predictably — at lower maintenance the cell already survives without
perchlorate; at higher maintenance the point falls below survivability
entirely:

| Δ NGAM | ATPM flux | Growth, no perchlorate | Growth, perchlorate available | PCR = CLD | Outcome |
|---|---|---|---|---|---|
| −20% | 7.200 | 0.01188 | 0.01562 | 0.2956 | both grow |
| −10% | 8.100 | 0.00520 | 0.00876 | 0.2975 | both grow |
| 0% (default) | 9.000 | infeasible | 0.00190 | 0.2995 | **rescue** |
| +10% | 9.900 | infeasible | infeasible | — | below threshold |
| +20% | 10.800 | infeasible | infeasible | — | below threshold |

The underlying *benefit* (perchlorate improves growth wherever the cell is
viable) is robust across the whole range; the specific O₂ value at which
benefit becomes life-or-death rescue tracks the maintenance assumption rather
than being fixed at 1.2.

**Flux variability (FVA).** The pathway is optional under oxygen-replete
conditions but required under limitation. Aerobically (O₂ cap 20, full feasible
range): PCR spans [0.0, 3.4], CLD [0.0, 20.4]. Under limitation (O₂ cap 1.3, at
90% of the growth optimum): both are pinned to [0.27, 0.36] — nonzero minimum
flux, i.e. no longer skippable.

---

## Supporting analysis: codon usage

A frequency-based adaptation score (documented reference table; not a
species-calibrated Codon Adaptation Index) for six candidate coding sequences,
with a synonymous recoding that preserves each encoded protein exactly. Effects
are modest, as expected from a single-reference heuristic; reported to document
the sequences and the recoding utility, not as an expression claim.

| Sequence | bp | GC % | Baseline | Optimized | Δ |
|---|---|---|---|---|---|
| *Deinococcus radiodurans* RecA | 300 | 65.33 | 0.442 | 0.459 | +4.1% |
| *Psychrobacter arcticus* cold-shock protein | 180 | 83.89 | 0.470 | 0.482 | +2.5% |
| *Methanocaldococcus jannaschii* carbonic anhydrase | 297 | 65.32 | 0.444 | 0.462 | +4.1% |
| *Caldalkalibacillus thermarum* alkaline protease | 180 | 86.67 | 0.470 | 0.475 | +1.1% |
| *Halobacterium salinarum* bacteriorhodopsin | 183 | 85.79 | 0.470 | 0.476 | +1.4% |
| *Thermus thermophilus* DNA polymerase | 180 | 85.56 | 0.470 | 0.478 | +1.8% |

---

## Model integrity

Two defects in an early pathway implementation were found and resolved before
any result was trusted:

- **Free-energy artifact (fixed).** An ATP-coupling reaction had no link to
  pathway flux and synthesized ATP from nothing; it inflated growth from 0.118
  to 0.267 h⁻¹ *even with pathway flux forced to zero*. Removed. The genuine
  benefit reported here is entirely oxygen-driven.
- **Mass balance (fixed).** The reductase and dismutase reactions were
  rebalanced to the correct stoichiometry (shown above); both now pass
  automated element- and charge-balance checks, enforced on every build.

---

## Scope and limitations

This is a constraint-based (flux balance) modeling study. It predicts what is
metabolically feasible and favourable under stated assumptions; it does not
measure enzyme kinetics, perchlorate or chlorite toxicity, expression burden,
or in-vivo behaviour. The oxygen benefit is a stoichiometric upper bound — the
cell *can* use CLD-derived oxygen to grow more, not a guarantee that a real
strain would. The codon-usage tool is a frequency-based heuristic against a
documented reference table, not a species-calibrated CAI.

The specific, testable prediction: **perchlorate availability should improve
growth of a PcrAB/Cld-expressing *B. subtilis* under microaerophilic but not
aerobic culture, in proportion to the severity of oxygen limitation.** Results
are a hypothesis to be tested experimentally, not a validated strain design.

---

## References

1. Hecht MH, Kounaves SP, Quinn RC, et al. Detection of perchlorate and the
   soluble chemistry of Martian soil at the Phoenix lander site. *Science*.
   2009;325(5936):64–67.
2. Coates JD, Achenbach LA. Microbial perchlorate reduction: rocket-fuelled
   metabolism. *Nature Reviews Microbiology*. 2004;2(7):569–580.
3. Lee AQ, Streit BR, Zdilla MJ, Abu-Omar MM, DuBois JL. Mechanism of and
   exquisite selectivity for O–O bond formation by the heme-dependent chlorite
   dismutase. *PNAS*. 2008;105(41):15654–15659.
4. Oh YK, Palsson BO, Park SM, Schilling CH, Mahadevan R. Genome-scale
   reconstruction of the metabolic network in *Bacillus subtilis*. *Journal of
   Biological Chemistry*. 2007;282(39):28791–28799.
5. Orth JD, Thiele I, Palsson BO. What is flux balance analysis?
   *Nature Biotechnology*. 2010;28(3):245–248.

## License

MIT — see `LICENSE`.
