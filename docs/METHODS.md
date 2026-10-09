# Methods

## Molecular Dynamics

**System preparation.** Four complex-only systems were built in CHARMM-GUI
Solution Builder against the ALDH3B2 holo-NAD receptor, covering both
consensus and selective/flagged hits from the ADMET triage:

| System | Compound | Hit class |
|---|---|---|
| Sys1 | CNP0220660.1 | Consensus |
| Sys2 | CNP0358193.1 | Consensus |
| Sys3 | CNP0384182.2 | Selective / flagged |
| Sys4 | CNP0088713.1 | Selective / borderline |

**MD production.** GROMACS 2026.3 (conda-forge, CUDA-enabled). Each system
ran position-restrained minimization and NVT equilibration before
production, with periodic checkpointing (`-noappend`, so each restart
produced a new numbered segment). Production did **not** run on Colab;
rented GPU infrastructure was used throughout, with providers changing
over the course of the project (see Project History & Decisions for the
full pivot log). Target length was 50 ns per system. Status: **complete**
— all four systems reached their production target:

| System | Length reached | Notes |
|---|---|---|
| Sys1 | ~49.4-49.9 ns | GPU.ai A100-SXM4-80GB production; earlier segment on Lium H100 |
| Sys2 | ~49.4-49.9 ns | GPU.ai A100-SXM4-80GB production; earlier segment on a contended Lium node |
| Sys3 | ~49.4-49.9 ns | GPU.ai A100-SXM4-80GB production; earlier partial segment on Kaggle |
| Sys4 | ~49.4-49.9 ns | GPU.ai A100-SXM4-80GB production; earlier partial segment on Kaggle |

Note for authors: GROMACS logs for the calibration/equilibration stages
record an NVIDIA Tesla T4 GPU rather than the H100/A100 used for final
production. Reconcile and state per-stage hardware precisely before
submission, rather than naming one GPU for the whole run.

**Binding free energy refinement.** Status: **complete**, not "planned."
`gmx_MMPBSA` (Generalized Born, igb=8, and Poisson-Boltzmann) was run on a
Hetzner CPX51 (Ashburn; 16 shared vCPUs), over the converged window of
each trajectory (~34.4-49.4 ns), finishing 2026-10-09. Results, from
`results/mmpbsa/Sys{1..4}/FINAL_RESULTS_MMPBSA_*.dat`:

| System | Frames | ΔG (GB), kcal/mol | ΔG (PB), kcal/mol |
|---|---|---|---|
| Sys1 | 157 | -2822.02 ± 4.43 (SEM) | -2847.35 ± 4.39 (SEM) |
| Sys2 | 122 | -2601.36 ± 5.12 (SEM) | -2638.37 ± 4.91 (SEM) |
| Sys3 | 146 | -2962.78 ± 4.97 (SEM) | -2961.36 ± 5.12 (SEM) |
| Sys4 | 113 | -2656.62 ± 5.95 (SEM) | -2668.26 ± 5.95 (SEM) |

These are total complex energies (BOND+ANGLE+...+GSOLV), as reported by
gmx_MMPBSA — not ligand-receptor interaction/binding deltas. Compute a
proper ΔG_bind = G_complex - G_receptor - G_ligand from the per-component
COM/REC/LIG decomposition before quoting a binding affinity.

**Trajectory analysis and visualisation.** Status: **complete**, not
"depends on production completion." MDAnalysis and PLIP were run on a
representative frame from the converged window per system; matplotlib
and seaborn produced figures; Py3Dmol produced structure views; GROMACS
command-line tools (`gmx rms`, `gmx rmsf`, etc.) and a cross-target
comparative figure were also generated. Same Hetzner CPX51 instance as
the MM-PBSA step (now deleted).
