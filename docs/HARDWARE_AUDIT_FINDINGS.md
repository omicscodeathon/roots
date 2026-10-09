# Hardware audit findings (from hardware_audit.sh, run 2026-10-09)

`hardware_audit.sh` was run against every GROMACS log under
`Sys1+Sys2_results/` and `Sys3+Sys4_results/` (including the duplicated
`roots_md/` backup copies of Sys1 and Sys2). Every single log segment,
for all four systems, reports:

| Field | Value |
|---|---|
| GROMACS version | 2026.3-conda_forge |
| CPU | Intel(R) Xeon(R) CPU @ 2.00GHz |
| GPU | NVIDIA Tesla T4, compute cap. 7.5 |
| Threading | 1 MPI rank x 4 OpenMP threads |

**This contradicts the Lium H100 / GPU.ai A100-SXM4-80GB narrative in the
current README and prior notes.** Observed `calib`/`equil` throughputs in
these logs are ~7.6-31.2 ns/day on a shared T4, consistent with a
CX33-class or similarly modest cloud instance — not consistent with the
70-79 ns/day A100 figures or the 481.5 ns/day H100 figure quoted
elsewhere. Possible explanations: the H100/A100 providers were used for
a different run or a short window not captured in these logs, or the
final production trajectories were generated elsewhere and only the
calibration/equilibration stages ran on T4 instances that got bundled
into these result folders.

**Action before publishing:** reconcile which provider actually produced
the frames analyzed in `results/*/mmpbsa/*/FINAL_RESULTS_MMPBSA_*`. Do
not quote the 481.5 ns/day or 70-79 ns/day figures as the hardware behind
the final trajectories until this is resolved. Report the GPU honestly
as Tesla T4 for any stage where that's what the logs show.
