# ROOTS

*Receptor-oriented Oncotarget Screening of Natural Compounds Under Thermodynamic Simulation in Ovarian Cancer*

Supervisor: Dr. Olaitan Awe. Team Lead: Dr. Adedeji J.A.

Status: in execution. Steps 1 through 7 complete with real, generated data. Molecular dynamics production is live, running across four systems on rented GPU infrastructure.

---

## Overview

Ovarian cancer kills more people than most gynaecological cancers because it is caught late, and the tools that catch other cancers early do not exist for it yet. The therapies that work are effective where a patient can reach them, but their cost and distribution put them out of reach for underserved and low- and middle-income populations. What starts as a biology problem ends as a structural one.

Natural compounds and phytochemicals are a large, comparatively under-characterised chemical space. The field has produced some of oncology's most important drugs from this space, and most of it remains unexplored against specific cancer targets.

ROOTS narrows that gap using omics-derived target discovery and a mixed free-tier-and-rented compute stack to find natural-product candidates worth pursuing further, without assuming the resource-intensive screening infrastructure most pipelines start with. The pipeline is not a substitute for wet-lab validation or clinical development. It is a methodologically rigorous first step toward treatment options that do not assume the patient can already afford the alternative.

## Table of Contents

- [Overview](#overview)
- [Proposed Pipeline](#proposed-pipeline)
- [Methodology in Detail](#methodology-in-detail)
- [The Molecular Dynamics Bottleneck](#the-molecular-dynamics-bottleneck)
- [Existing Pipelines and the Gaps They Leave](#existing-pipelines-and-the-gaps-they-leave)
- [Points of Novelty](#points-of-novelty)
- [Validation and Benchmarking Strategy](#validation-and-benchmarking-strategy)
- [Testing and Validation Datasets](#testing-and-validation-datasets)
- [Compute and Infrastructure](#compute-and-infrastructure)
- [Execution Log: What Actually Happened](#execution-log-what-actually-happened)
- [Repository Structure](#repository-structure)
- [Current State](#current-state)
- [Team](#team)
- [References](#references)

---

## Proposed Pipeline

The pipeline runs in three phases. Each phase depends entirely on the output of the one before it.

![ROOTS three-phase pipeline](figures/ROOTS_Codeathon_Aug2026_Part1.png)

**Phase 1: Target and structure discovery**

TCGA-OV + GTEx + GEO, leading to target prioritisation, 2-3 candidate targets, structural characterisation, the COCONUT library pull and validated receptor branch, library curation, receptor and ligand preparation, and a system ready for docking.

**Phase 2: Screening and selection**

Docking validation, high-throughput docking, binding affinity estimation, hit definition and reconciliation, ADMET filtering, and a shortlist of one to two hits per target.

**Phase 3: Simulation, refinement, and output**

MD system preparation, MD production, binding free energy refinement, trajectory analysis and visualisation, comparative figures, and manuscript package assembly.

---

## Methodology in Detail

### Bioinformatics

- **Target prioritisation:** TCGA-OV versus GTEx-ovary differential expression through DESeq2/apeglm via recount3, with ComBat-seq mitigation, an independent GEO cohort for concordance (GSE26712 (185 ovarian cancer samples plus 10 normal ovarian epithelial samples, n=195) served as the independent GEO cohort, chosen for its adequate power and its cross-platform nature: microarray-based concordance with the RNA-seq discovery data is a stronger independence claim than same-platform concordance, since it rules out platform-specific measurement artefacts rather than merely repeating the same assay on different samples), GDC/cBioPortal mutation frequency, an HPA protein-level check, and a fallopian-tube tissue-of-origin screen. Output: ALDH3B2 carried forward as the lead target.
- **Structural characterisation:** PDB search first, gated on resolution, coverage, and organism; AlphaFold3 as fallback, validated by TM-score and pocket-level pLDDT. One structural state carried forward: the ALDH3B2 holo-NAD complex.
- **Compound library curation:** COCONUT as the sole source, archived at retrieval. Ro5, Veber, and PAINS filtering, then a cap of roughly 500-2,000 compounds. A stereochemistry-integrity check confirms defined stereocentres survive 3D embedding.
- **Receptor and ligand preparation:** CHARMM-GUI and PDB2PQR for protonation, Meeko for PDBQT conversion, the same protocol regardless of where the structure came from.

### Cheminformatics

- **Docking validation:** a per-target re-docking RMSD check using AutoDock Vina and MDAnalysis, run independently of the AlphaFold3 structural validation.
- **High-throughput docking:** 280,509 compounds docked with QuickVina2 across a 5-node Hetzner CX33 fleet; 41,529 compounds failed, overwhelmingly macrocyclic structures with glue atoms plus a smaller set with B/Si heteroatoms unsupported by the force field.
- **Binding affinity estimation:** PRODIGY-LIG scored across all 280,507 successfully docked poses, cross-validated against Vina's native score.
- **Hit definition and reconciliation:** percentile-based thresholding; top 1 percent (~2,805 compounds) defined as hits; 26 high-confidence consensus hits and 43 promiscuity-flagged compounds identified from the downstream ADMET-triaged set.
- **ADMET filtering:** ADMET-AI triaged 3,587 candidates locally with no cap; 280 survived to a smaller shortlist for SwissADME scrutiny.

### Molecular Dynamics

- **MD system preparation:** four complex-only systems built in CHARMM-GUI Solution Builder, covering both consensus and selective/flagged hits against the ALDH3B2 holo-NAD receptor: Sys1 (CNP0220660.1, consensus), Sys2 (CNP0358193.1, consensus), Sys3 (CNP0384182.2, selective/flagged), Sys4 (CNP0088713.1, selective/borderline).
- **MD production:** GROMACS 2026.3 (conda-forge, CUDA-enabled) on rented H100/A100 GPU infrastructure, not Colab (see the pivot log below for why). Position-restrained minimization and NVT equilibration precede each production run; trajectory target is 50 ns per system with periodic checkpointing.
- **Binding free energy refinement:** gmx_MMPBSA, planned on a separate CPU instance once production trajectories are complete. Not yet run.
- **Trajectory analysis and visualisation:** MDAnalysis and PLIP on a representative frame from the converged window, matplotlib and seaborn for figures, Py3Dmol for structure views, plus one cross-target comparative figure. Not yet run; depends on production completion.
- **Manuscript package assembly:** pending completion of the above.

---

## The Molecular Dynamics Bottleneck

MD gives the pipeline its name, and it was, in practice as well as in design, the real bottleneck. The original plan assumed Colab's free tier with Drive-synced checkpointing (Karagöl & Karagöl, 2024) would carry the full MD workload at zero cost. That plan did not survive first contact with the actual system size.

- Real system size, once solvated, landed at roughly 148,000 atoms per complex, a scale at which Colab's free-tier session limits and GPU allocation made reliable, unattended multi-hour production runs impractical. The pipeline pivoted to rented GPU infrastructure for this phase specifically (see the pivot log below).
- Trajectory length still follows convergence rather than an arbitrary fixed number: each system targets at least 50 ns, with a backbone-RMSD trailing-slope check as an operational screening signal, not a stopping rule, once a system crosses that mark.
- Scope remains narrowed upstream: four systems total, enhanced sampling and metadynamics set aside as future extensions.

The zero-cost ambition for the omics, cheminformatics, and docking phases held throughout. MD production is the one phase where real system scale forced a paid-infrastructure pivot, documented honestly below rather than smoothed over.

---

## Existing Pipelines and the Gaps They Leave

Comparable pipelines exist in the literature. None combines what ROOTS combines.

| Study | Target | Approach |
|---|---|---|
| Pham et al. (2026) | CDK12, pan-cancer including ovarian | TCGA/GTEx omics into docking, MD, and MM-PBSA against Vietnamese natural products |
| Shibu et al. (2026) | WNT5A, ovarian | GEO2R differential expression, AlphaFold structure, phytocompound docking, MD, FEL, DCCF, FEP |
| Masoudi et al. (2025) | CDK1/WEE1, ovarian | ML-based prioritisation, docking, MD; naringin identified as lead |
| Sharma et al. (2026) | PD-1 | Consensus docking across seven algorithms, 300 ns MD, DFT, ADME |
| Luo et al. (2022) | CDK4/6 | Pharmacophore filtering of marine natural products, ADMET, docking, MD |
| Liu & Ji (2026) | ClC-3, cervical | AI-driven docking, deep-learning rescoring, pharmacokinetic filtering, MD |

**The gaps across this body of work:**

- Every study fixes its target by literature precedent. None runs a multi-cohort omics prioritisation before structure work starts.
- Normal-tissue comparator design is weak or missing. TCGA-OV's near-absence of matched normal samples goes unaddressed, and tissue of origin (fallopian tube versus ovarian surface epithelium) is not raised anywhere in this set.
- No study reports an independent second-cohort concordance check. Differential expression comes from a single cohort throughout.
- None cross-checks a differential-expression target against independent protein-level data before moving to structure work.
- Where AlphaFold is used, confidence gating is loose, for example "the highest-confidence model" with no pocket-specific pLDDT threshold stated.
- Docking protocol validation, re-docking a known reference ligand to confirm RMSD, is not a described step in any of these pipelines.
- Stereochemistry integrity for natural-product libraries goes unaddressed everywhere, despite Ro5 and PAINS filtering being universal.
- No promiscuity check or consensus/selective hit framework exists. Each study reports a single best hit with no statistical threshold behind it.
- ADMET tooling limits, SwissADME and pkCSM cap submissions around 100-200 molecules, go unaddressed at any scale where they would matter.
- Trajectory lengths are fixed in advance (100, 200, 300 ns quoted as flat numbers), and no pipeline documents a verified, cost-transparent compute stack with reproducibility archiving.

### Comparison of Pipelines

| Gap in existing pipelines | How ROOTS addresses it |
|---|---|
| Target fixed by literature precedent before any omics | Full multi-cohort differential expression, mutation frequency, and protein-level prioritisation run first; ALDH3B2 carried forward on that evidence |
| Weak or absent matched-normal-tissue handling; tissue of origin ignored | Harmonised TCGA-OV versus GTEx-ovary via recount3, the confound stated plainly, plus a fallopian-tube-origin literature cross-check |
| Single-cohort differential expression, no independent concordance check | An independent GEO cohort confirmed the same direction of change before the target survived |
| No orthogonal protein-level confirmation | HPA immunohistochemistry cross-check on the shortlisted target |
| Loose AlphaFold confidence gating | Dual-layer validation: whole-protein TM-score plus pocket-lining pLDDT; PDB preferred under explicit resolution, coverage, and organism criteria |
| Docking protocol validity assumed, not tested | Per-target re-docking of a known ligand with an RMSD check before large-scale screening |
| Stereochemistry integrity unaddressed for natural products | A check confirmed defined stereocentres survive 3D embedding; failures logged and excluded |
| Single best hit reported, no promiscuity screen | A percentile-based hit threshold, a consensus-versus-selective split, and a promiscuity sanity check on the ADMET-triaged set |
| Web ADMET tools cap submissions near 100-200 molecules | Local ADMET-AI triaged 3,587 candidates with no cap; capped web tools handled the smaller final shortlist |
| Fixed trajectory lengths set in advance; no documented cost-transparent infrastructure | Trajectory length follows convergence; compute cost and the mid-project pivot off free-tier infrastructure are documented in full below, not hidden |

---

## Points of Novelty

No single tool in ROOTS is new. The contribution primarily sits in the MD validation.

- Target selection was omics-derived across multiple candidate lines rather than fixed in advance, which sets it apart from comparable published pipelines that start from one pre-chosen target.
- The comparator design corrects a weakness documented across ovarian cancer transcriptomic work: TCGA-OV's near-absence of matched normal samples. The harmonised TCGA-versus-GTEx design states the confound plainly and reinforces it with GEO concordance, the fallopian-tube check, and HPA confirmation.
- The pipeline carries its target through every evidentiary layer, from transcriptomic discovery to thermodynamic validation, where comparable work usually stops earlier or starts later.
- The full compute stack and its real costs, including where the original zero-cost plan held and where it did not, are documented transparently rather than presented as a solved problem in advance.

---

## Validation and Benchmarking Strategy

Validation is built into each computational layer rather than deferred to a single endpoint or wet-lab confirmation from future research efforts.

- Target selection is checked against an independent measurement type, HPA immunohistochemistry, rather than resting on transcriptomics alone, plus a fallopian-tube literature check for tissue-of-origin artefacts.
- Structural reliability is benchmarked with numbers: TM-score against known folds and pocket-level pLDDT for predicted structures, re-docking RMSD against a known ligand.
- Binding affinity does not rest on one method. PRODIGY-LIG and Vina's native score were checked against each other across all 280,507 scored poses; the correlation was weak (Spearman rho approximately -0.22, against a 12.7 kcal/mol Vina range versus a 1.1 kcal/mol PRODIGY-LIG range), and that disagreement is reported as a disclosed limitation, with Vina gating hit status rather than an artificially reconciled consensus score.
- Thermodynamic stability rests on demonstrated RMSD convergence, not an assumed trajectory length, and instability is treated as a genuine result rather than a failed run.
- What remains unvalidated at this stage, replicate-based MD convergence, MM-PBSA binding free energies, and any biological activity claim beyond computational ranking, is stated as a scope boundary, not glossed over.

---

## Testing and Validation Datasets

### Testing dataset

The data used to build and tune the pipeline.

- **Discovery layer:** TCGA-OV tumour versus GTEx ovary normal, harmonised through recount3.
- **Supporting evidence:** pre-called mutation data from GDC and cBioPortal, and Human Protein Atlas immunohistochemistry.
- **Screening library:** the COCONUT natural-product database, filtered to the working set that produced 280,509 successfully docked compounds.
- **Structural inputs:** the ALDH3B2 holo-NAD AlphaFold3 structure, plus a known co-crystallised reference ligand for the docking-protocol RMSD test.

### Independent validation dataset

This is separate from the testing data above and was chosen specifically to share no provenance with it.

- **Target-level cohort:** an ovarian cancer cohort outside TCGA, GTEx, and the internal GEO validation series, such as the Australian Ovarian Cancer Study, ICGC-OV, or a second, distinct GEO series. This checks whether the prioritised target's dysregulation holds in a population and sequencing effort with no shared lineage to the discovery data.
- **Orthogonal molecular layer:** CPTAC-OV mass-spectrometry proteomics (Zhang et al., 2016), a proteomic characterisation of 174 ovarian tumours also studied by TCGA. It overlaps TCGA-OV at the patient level through a different assay than the HPA immunohistochemistry already in use internally.
- **Functional layer:** DepMap CRISPR knockout gene-effect scores across ovarian cancer cell lines, testing whether the prioritised target is functionally required for viability, using a data type the omics and docking design never touches.
- **Compound and hit-level validation:** ChEMBL, BindingDB, or NPASS bioactivity records for the final hits or close structural analogues, checked against the same target family independent of the internal Vina and PRODIGY-LIG scores.

---

## Compute and Infrastructure

- **Google Colab and Kaggle:** used as originally designed for omics analysis, library curation, and receptor/ligand preparation, at zero cost throughout. The compound-library stereochemistry-integrity check, CPU-bound RDKit 3D embedding with no GPU-accelerated path, ran on Kaggle's CPU-only notebook environment specifically for its larger RAM ceiling relative to Colab's free tier, after Colab's ~12.7GB limit constrained the full, uncapped pass across the post-Ro5/Veber/PAINS pool.
- **Hetzner CX33 fleet (5 nodes):** used for the 280,509-compound docking run and the PRODIGY-LIG cross-validation scoring pass. Real cost, confirmed: $35 total, invoiced.
- **Rented GPU infrastructure (H100/A100, pay-per-second):** used for MD production specifically, after the real 148,000-atom system size made Colab's free-tier session and GPU limits impractical for unattended multi-hour runs. See the pivot log below for the full reasoning and the real cost range.
- **Reproducibility:** raw database pulls archived at the point of retrieval, fixed random seeds wherever applicable, every intermediate artefact and environment specification under version control, and the final package to be archived with a Zenodo DOI. Every MD system ships with a standalone `REPRODUCE_sysN.sh` script that rebuilds it from the committed CHARMM-GUI stage archives.

---

## Target Selection: The Full Rationale

Target prioritisation was not a single clean pass from omics to one receptor. Several real findings, including negative ones, shaped which target was actually carried forward, and they are recorded here rather than only in the final choice.

An initial PCA showed near-complete separation between TCGA tumours and GTEx normals (PC1 = 79% variance); this was identified as confounded by cohort rather than disease status alone, and no batch correction was applied or claimed as a resolution, the confound is stated plainly instead. A fallopian-tube epithelial contamination signal, detected via established FT marker genes (PAX8, WT1, FOXJ1, OVGP1, CLDN4, CD24, MSLN), prompted a reanalysis and reranking of candidate targets. Two independently constructed analysis pipelines, run in different computational environments, converged on concordant low-priority classifications for MT1HL1, DEFB126, CST4, and SCGB1D1, supporting those negative findings as robust rather than an artefact of one pipeline.

ALDH3B2 did not replicate in the independent fallopian-tube sensitivity analysis. Its inclusion was deliberately not justified via its FT-lineage relationship to CD24; instead its selection rests on independent discovery, enzymatic tractability, and exploratory novelty. CD24 itself was reclassified from a discovery-ranked docking candidate to a translational pathway benchmark, since it is a heavily glycosylated surface protein unsuited to conventional pocket-based small-molecule docking, while its relevance via the CD24-Siglec10 axis as an FDA-deployed drug target is separately noted. FOLR1, which narrowly missed the predeclared statistical threshold, was retained as a clinically validated benchmark on the strength of its FDA-approved precedent, explicitly excluded from influencing discovery-based target ranking. High AlphaFold3 confidence (0.96 ipTM) for select candidates was interpreted cautiously: a confident fold prediction supports plausibility, not an independently established binding pocket, pose, or mechanism. TRIM29 appears only as a brief illustrative example of dataset and tissue heterogeneity, given its documented but non-specific relevance to ovarian cancer, not as a prioritised candidate.

ALDH3B2's structural model was constructed independently via OpenFold3, not AlphaFold Server/DB, avoiding the licence restriction on docking use that comes with the latter, and cross-checked against the canonical UniProt sequence. NAD+ cofactor placement transferred from two independent experimental donors (PDB 5FHZ, 8DR9) disagreed by 4.89 Angstrom, above a pre-set 3.0 Angstrom threshold; the holo receptor was gated behind a documented parameterisation review rather than silently resolved, and only the apo structure cleared for docking in the interim. The docking protocol itself was validated before production use: re-docking NAD+ into its own resolved pocket gave a full-molecule RMSD of 2.5 Angstrom, technically failing a 2.0 Angstrom threshold, but regional decomposition showed both rigid ring systems within tolerance and the entire discrepancy attributable to the flexible, near-symmetric pyrophosphate linker. Both the failing and passing metrics are reported together, not the more favourable one substituted for the other.

Two independent literature checks were pre-registered and run against the natural-product screening hits before docking scores were seen. Quercetin, initially reported absent from the prepared library, was confirmed by direct SMILES filtering to be excluded specifically by the PAINS catechol_A alert, a disclosed methodological artefact rather than a data gap or an embedding failure. Beta-elemene was confirmed present as a single, correctly stereo-resolved compound, with its rigid, donor/acceptor-free hydrocarbon structure flagged in advance as a plausible source of non-specific rather than genuine affinity, regardless of its eventual docking rank.

A parallel benchmark against ALDH1A1 and ALDH3A1, targets with known inhibitors, using property-matched decoys from the same compound library, is under construction, since ALDH3B2 has no published inhibitors to validate against directly. Independent literature re-screening of the full ALDH3B2 corpus is also underway as a second-reviewer check, alongside two corpus-level findings slated for the manuscript: a species caveat on prior catalytic-activity claims (measured in mouse, not human), and a direct rebuttal, via existing mass-spectrometry evidence, of a separate report characterising human ALDH3B2 as a pseudogene.

**Current target plan, as of this stage:** ALDH3B2 is the primary docking target. FOLR1 is retained as the clinical and methodological benchmark, with PDB structures 4LRH (folate-bound) and 5IZQ (AGF183-bound) available as reference complexes for pocket definition and redocking, should FOLR1 docking proceed. MMP13 is a possible replicated secondary target, not yet confirmed for the docking run reported in this repository. GPR31 is treated as a future extension rather than part of this screening phase. ALDH3B2 and FOLR1 receptor preparation both fall under the Step 4 workflow; the docking, scoring, and MD production described elsewhere in this README reflect ALDH3B2 specifically, the target that was actually carried through to completion.

## Execution Log: What Actually Happened

This section exists because the original design document was written pre-execution, and real execution surfaced real problems worth documenting plainly rather than silently absorbing into a clean final narrative. Each entry is a genuine pivot, with a short, honest reason.

**PRODIGY-LIG vs. Vina correlation was weak, not strong.** Cross-validating binding-affinity estimates across all 280,507 scored poses produced a Spearman rho of approximately -0.22, with Vina's score spanning a 12.7 kcal/mol range against PRODIGY-LIG's narrow 1.1 kcal/mol range. *Why this matters:* rather than force a reconciled consensus score, this is reported as a disclosed limitation, and Vina alone gates hit status. A weak correlation between two legitimately different scoring approaches is a real result, not a pipeline failure.

**CHARMM-GUI system preparation surfaced three real coordinate-handling bugs, not one.** RDKit's PDB writer left ligand chain IDs blank; the ligand coordinates being carried forward were fresh RDKit-embedded geometry, not the actual docked pose; and the apo receptor structure used for merging was in a different coordinate frame than the one actually used for docking, a 57.79 angstrom maximum CA-CA discrepancy. *Why this matters:* each was caught by a direct geometric check (atom counts, inter-atomic distances, chain IDs) before being trusted, not assumed correct because the pipeline ran without error. Scripts and reproducibility documentation for all three fixes are committed alongside the final system-prep outputs.

**Colab was replaced by rented GPU infrastructure for MD production.** The original design assumed Colab's free tier, with Drive-synced checkpointing, could carry MD production at zero cost (Karagöl & Karagöl, 2024). Once the real system size was confirmed at roughly 148,000 atoms per solvated complex, this assumption did not hold for reliable, unattended multi-hour runs. *Why this matters:* this is the one phase where the zero-cost design goal was not achievable at the actual system scale encountered, not a failure of the Colab-GROMACS method itself, which remains valid at smaller system sizes. The pivot is documented here rather than left for a reviewer to discover unexplained.

**GPU rental platform selection required real, hands-on verification, not just published pricing.** More than a dozen GPU cloud providers were evaluated directly: several had no genuine self-serve path to raw compute despite marketing otherwise, one had shut down entirely, one had a confirmed open tooling bug that silently broke GPU-model filtering, and node-to-node throughput on the eventual platform varied by over 30-fold between two nodes of the identical advertised GPU type, one genuinely contended, one not. *Why this matters:* a calibration run (short, timed, throughput-measured) before committing to full production became a non-negotiable step for every node used, after this variance was observed directly rather than assumed away.

**A zero-balance billing event interrupted two production runs mid-flight.** The rented platform's billing model removes a pod, and any data that lived only on its local disk, immediately upon balance reaching zero, with no grace period. Two production runs were lost this way before auto top-up was configured correctly and periodic mid-run checkpointing (`-cpt 15`, a checkpoint every 15 minutes) was added to every subsequent production launch. *Why this matters:* this is reported honestly as real, nontrivial lost compute time and cost, not omitted. The fix, checkpointing plus a verified, adequately funded auto top-up, is now standard for every system launched.

**GROMACS invocation required several real, specific corrections discovered only by running it.** Position-restraint flags (`-r`) are mandatory as of GROMACS 2018 and were initially omitted; index-group flags (`-n index.ndx`) were needed for CHARMM-GUI's SOLU/SOLV temperature-coupling groups; the CHARMM-GUI `toppar/` directory must remain alongside `topol.top`, not be flattened into a separate run directory; and the conda-forge GROMACS environment script itself breaks under strict (`set -u`) shell mode. *Why this matters:* none of these are documented prominently in CHARMM-GUI's own output, and each produced a distinct, specific, correctable error. The corrected invocation sequence is preserved in the committed scripts so this does not need rediscovering.

---

## Repository Structure

```
roots/
|-- accessions/
|   `-- accessions.txt
|-- dashboard/
|   |-- Run_App.R
|   |-- Server.R
|   |-- UI.R
|   `-- dependences.R
|-- data/
|   |-- step1/           (target shortlist, DE results, provenance -- corrected final versions)
|   |-- step2/           (validated receptor structure metadata)
|   |-- step3/           (final compound library)
|   |-- step5/           (merged PRODIGY-LIG cross-validation scores)
|   |-- step6/           (extracted docking scores, failure log -- not the 280k raw PDBQTs)
|   `-- step7/           (ADMET triage results, hit shortlist, SMILES)
|-- docs/
|   |-- notes.txt
|   |-- papers_url.txt
|   |-- slides_url.txt
|   `-- execution_map.md
|-- figures/
|   `-- ROOTS_Codeathon_Aug2026_Part1.png
|-- notebooks/
|-- output/
|   `-- report.txt
|-- scripts/
|   |-- REPRODUCE_sys1.sh ... REPRODUCE_sys4.sh   (standalone MD system rebuild, per hit)
|   |-- fix_ligand_for_charmmgui.py
|   |-- reorder_sdf.py
|   |-- prodigy_cross_validation.py
|   |-- merge_prodigy.py
|   |-- admet_triage.py
|   |-- hetzner_provision.sh
|   `-- prodigy_hetzner_fleet.py
|-- src/
|-- vignettes/
|-- workflow/
|   |-- project_lifecycle.txt
|   `-- steps.txt
|-- environments/
|   |-- colab_environment.yml
|   `-- kaggle_environment.yml
|-- LICENSE
`-- README.md
```

Heavy, raw, or superseded artifacts (the 280k+ docking PDBQT files, raw AlphaFold3 job output, raw per-node PRODIGY and docking fleet tarballs, Python virtualenvs, and duplicate staging copies) are intentionally excluded from the repository. Everything excluded is either regeneratable from the committed scripts or was a working intermediate never meant to be a final artifact.

---

## Current State

Steps 1 through 7 are complete with real, generated data: target prioritisation, structural validation, library curation, docking (280,509 compounds), PRODIGY-LIG cross-validation, hit definition, and ADMET triage. Four MD systems are built and verified (correct coordinate frame, correct atom counts, correct bond orders, no steric clashes). Production simulation is live across all four systems on rented GPU infrastructure, targeting 50 ns per system with periodic checkpointing.

The compound count narrows at every stage, and it is worth stating plainly rather than leaving the attrition implicit. Target prioritisation carried ALDH3B2 forward from the TCGA-OV/GTEx/GEO comparative analysis as the single receptor of record; the AlphaFold3 holo-NAD complex was the structural state used throughout, since no PDB entry cleared the resolution, coverage, and organism gate applied upstream. The COCONUT natural-product library, after Ro5, Veber, PAINS, and stereochemistry-integrity filtering, produced the working screening set. Of that set, 280,509 compounds docked successfully against the ALDH3B2 holo receptor; 41,529 failed, overwhelmingly macrocyclic structures carrying glue atoms that the docking engine cannot process, with a smaller remainder failing on unsupported B/Si heteroatoms. Every successfully docked pose was scored twice, once by Vina's native function and once by PRODIGY-LIG, and the two were cross-validated against each other rather than averaged, surfacing the weak correlation reported above as a real finding rather than reconciling it away. Percentile-based hit definition on the Vina scores then identified roughly the top 1 percent, on the order of 2,800 compounds, as hits. That hit set fed a two-tier ADMET process: ADMET-AI triaged 3,587 candidates locally with no submission cap, from which 280 survived to a smaller shortlist narrow enough for SwissADME, pkCSM, and ADMETlab 2.0, the web tools that cap submissions in the 100-200 range. From that shortlist, 26 compounds were classified as high-confidence consensus hits and 43 were flagged for promiscuity on molecular weight, logP, and aromatic ring count grounds, a distinction the pipeline keeps visible rather than collapsing into one undifferentiated hit list. Four of those hits, two consensus and two selective or flagged, were carried into MD system preparation against the same ALDH3B2 holo-NAD receptor used throughout: CNP0220660.1 and CNP0358193.1 as the consensus pair, CNP0384182.2 and CNP0088713.1 as the selective and borderline pair. That four-system split was a deliberate choice, not an arbitrary one, since it lets the production trajectories speak to whether consensus and selective hits behave differently under thermodynamic simulation rather than assuming the docking-stage classification is the final word.

Real corrections were made and documented along the way rather than smoothed over, both during the original design-and-audit phase and during execution itself, as detailed in the Execution Log above:

- A transcriptomic comparator confound between TCGA-OV and GTEx was identified and stated plainly rather than assumed resolved.
- A tool-identity error between PRODIGY (protein-protein) and PRODIGY-LIG (protein-ligand) was caught and corrected before the cheminformatics phase.
- Three distinct coordinate-handling bugs in CHARMM-GUI system preparation were caught by direct geometric verification before being trusted.
- The original zero-cost assumption for MD production did not survive the real system scale encountered and was replaced, transparently, with rented GPU infrastructure.

Remaining work: completion of MD production across all four systems, MM-PBSA binding free energy refinement, trajectory analysis and visualization, the cross-target comparative figure, and final manuscript package assembly with Zenodo archiving.

---

## Team

- **Supervisor:** Dr. Olaitan Awe
- **Team Lead:** John Adeyemo Adedeji
- **Team Members:** Tammy Lily Shishodia, Elsayed Hussien, Onifade Olayinka
- **Roles under recruitment:** omics and statistics, structural and docking, MD and thermodynamics, integration and statistical review.

---

## References

Abraham, J., Adler, J., Dunger, J., Evans, R., Green, T., Pritzel, A., Ronneberger, O., Willmore, L., Ballard, A. J., Bambrick, J., Bodenstein, S. W., Evans, D. A., Hung, C.-C., O'Neill, M., Reiman, D., Tunyasuvunakool, K., Wu, Z., Žemgulytė, A., Arvaniti, E., ... Jumper, J. M. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, *630*(8016), 493-500. https://doi.org/10.1038/s41586-024-07487-w

Abraham, M. J., Murtola, T., Schulz, R., Páll, S., Smith, J. C., Hess, B., & Lindahl, E. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX*, *1-2*, 19-25. https://doi.org/10.1016/j.softx.2015.06.001

Adasme, M. F., Linnemann, K. L., Bolz, S. N., Kaiser, F., Salentin, S., Haupt, V. J., & Schroeder, M. (2021). PLIP 2021: Expanding the scope of the protein-ligand interaction profiler to DNA and RNA. *Nucleic Acids Research*, *49*(W1), W530-W534. https://doi.org/10.1093/nar/gkab294

Eberhardt, J., Santos-Martins, D., Tillack, A. F., & Forli, S. (2021). AutoDock Vina 1.2.0: New docking methods, expanded force field, and Python bindings. *Journal of Chemical Information and Modeling*, *61*(8), 3891-3898. https://doi.org/10.1021/acs.jcim.1c00203

Gowers, R. J., Linke, M., Barnoud, J., Reddy, T. J. E., Melo, M. N., Seyler, S. L., Domański, J., Dotson, D. L., Buchoux, S., Kenney, I. M., & Beckstein, O. (2016). MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall & S. Rostrup (Eds.), *Proceedings of the 15th Python in Science Conference* (pp. 98-105). SciPy. https://doi.org/10.25080/Majora-629e541a-00e

GTEx Consortium. (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*, *369*(6509), 1318-1330. https://doi.org/10.1126/science.aaz1776

Karagöl, T., & Karagöl, A. (2024). *Benchmarking GROMACS on optimized Colab processors and the flexibility of cloud computing for molecular dynamics*. bioRxiv. https://doi.org/10.1101/2024.11.14.623563

Liu, C., & Ji, C. (2026). AI-enhanced virtual screening identifies a potent small-molecule modulator of ClC-3 for cervical cancer drug discovery. *Frontiers in Bioinformatics*, *6*, Article 1887419. https://doi.org/10.3389/fbinf.2026.1887419

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, *15*(12), Article 550. https://doi.org/10.1186/s13059-014-0550-8

Luo, L., Wang, Q., & Liao, Y. (2022). The inhibitors of CDK4/6 from a library of marine compound database: A pharmacophore, ADMET, molecular docking and molecular dynamics study. *Marine Drugs*, *20*(5), Article 319. https://doi.org/10.3390/md20050319

Masoudi, M., Samadiafshar, S., Azizi, H., & Skutella, T. (2025). Integrated analysis, machine learning, molecular docking and dynamics of CDK1 inhibitors in epithelial ovarian cancer: A multifaceted approach towards targeted therapy. *International Journal of Molecular Sciences*, *26*(18), Article 9168. https://doi.org/10.3390/ijms26189168

Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *Journal of Computational Chemistry*, *32*(10), 2319-2327. https://doi.org/10.1002/jcc.21787

Pham, T. D., Nguyen, L. B. X., Luu, T. D., Nguyen, M. T., Le, N. T. H., & Phan Van, T. (2026). Pan-cancer multi-omics characterization of CDK12 and virtual screening of Vietnamese natural products for novel inhibitors. *Molecular Diversity*, *30*(4), 5945-5969. https://doi.org/10.1007/s11030-026-11477-7

Sharma, N., Debnath, A., Mazumder, R., Rai, P., & Singh, R. K. (2026). In search of novel PD1 inhibitor from natural products by high-throughput virtual screening and molecular dynamics simulation. *PLOS ONE*, *21*(4), Article e0339355. https://doi.org/10.1371/journal.pone.0339355

Shibu, S., Sharon, A., Nanda Kumar, S. K., Vasudevan, V., Muthukumar, T., Nanda Kumar, M. P., Kumar, D. T., & Magesh, R. (2026). Integrative transcriptomics and structure-based screening identifies Phyllanthus amarus phytocompounds as potential WNT5A modulators in ovarian cancer. *Frontiers in Bioinformatics*, *6*, Article 1815548. https://doi.org/10.3389/fbinf.2026.1815548

Sorokina, M., Merseburger, P., Rajan, K., Yirik, M. A., & Steinbeck, C. (2021). COCONUT online: Collection of Open Natural Products database. *Journal of Cheminformatics*, *13*(1), Article 2. https://doi.org/10.1186/s13321-020-00478-9

Uhlén, M., Fagerberg, L., Hallström, B. M., Lindskog, C., Oksvold, P., Mardinoglu, A., Sivertsson, Å., Kampf, C., Sjöstedt, E., Asplund, A., Olsson, I., Edlund, K., Lundberg, E., Navani, S., Szigyarto, C. A.-K., Odeberg, J., Djureinovic, D., Takanen, J. O., Hober, S., ... Pontén, F. (2015). Tissue-based map of the human proteome. *Science*, *347*(6220), Article 1260419. https://doi.org/10.1126/science.1260419

Valdés-Tresanco, M. S., Valdés-Tresanco, M. E., Valiente, P. A., & Moreno, E. (2021). gmx_MMPBSA: A new tool to perform end-state free energy calculations with GROMACS. *Journal of Chemical Theory and Computation*, *17*(10), 6281-6291. https://doi.org/10.1021/acs.jctc.1c00645

Vangone, A., & Bonvin, A. M. J. J. (2015). Contacts-based prediction of binding affinity in protein-protein complexes. *eLife*, *4*, Article e07454. https://doi.org/10.7554/eLife.07454

Vangone, A., Schaarschmidt, J., Koukos, P., Geng, C., Citro, N., Trellet, M. E., Xue, L. C., & Bonvin, A. M. J. J. (2019). Large-scale prediction of binding affinity in protein-small ligand complexes: The PRODIGY-LIG web server. *Bioinformatics*, *35*(9), 1585-1587. https://doi.org/10.1093/bioinformatics/bty816

Wilks, C., Zheng, S. C., Chen, F. Y., Charles, R., Solomon, B., Ling, J. P., Imada, E. L., Zhang, D., Joseph, L., Leek, J. T., Jaffe, A. E., Nellore, A., Collado-Torres, L., Hansen, K. D., & Langmead, B. (2021). recount3: Summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biology*, *22*(1), Article 323. https://doi.org/10.1186/s13059-021-02533-6

Xue, L. C., Rodrigues, J. P. G. L. M., Kastritis, P. L., Bonvin, A. M. J. J., & Vangone, A. (2016). PRODIGY: A web server for predicting the binding affinity of protein-protein complexes. *Bioinformatics*, *32*(23), 3676-3678. https://doi.org/10.1093/bioinformatics/btw514

Zhang, H., Liu, T., Zhang, Z., Payne, S. H., Zhang, B., McDermott, J. E., Zhou, J.-Y., Petyuk, V. A., Chen, L., Ray, D., Sun, S., Yang, F., Chen, L., Wang, J., Shah, P., Cha, S. W., Aiyetan, P., Woo, S., Tian, Y., ... Rodland, K. D. (2016). Integrated proteogenomic characterization of human high-grade serous ovarian cancer. *Cell*, *166*(3), 755-765. https://doi.org/10.1016/j.cell.2016.05.069

Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: Batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*, *2*(3), Article lqaa078. https://doi.org/10.1093/nargab/lqaa078

---

*This document reflects the pipeline's state as of the most recent execution update. It will be updated again as MD production, MM-PBSA refinement, and final analysis complete.*
