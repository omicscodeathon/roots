<div align="center">

# ROOTS

**Receptor-oriented Oncotarget Screening of Natural Compounds Under Thermodynamic Simulation in Ovarian Cancer**

![License](https://img.shields.io/badge/license-MIT-blue) ![Python](https://img.shields.io/badge/python-3.12-blue) ![GROMACS](https://img.shields.io/badge/GROMACS-2026.3-orange) ![Status](https://img.shields.io/badge/status-MD%20%26%20MM--PBSA%20complete-success)

Supervisor: Dr. Olaitan Awe · Team Lead: John Adeyemo Adedeji

**[Project History & Decisions](https://github.com/omicscodeathon/roots/blob/main/docs/PROJECT_HISTORY_AND_DECISIONS.md)** — the full pivot log, target-selection rationale, and execution record behind every number on this page.

</div>

---

## Table of Contents

- [Overview](#overview)
- [Objectives](#objectives)
- [Background: Ovarian Cancer and the Case for ALDH3B2](#background-ovarian-cancer-and-the-case-for-aldh3b2)
- [Pipeline](#pipeline)
- [Data Sources](#data-sources)
- [Methods](#methods)
- [Results](#results)
- [What Sets ROOTS Apart](#what-sets-roots-apart)
- [Reproducibility](#reproducibility)
- [Repository Structure](#repository-structure)
- [Team](#team)
- [Acknowledgments](#acknowledgments)
- [References](#references)

---

## Overview

Ovarian cancer kills more people than most gynaecological cancers because it is caught late, and the screening tools that catch other cancers early do not exist for it yet. The therapies that do work are effective where a patient can reach them, but cost and distribution put them out of reach for underserved and low- and middle-income populations. A biology problem ends as a structural one.

Natural compounds and phytochemicals remain a large, under-characterised chemical space relative to their track record: the field has produced some of oncology's most important drugs from this space, and most of it has never been screened against a specific cancer target. ROOTS narrows that gap for ovarian cancer using omics-derived target discovery, structure-based virtual screening, and molecular dynamics, run on a mixed free-tier-and-rented compute stack rather than assuming the resource-intensive infrastructure most pipelines start with. The pipeline is not a substitute for wet-lab validation or clinical development. It is a methodologically rigorous first step toward treatment candidates that do not assume the patient can already afford the alternative.

## Objectives

- Prioritise a druggable ovarian cancer target through multi-cohort differential expression, independent-cohort concordance, mutation frequency, and protein-level confirmation, rather than fixing a target by literature precedent.
- Screen a large natural-product library against that target's structure and reduce it to a reconciled, promiscuity-screened hit set.
- Carry the surviving hits through molecular dynamics and binding free energy estimation to determine which are worth pursuing experimentally.

## Background: Ovarian Cancer and the Case for ALDH3B2

Ovarian cancer was the eighth most common cancer in women globally in 2020 (Webb & Jordan, 2024), and GBD 2019 attributes 294,422 incident cases, 198,412 deaths, and 5.36 million DALYs to it, with the absolute burden rising since 1990 (Zhang et al., 2022). Roughly 80% of patients present at stage III-IV (Caruso et al., 2025), and the five-year survival rate for advanced disease sits between 10% and 40% despite platinum-based chemotherapy remaining standard of care for decades (Richardson et al., 2023). High-grade serous carcinoma, the dominant histotype, resolves into at least four transcriptional subtypes on its own ("Integrated genomic analyses of ovarian carcinoma," 2011), and single-cell and spatial work keeps finding more heterogeneity within individual tumours, not less (Luo et al., 2024). A single-target, single-mechanism approach is a poor fit for a disease that looks this different from one cohort, and one cell, to the next.

Target discovery for ovarian cancer carries two recurring, correctable weaknesses in the published literature. First, TCGA-OV has almost no matched normal tissue, and the standard substitute, GTEx ovarian stroma, is not the tissue most high-grade serous disease is now thought to originate from; O'Shannessy et al. (2013) showed that genes reported as overexpressed against ovarian stroma lost that signal entirely against fallopian tube, and more recent single-cell work finds fallopian-tube progenitor states with direct transcriptional continuity into HGSOC (Li et al., 2026). Second, differential expression from a single cohort is routinely carried straight into structure-based work with no independent concordance check and no protein-level confirmation. ROOTS' own target-prioritisation pipeline is built against both of these points directly, documented in full in the [Project History](https://github.com/omicscodeathon/roots/blob/main/docs/PROJECT_HISTORY_AND_DECISIONS.md): a harmonised TCGA-OV-versus-GTEx design that states the cohort confound plainly rather than batch-correcting it away, an independent GEO microarray cohort for concordance, and an HPA protein-level cross-check.

ALDH3B2 is the target that survived that process, and it is worth being precise about what the literature does and does not establish for it. The aldehyde dehydrogenase 3 family detoxifies lipid-peroxidation aldehydes such as 4-HNE and contributes to redox homeostasis (Liu et al., 2025; Xanthis et al., 2023), and individual ALDH3 isoforms are not interchangeable: ALDH3A2 loss sensitises ovarian cancer cells to ferroptosis by depleting glutathione (Dong et al., 2023), and ALDH3B2 overexpression drives proliferation and metastasis through ITGB1 signalling in cholangiocarcinoma (Wang et al., 2021). But the direct biochemical characterisation of ALDH3B2 enzymatic activity, substrate preference, and lipid-droplet localisation comes from mouse protein (Kitamura et al., 2015); the human gene carries a premature termination codon, and the only direct evidence for a full-length human protein is a placental study proposing stop-codon readthrough as the mechanism (Michorowska et al., 2019), which stands as a direct rebuttal to a separate report characterising human ALDH3B2 as a pseudogene rather than as settled confirmation of its own. The ovarian-cancer-specific evidence is transcriptomic and survival-association only: ALDH3B2 was among genes flagged as ovarian-carcinoma-specific in early expression profiling but was not one of the seven genes carried to immunohistochemical validation (Hibbs et al., 2004), and its significant univariable survival association (HR 0.983, 95% CI 0.970-0.996) did not survive into the final eight-gene prognostic signature it was tested alongside (Bi et al., 2021). No published inhibitor exists for ALDH3B2 specifically, unlike ALDH3A1, which has a characterised submicromolar selective inhibitor (Parajuli et al., 2014), so there is no positive control to benchmark docking against; ROOTS addresses this with a property-matched decoy benchmark against ALDH1A1 and ALDH3A1 inhibitors, documented as in progress in the Project History. And the structural model used throughout this pipeline is a predicted one, built via OpenFold3 against the canonical UniProt sequence, not an experimentally resolved ALDH3B2 structure; Michorowska et al. (2019) show that homology modelling against the rat class 3 ALDH template (50% sequence identity) is the available alternative, not a settled replacement.

None of this disqualifies ALDH3B2. It is exactly the profile that makes exploratory, well-validated computational screening worth doing before committing wet-lab resources: a biologically plausible, context-dependent cancer gene with real but incomplete human evidence, no existing chemical starting point, and a structure that has to be predicted rather than looked up. ROOTS treats each of those gaps as a validation requirement rather than an inconvenience to write around, which is the thread running through the re-docking RMSD check, the dual-layer AlphaFold confidence gating, and the consensus-versus-selective hit split described below.

## Pipeline

The pipeline runs in three phases, each depending entirely on the output of the one before it.

![ROOTS three-phase pipeline](figures/ROOTS_Codeathon_Aug2026_Part1.png)

**Phase 1, target and structure discovery:** TCGA-OV + GTEx + GEO comparative omics, target prioritisation, structural characterisation, COCONUT library curation, receptor and ligand preparation.

**Phase 2, screening and selection:** docking protocol validation, high-throughput docking, binding affinity estimation, hit definition and reconciliation, ADMET filtering.

**Phase 3, simulation, refinement, and output:** MD system preparation, MD production, binding free energy refinement, trajectory analysis, comparative figures, manuscript package assembly.

## Data Sources

| Layer | Source | Role |
|---|---|---|
| Discovery cohort | TCGA-OV vs. GTEx-ovary, via recount3 | Primary differential expression |
| Independent concordance | GEO GSE26712 (185 tumour, 10 normal, n=195) | Cross-platform confirmation of the discovery signal |
| Mutation frequency | GDC, cBioPortal | Orthogonal target support |
| Protein level | Human Protein Atlas | Immunohistochemistry cross-check |
| Compound library | COCONUT natural-product database | Sole screening library, archived at retrieval |
| Receptor structure | PDB (gated on resolution/coverage/organism), OpenFold3 as fallback | ALDH3B2 holo-NAD structural model |

## Methods

**Bioinformatics.** Target prioritisation through DESeq2/apeglm with ComBat-seq mitigation, independent GEO concordance, GDC/cBioPortal mutation frequency, an HPA protein-level check, and a fallopian-tube tissue-of-origin screen. Structural characterisation via PDB-first, OpenFold3-fallback, gated by TM-score and pocket-level pLDDT. Library curation through Ro5, Veber, and PAINS filtering plus a stereochemistry-integrity check confirming defined stereocentres survive 3D embedding.

**Cheminformatics.** A per-target re-docking RMSD check precedes any production docking. 280,509 compounds were docked with QuickVina2 across a 5-node Hetzner CX33 fleet (41,529 failed, overwhelmingly macrocyclic structures or unsupported B/Si heteroatoms). Every successfully docked pose was scored twice, once by Vina and once by PRODIGY-LIG, cross-validated against each other rather than averaged. Hits were defined by percentile thresholding, then split into consensus and promiscuity-flagged tiers. ADMET-AI triaged the hit set locally with no submission cap before a smaller shortlist went to SwissADME.

**Molecular dynamics.** Four complex-only systems were built in CHARMM-GUI Solution Builder against the ALDH3B2 holo-NAD receptor, covering both consensus and selective/flagged hits: Sys1 (CNP0220660.1, consensus), Sys2 (CNP0358193.1, consensus), Sys3 (CNP0384182.2, selective/flagged), Sys4 (CNP0088713.1, selective/borderline). Production ran on GROMACS 2026.3 (CUDA-enabled) on rented GPU infrastructure after the ~148,000-atom solvated system size made Colab's free tier impractical for unattended multi-hour runs. Binding free energy refinement (gmx_MMPBSA, Generalized Born and Poisson-Boltzmann) and trajectory analysis (MDAnalysis, PLIP, per-residue decomposition, representative frames) are both complete; see Results below for what that output does and does not yet establish.

## Results

**Docking and hit funnel**

| Stage | Count |
|---|---|
| Compounds docked (QuickVina2) | 280,509 (41,529 failed, mostly macrocycles) |
| Poses cross-scored (Vina + PRODIGY-LIG) | 280,507 |
| Hits (top 1 percentile, Vina) | ~2,805 |
| ADMET-AI triaged locally, no cap | 3,587 |
| Survived to SwissADME shortlist | 280 |
| High-confidence consensus hits | 26 |
| Promiscuity-flagged | 43 |

Vina and PRODIGY-LIG agreed weakly (Spearman rho approximately -0.22) across all 280,507 scored poses, against very different score ranges (12.7 kcal/mol for Vina, 1.1 kcal/mol for PRODIGY-LIG). That disagreement is reported rather than reconciled; Vina alone gates hit status.

**Molecular dynamics and MM-PBSA**

All four systems reached a converged window (approximately 34.4-49.4 ns) before production stopped under the project's compute-budget guardrails. That window, not an arbitrary fixed length, is what was carried into analysis:

| System | Compound | Hit tier | Frames analysed | ΔG complex, GB (kcal/mol) | ΔG complex, PB (kcal/mol) |
|---|---|---|---|---|---|
| Sys1 | CNP0220660.1 | Consensus | 157 | -2822.02 ± 4.43 (SEM) | -2847.35 ± 4.39 (SEM) |
| Sys2 | CNP0358193.1 | Consensus | 122 | -2601.36 ± 5.12 (SEM) | -2638.37 ± 4.91 (SEM) |
| Sys3 | CNP0384182.2 | Selective/flagged | 146 | -2962.78 ± 4.97 (SEM) | -2961.36 ± 5.12 (SEM) |
| Sys4 | CNP0088713.1 | Selective/borderline | 113 | -2656.62 ± 5.95 (SEM) | -2668.26 ± 5.95 (SEM) |

These are total complex energies (BOND+ANGLE+...+GSOLV) as reported by gmx_MMPBSA, not ligand-receptor binding deltas, and they are not quoted here as binding affinities. A proper ΔG_bind = G_complex - G_receptor - G_ligand from the per-component COM/REC/LIG decomposition is the next step before any binding-affinity claim is made. Per-residue decomposition, RMSD/RMSF, PLIP interaction profiles, and representative frames are complete for all four systems, and a six-panel cross-target comparison figure has been assembled from them.

## What Sets ROOTS Apart

- Target selection is omics-derived across multiple candidate lines rather than fixed by precedent before any analysis.
- The comparator design states the TCGA-vs-GTEx cohort confound plainly instead of batch-correcting over it, and reinforces target choice with independent GEO concordance, a fallopian-tube tissue-of-origin check, and HPA protein confirmation.
- Docking protocol validity is tested with a re-docking RMSD check before any production screening run, not assumed.
- Two independent binding-affinity methods are cross-validated against each other and the disagreement is disclosed, not smoothed into a single reconciled score.
- Hits are split into consensus and promiscuity-flagged tiers rather than reported as one undifferentiated list.
- The full compute stack, including where the original zero-cost design held and where real system scale forced a paid-infrastructure pivot, is documented rather than presented as solved in advance.

Full detail: [Project History & Decisions](https://github.com/omicscodeathon/roots/blob/main/docs/PROJECT_HISTORY_AND_DECISIONS.md).

## Reproducibility

```bash
git lfs install
git clone https://github.com/omicscodeathon/roots.git
cd roots
git lfs pull
mamba env create -f environments/colab_environment.yml  # or kaggle_environment.yml
```

Every MD system ships with a standalone `scripts/REPRODUCE_sys1.sh` through `REPRODUCE_sys4.sh` that rebuilds it from the committed CHARMM-GUI stage archives. Raw database pulls are archived at retrieval, random seeds are fixed wherever applicable, and the final package is planned for Zenodo DOI archiving alongside the GitHub release.

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
|   |-- step1/     (target shortlist, DE results, provenance -- corrected final versions)
|   |-- step2/     (validated receptor structure metadata)
|   |-- step3/     (final compound library)
|   |-- step5/     (merged PRODIGY-LIG cross-validation scores)
|   |-- step6/     (extracted docking scores, failure log -- not the 280k raw PDBQTs)
|   `-- step7/     (ADMET triage results, hit shortlist, SMILES)
|-- docs/
|   |-- PROJECT_HISTORY_AND_DECISIONS.md
|   |-- METHODS.md
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

Heavy, raw, or superseded artifacts (the 280k+ docking PDBQT files, raw AlphaFold3/OpenFold3 job output, raw per-node PRODIGY and docking fleet tarballs, Python virtualenvs, duplicate staging copies) are intentionally excluded from the repository; everything excluded is either regeneratable from the committed scripts or was a working intermediate never meant to be a final artifact. Full rationale: Project History & Decisions.

## Team

| Name | Affiliation(s) | Role(s) |
|---|---|---|
| John Adeyemo Adedeji | Osun State University Teaching Hospital | Team Leader |
| Tammy Lily Shishodiya | — | Team Co-Lead |
| Elsayed Hussien | — | Team Member |
| Onifade Olayinka | — | Team Member |
| Olaitan I. Awe | Institute for Genomic Medicine Research (IGMR); African Society for Bioinformatics and Computational Biology (ASBCB) | Project Advisor |

## Acknowledgments

We thank the maintainers of TCGA, GTEx, GEO, COCONUT, and the Human Protein Atlas for open biomedical data, and the GROMACS, AutoDock Vina, and gmx_MMPBSA communities for open-source tooling. We also thank Dr. Awe for ongoing guidance throughout target selection, infrastructure pivots, and simulation design.

## References

Full citation list and methodology sourcing: `docs/PROJECT_HISTORY_AND_DECISIONS.md` and `docs/METHODS.md`.
