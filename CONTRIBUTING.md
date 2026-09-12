# Contributing to ROOTS

This file states the working standard for this project. Read it before your first pull request.

## Verification standard

Every tool, API, accession, and dataset used in this pipeline is verified against its primary source before being committed to code, documentation, or results. This is not a suggestion, it is the standard the whole design has been held to so far, and it is the standard every contribution is expected to meet.

Concretely:

- **Check tool identity before use.** Tools with similar names can do genuinely different things. PRODIGY (protein-protein binding affinity) and PRODIGY-LIG (protein-ligand binding affinity) are not interchangeable, and this distinction was caught and corrected during the pipeline's design. Confirm you are using the tool that matches the interaction type in question, not the one that matches by name alone.
- **Check API and platform behaviour against current documentation, not assumption.** An inherited assumption that Google Colab could not reliably run molecular dynamics was overturned by a published, benchmarked checkpoint-persistence method, found only by checking rather than accepting the earlier limitation as fact. Treat any stated infrastructure or tool limitation as something to verify, especially if it is load-bearing for a design decision.
- **Check accessions, database versions, and reference releases explicitly.** Databases update without versioning (COCONUT, GDC, GTEx, Ensembl/GENCODE). Record the exact release or retrieval date used, and archive raw pulls before filtering, so results can be reproduced or audited later.
- **Navigate to primary sources directly.** Do not rely on secondhand summaries, forwarded links, or search snippets for anything that will appear in code, documentation, or a result. Go to the publisher, the database, or the tool's own documentation.

## When you find an error

Document it, do not quietly absorb it. If a design decision, a tool choice, or a prior result turns out to be wrong, correct it openly: state what was assumed, what the check found, and what changed as a result. Several real corrections were made this way during the pipeline's design (a transcriptomic comparator confound, a tool-identity error, an infrastructure assumption), and each one made the pipeline stronger, not weaker, for having been caught and shown rather than fixed silently.

Null or negative results are reported honestly, not hidden. A statistically non-significant result, a docking pose that fails to converge, or a target that does not survive validation is a real finding and belongs in the record.

## Before opening a pull request

- Confirm any tool, accession, or dataset reference against its primary source.
- If your change touches the omics layer (target prioritisation, differential expression, comparator design), flag it in the pinned "Statistical review: open concerns" thread in `#roots-omics` before merging, since four other pipeline stages build on this layer's output.
- State clearly in the PR description what was checked and how, not just what changed.

## Scope boundaries

Some design choices in this pipeline are deliberate scope boundaries, not oversights: CNV exclusion from target prioritisation, single structural state per target, no enhanced sampling or metadynamics in this pass. If you believe one of these should change, raise it as an issue with your reasoning rather than changing it silently in a PR, since these boundaries were set for stated reasons documented in the README.
