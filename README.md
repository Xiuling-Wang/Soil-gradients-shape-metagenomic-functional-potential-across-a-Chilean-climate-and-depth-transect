# Chilean soil metagenomics: resource-use potential

Code for “Regional differentiation and genomic context of microbial resource-use potential along the Chilean Coastal Cordillera.” The study was submitted to Environmental Microbiome on 19 October 2026; V16 is the current revised source manuscript, with no acceptance asserted.

## Current analysis and data

The [5 November 2026 release](release_2026-11-05/README.md) adds the current 15-library analysis code. Original C03 and C18 are excluded without remapping. C01 is descriptive; SG/LC/NA surface inference uses n=9, NA depth uses n=8, and complete-pit sensitivity uses n=6. The fixed gene and 228-MAG references are unchanged.

- Current analysis entry points, input requirements and limits: [release guide](release_2026-11-05/README.md).
- Current supplementary workbooks and figures: [Zenodo 1.0.6](https://doi.org/10.5281/zenodo.23958664). Source-paper V16; Data C1 has 18 sheets.
- Raw reads, all 17 deposited runs: [ENA PRJEB123235](https://www.ebi.ac.uk/ena/browser/view/PRJEB123235). Deposited C03 and C18 are excluded from current comparisons.
- Current code/source checksums: [manifest](release_2026-11-05/CODE_MANIFEST.json).

## Historical code

The older [figure guide](figs.md), script/figures.py runner and root CODE_MANIFEST.json describe the earlier submitted snapshot. They remain available for provenance and are not the current 15-library analysis entry points. The current release guide distinguishes numerical analysis, released result workbooks and final figure reproduction.

## Reproduction scope

The repository distributes code, not the full gene-count/annotation input catalogue or every intermediate. A complete rerun requires the documented external input layout. Current scripts were syntax-checked; models were not rerun for this update. Sequence abundance describes encoded potential, and recovered MAG gene combinations do not establish transcription, activity or physiological function. See the release guide for exact sample sets and method boundaries.
