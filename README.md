# Chilean soil metagenomics: resource-use potential

Analysis and figure code for a Chilean soil climate transect and paired depth profiles. No measurements, sample metadata, sampling coordinates, annotations, abundance matrices, sequence files, manuscripts or internal reviews are included.

Current submitted title: **Regional differentiation and genomic context of microbial resource-use potential along the Chilean Coastal Cordillera.** Submitted to Environmental Microbiome on 19 October 2026; no acceptance is asserted.

Raw metagenomic reads (17 runs): [ENA PRJEB123235](https://www.ebi.ac.uk/ena/browser/view/PRJEB123235). Current result workbooks and thesis supplementary figures: [Zenodo 1.0.5](https://doi.org/10.5281/zenodo.23225015). Zenodo maps thesis Appendix C labels to the journal S labels retained inside the submitted workbooks.

The 23 October synchronization updates the study identity, public links, input boundaries and submitted-file provenance. The 27-file scientific-code snapshot is unchanged. The figure runner reproduces the earlier approved selection, not the final 19 October typography and Figure 6 finishing steps; final submitted-output hashes are recorded separately and are not a claim of full end-to-end reproduction.

## Figures

[Figure guide](figs.md) | [Runner](script/figures.py)

```sh
python script/figures.py --inputs /path/to/local-inputs --work /path/to/new-workspace
```

Builds the earlier approved seven-main/three-supplementary figure selection, including the larger-bubble Figure 5, graphical Figure 7 and complete two-page S2. Frozen results are used without model refitting. Unchanged panels come from source PDFs. These PDFs and underlying data must be supplied separately. Existing inputs are not overwritten. See the synchronization note above for the final submitted typography boundary.

## Tables

[Formatting code](script/tables.py)

```sh
python script/tables.py --inputs /path/to/original-workbooks --out /path/to/new-tables
```

Applies S1A widths and displayed precision, S4 units in headers, and the supplementary-data table style. Measurements, formulas and other sheets are retained. Input names: `Tables_S1-S4.xlsx` and `Supplementary_Data_S1.xlsx`.

## Analysis

Machine-specific paths have been replaced with portable roots. Run R scripts from the project root, or set `CHILE_PROJECT_ROOT`; the phylogenetic Python script also accepts `CHILE_PHYLO_BIN`. Mathematical operations, sample selection and statistical models were not changed. The original script numbers preserve analytical provenance rather than manuscript numbering.

### Entry Points

| Component | Script |
|---|---|
| Curated marker registry and TPM construction | `38_curated_module_validation.R`, `_functional_marker_registry.R` |
| Soil axes and cross-resolution matrices | `148_rerun_figure2_comments_environmental_axes.R`, `150_rebuild_figure3_cross_resolution.R` |
| Gene/KO/MAG composition | `158_audit_v06_ecological_spine.R` |
| Exact within-pit permutations, missingness and module tests | `261_v40_phase1_exact_permutation_sensitivity.R` |
| rCLR and MAG quality sensitivity | `281_v43_targeted_robustness_analysis.R` |
| KO differential abundance | `333_prepare_exploratory_KO_DA_inputs.R`, `332_exploratory_KO_differential_abundance.R` |
| AZ influence diagnostic | `420_v19_surface_rclr_az_sensitivity.R` |
| Bacterial MAG gene combinations | `166_explore_bacterial_MAG_candidates.R` |
| AmoA sequence validation | `123_build_amoA_targeted_validation.py` |
| Shared figure theme | `_theme_sbb.R` |

## Current Versus Historical Analyses

- Surface comparisons use soil profiles as biological replicates within four observational sites; climate classes are not independently replicated. Depth comparisons are paired within NA soil pits.
- Script 261 supplies the current exact-design tests. Unrestricted depth tests in earlier scripts are historical, not current inference.
- In script 332, the pit-adjusted model `~ pit + depth_display` is the reported depth comparison. Its older unblocked-depth labels and generated summaries are retained for source traceability, not promoted as current primary results.
- Gene, KO and MAG representations are complementary views of the same metagenomes. Sequence abundance measures encoded potential, not transcription or process rates.

## Inputs and Software

Input files are intentionally not distributed in this code-only repository. Their required relative paths are declared in each script. In particular, the exact-design and robustness workflows also require historical reference TSVs; the AmoA workflow requires MAG contigs and annotations. Code availability alone does not establish data availability or a complete reproducible release.

Required R packages include `data.table`, `tidyverse`, `vegan`, `permute`, `edgeR`, `rdacca.hp`, `patchwork`, `cowplot`, `ggrepel`, `ggtext`, `ComplexHeatmap`, `circlize`, `ggplotify` and `scales`. The Python workflow uses the standard library plus external MAFFT, trimAl and IQ-TREE executables. Live UniProt queries can return updated records.

The current figure/table Python dependencies are in [requirements.txt](requirements.txt). Figure rendering also uses R packages `ape`, `ggtree`, `ggnewscale`, `gridtext` and `jsonlite`. It is tested on macOS with Arial. Set `CHILE_FONT_DIR` when `Arial.ttf`, `Arial Bold.ttf` and `Arial Italic.ttf` are outside the system font directory. The inherited Figure 6 renderer uses Quartz; cross-platform rendering is not certified. The figure runner configures its own isolated project root.

The distributed code was syntax-checked and scanned for identifying strings, home-directory paths, email addresses and credential patterns. Statistical models were not rerun. This is not a complete end-to-end reproducibility certification, and no new reuse license is assigned by this snapshot.

`CODE_MANIFEST.json` records distributed code hashes, historical rendered-output hashes and separately identified final submitted-file hashes. [checks.md](checks.md) distinguishes the original renderer tests from this metadata/provenance synchronization. Generated outputs and datasets are excluded from Git. Commit attribution uses an anonymous placeholder; the hosting account remains visible as the repository owner.
