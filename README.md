# Soil Metagenomics Analysis Code

Code-only snapshot for the V35 manuscript. No measurements, sample metadata, sampling coordinates, annotations, abundance matrices, sequence files, manuscripts or internal reviews are included.

Machine-specific paths have been replaced with portable roots. Run R scripts from the project root, or set `CHILE_PROJECT_ROOT`; the phylogenetic Python script also accepts `CHILE_PHYLO_BIN`. Mathematical operations, sample selection and statistical models were not changed. The original script numbers preserve analytical provenance rather than manuscript numbering.

## Entry Points

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

- Surface comparisons use independent sites. Depth comparisons are paired within soil pits.
- Script 261 supplies the current exact-design tests. Unrestricted depth tests in earlier scripts are historical, not current inference.
- In script 332, the pit-adjusted model `~ pit + depth_display` is the reported depth comparison. Its older unblocked-depth labels and generated summaries are retained for source traceability, not promoted as current primary results.
- Gene, KO and MAG representations are complementary views of the same metagenomes. Sequence abundance measures encoded potential, not transcription or process rates.

## Inputs and Software

Input files are intentionally not distributed in this code-only repository. Their required relative paths are declared in each script. In particular, the exact-design and robustness workflows also require historical reference TSVs; the AmoA workflow requires MAG contigs and annotations. Code availability alone does not establish data availability or a complete reproducible release.

Required R packages include `data.table`, `tidyverse`, `vegan`, `permute`, `edgeR`, `rdacca.hp`, `patchwork`, `cowplot`, `ggrepel`, `ggtext`, `ComplexHeatmap`, `circlize`, `ggplotify` and `scales`. The Python workflow uses the standard library plus external MAFFT, trimAl and IQ-TREE executables. Live UniProt queries can return updated records.

The distributed code was syntax-checked and scanned for identifying strings, home-directory paths, email addresses and credential patterns. Statistical models were not rerun. This is not a complete end-to-end reproducibility certification, and no new reuse license is assigned by this snapshot.

`CODE_MANIFEST.json` records the hashes of the distributed scripts. Generated outputs and datasets are excluded from Git. Commit attribution uses an anonymous placeholder; the hosting account remains visible as the repository owner.
