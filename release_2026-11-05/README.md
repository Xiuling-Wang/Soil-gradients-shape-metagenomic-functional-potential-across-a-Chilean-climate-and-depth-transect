# Current metagenomic analysis code

This release accompanies source-paper V16 and the accepted thesis V17. The current analysis excludes original C03 and C18 because their identities were unresolved, retains C01 and C04–C17, and performs no remapping. The fixed gene catalogue and 228-MAG reference are retained, including genomes originally assembled from excluded libraries as reference members only. The raw ENA archive still contains 17 runs.

## Current entry points

Run from a complete project root with the input layout in INPUT_MANIFEST.json. The repository itself is not a complete data bundle. Use a separate writable project copy for a rerun.

1. script/c03_c18_exclude_both_matrices_20261028.R — retained analysis sets, full-catalog library sizes, TMM factors and KO/MAG matrices.
2. script/c03_c18_exclude_both_curated_20261028.R — frozen annotation-name filters and within-module fractional allocation. Uses source_code/38_curated_module_validation.R in this release; only its quantification expressions are evaluated.
3. script/c03_c18_exclude_both_statistics_20261028.R — composition, restricted permutations, soil-axis/missing-data sensitivity and category/module summaries.
4. script/c03_c18_exclude_both_edgeR_20261028.R — full-family KO omnibus and pairwise tests, robust pit-adjusted depth models and the complete-pit sensitivity.
5. script/c03_c18_exclude_both_outputs_20261028.R — adopted result export and regression checks. This provenance/export stage additionally requires the historical comparison outputs named in the script.
6. script/815_v55_leave_C01_out_DIAGNOSTIC.R — separately labelled leave-C01-out soil-axis diagnostic. It does not replace the primary soil-PCA analysis.

The replicated surface comparison is SG/LC/NA n=9. AZ C01 is descriptive. Primary NA depth n=8; complete-pit sensitivity n=6. No climate-class replication is inferred from one site per climate class. Curated-marker and whole-KO representations have distinct assignment rules. The depth KO model includes pit; depth distance tests restrict permutations within pit.

## Inputs and result files

The retained-library metadata and current results are in the V16 workbooks in [Zenodo 1.0.6](https://doi.org/10.5281/zenodo.23958664). Tables C1–C4 retain source S sheet names; Data C1 has 18 sheets, including NA_marker_pairwise, Selected_MAG_neighbourhoods and Analysis_notes. The latter records model definitions and sources. These workbooks are the released result tables, not every input matrix.

INPUT_MANIFEST.json records required frozen input paths/hashes without distributing them. The scripts also reference the fixed gene catalogue, MAG annotations, marker registry and historical diagnostic outputs; inspect each entry point before execution. For later selected-MAG neighbourhood annotations and display construction, the current workbooks provide the released evidence; this code release does not claim to reproduce all final figures or automatically rebuild every workbook sheet.

## Changes and validation

The distributed R files were copied from the current project and syntax-checked. Only project-root paths and the frozen-source location were adapted. Line counts and scientific expressions were preserved. No model was rerun for this publication update. CODE_MANIFEST.json records both source and distributed hashes.

Earlier numbered scripts and renderers remain historical provenance. In particular, scripts 261, 281, 332 and 420 and the old figure runner predate the C03/C18 exclusion and do not define the current primary analysis. The 23 October release documented the earlier submission. No new software/data licence is assigned by this code synchronization.
