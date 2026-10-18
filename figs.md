# Figures

Use `python script/figures.py --inputs INPUTS --work NEW_WORKSPACE`.
Add `--check` to check required input paths without writing files. The output
directory must not already exist; source inputs are linked only for reading.

| Figure | Approved rendering source |
|---|---|
| 1 | Original source PDF plus square-root scale label |
| 2 | Original vector panels with revised spacing |
| 3, 4 | `760` wrapping `465`: lighter category strips and refined heatmap geometry |
| 5 | Original A plus `766` wrapping `553`/`467`: original layout, 1.25x bubble diameters |
| 6 | `760`/`695`/`470`: corrected split support, lighter shading and panel spacing |
| 7 | `767`: graphical synthesis, italic statistical p and explicit gene combinations |
| S1 | Original source PDF plus revised C from `760`/`446` |
| S2 | `766` using `212`/`471`: complete native A/B pages, original large legends |
| S3 | Unchanged source PDF |

The runner is the current entry point. Do not execute old renderers independently
to reconstruct this selection: their historical default outputs differ. Extracted
sections of `467` and `553` are used for Figure 5, not their whole legacy pipelines.

## Local Input Tree

The input root must retain these relative paths (also checked by the runner):

```text
data/metadata/env_2024_meta.tsv
data/processed/metagenome_output/MAGs/
Environmental_Microbiome/Coauthor_Review_2026-10-01_XW/Main_Figures_PDF/
Environmental_Microbiome/Coauthor_Review_2026-10-01_XW/Supplementary_Figures_PDF/
Environmental_Microbiome/Coauthor_Review_2026-10-01_XW/Supplementary_Data/
Environmental_Microbiome/V35/Supplementary_Data/AmoA_sequences/
analysis/V42_ecological_story_revision_2026-08-26/source_data/
analysis/V53_author_comments_revision_2026-08-28/source_data/
analysis/V59_KO_ordination_architecture_2026-08-29/source_data/
analysis/V05_comment_review_2026-08-17/figure2_reanalysis/
analysis/exploratory_KO_differential_abundance_2026-08-29/KO_DA_README.md
```

These are provenance paths, not downloads supplied by this code-only release.
Specific file names and column checks remain in the source renderers. Incomplete
input trees fail with a missing-input message. No model or sequence workflow is
launched automatically to fill missing inputs.

The resulting ten complete PDFs are under
`NEW_WORKSPACE/Environmental_Microbiome/Figure_table_review_2026-10-18_round2/`.
Intermediate exports and invariant checks stay inside the new workspace.
