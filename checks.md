# Release Checks

## 23 October 2026 public-link and provenance synchronization

- The 27 scientific source files remain byte-identical to the existing code snapshot; no models or sequence analyses were rerun.
- Current submitted title, 19 October submission date, ENA accession and Zenodo 1.0.5 links were synchronized.
- Final submitted workbooks and four finalized figure PDFs are fingerprinted separately in `CODE_MANIFEST.json`; they are not uploaded to this code-only repository.
- The following renderer certification belongs to the earlier approved selection. It does **not** certify the subsequent Figure 6 finishing and journal-specific taxonomic typography changes.

## Original code-snapshot checks

Snapshot: `v49`, 19 October 2026. Scope: code only.

- All 27 Python/R source files passed syntax checks.
- The 13 files from the original analysis snapshot are byte-identical.
- The current figure runner completed in a new workspace using local frozen inputs.
- All 10 complete figure PDFs (11 pages) match the approved selection at 216 dpi,
  with zero differing rendered pixels and identical extracted page text.
- Both table-workbook outputs have identical ZIP entry names and uncompressed
  entry bytes to the approved workbooks. Values, formulas, formats, column widths,
  row heights, filters, freeze panes and table styles match (170,289 cells checked).
- Manuscripts and approved source outputs remained unchanged.
- Code was checked for home-directory paths, email addresses and token patterns.
- No raw data, processed data, source PDFs, workbooks, manuscripts or internal
  review records are included in the published tree.

Rendering was tested on macOS using the existing Arial/R environment. Tests
required separately held inputs and do not establish a public, data-complete
reproducible archive. Statistical models and sequence analyses were not rerun.
No software versions are inferred from historical dates.

`CODE_MANIFEST.json` binds each distributed source to its SHA-256 and retains
the original approved output hashes as historical provenance. PDF binary hashes may differ after rerendering
because of creation metadata; the comparison above uses rendered pages and text.
