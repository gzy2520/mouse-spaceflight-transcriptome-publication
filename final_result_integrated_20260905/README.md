# Integrated final-result refresh (2026-09-05)

`final_result_integrated_20260905/` is a reproducible, dated integration of
the approved `final_result/` snapshot and the three later, source-rendered
updates below.  The approved snapshot is never modified.

## Integrated updates

- **Fig. 1b** is re-rendered with enlarged legend typography and point-size
  keys.  Its 15 discrete age colours match the supplied final-reference image:
  chronological yellow -> green -> blue -> purple.
- **Fig. 2 and Fig. S6** use mouse/sample-level, ground-referenced log2FC;
  within-tissue Spearman correlations are aggregated across eligible tissues
  by equal-weight Fisher-z.  The compact overall Figure 2 uses the same
  compact GO names on both axes; per-tissue and merged Fig. S6 panels are
  included.
- **Fig. 4a and Fig. 5a** retain the frozen YARN qsmooth matrix, Spearman tree
  distances, tissue order and palette, but omit numeric cell labels.  Their
  right-side tissue labels are English.

All other main figures, supplementary figures, provenance files and 27 tables
are copied unchanged from `final_result/`.

## Reproduce

From the `publication_release` repository root:

```sh
PYTHON=.venv/bin/python bash workflow/build_integrated_final_result_20260905.sh final_result_integrated_20260905
```

The builder refuses to overwrite an existing output directory.  It renders
the three updates from their source scripts, verifies the original table set
and qsmooth tree audit, checks that all 26 per-tissue PNG/PDF pairs exist, and
writes SHA-256 records for the integrated replacements.

`PYTHON` must point to an interpreter with the packages in `requirements.txt`.
After the repository setup in `README.md`, use `.venv/bin/python`; an existing
environment may be supplied instead.
