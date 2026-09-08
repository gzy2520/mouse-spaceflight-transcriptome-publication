# Post-landing tissue-collection timing audit — API-first review (2026-09-02)

This directory contains the publication-facing timing metadata for the 48
analysis accessions and the two 761-row mouse metadata tables.  The final metric
is **how long after landing the analyzed tissue was obtained**.  The earlier
euthanasia-only interpretation is retained in the accession audit for
provenance, but it is not the field copied into the final augmented tables.

No frozen input or approved `results/` file was modified.  The two final tables
are new copies of the source tables with three tissue-timing columns appended.

## Primary source and rules

The saved JSON payloads in `repo_study_fresh_20260901/` are the primary source;
they are the 2026-09-01 snapshots fetched from the OSDR repository-study API.
`study_descriptions_20260901.csv` is only a description cross-check.  The audit
builder applies the following order:

1. Use explicit `samples.table` dissection-date fields when present and compute
   date minus the API `missionEnd` **for the flight group only**.  Ground, basal,
   vivarium, and cohort-control dates are retained as raw API provenance but do
   not have a landing interval.
2. Use an explicit API description statement about collection/dissection time
   when no usable sample date is available.
3. When the description gives a live-animal-return recovery before sacrifice or
   euthanasia but no separate dissection timestamp, report that interval as a
   `post_landing_euthanasia_proxy`.  This follows the requested convention that
   the terminal tissue procedure is effectively the animal endpoint, while
   keeping the proxy status visible.
4. Keep on-orbit tissue as `on_orbit` (no post-landing interval).  Mixed cohorts
   retain their scope.  Mixed or conflicting API date text keeps the raw value
   and leaves the numeric interval blank rather than forcing a number.

## Canonical outputs

| File | Role |
|---|---|
| `post_landing_tissue_collection_audit_20260902.csv` | **Canonical 48-row audit**; final tissue timing, raw API fields, exact description evidence, and the retained strict euthanasia fields |
| `09_mouse_level_metadata_with_tissue_collection_20260902.csv` | `03_mouse_level_metadata.csv` + 3 tissue-timing columns (761 rows) |
| `10_fig6_mouse_sample_metadata_audit_with_tissue_collection_20260902.csv` | Fig. 6 audit table + the same 3 tissue-timing columns (761 rows) |
| `build_euthanasia_audit_reviewed_20260902.py` | Rebuilds the canonical API-first accession audit |
| `build_mouse_tables_with_euthanasia_reviewed_20260902.py` | Joins the canonical tissue fields onto both frozen source tables (filename retained for backward compatibility) |

The files ending in `_first_pass_20260901` or the unsuffixed augmented-table
names are historical first-pass outputs and are not the final tissue metric.
`euthanasia_after_landing_audit_20260901_reviewed.csv` is the earlier strict
euthanasia review; it is preserved separately from the canonical tissue audit.

## Important reviewed cases

- **OSD-714**: the final tissue field is `2`; the API says skeletal muscle was
  collected two days after landing.  The strict euthanasia field remains blank
  because the API does not separately state euthanasia timing.
- **OSD-240/241, OSD-379, OSD-580 and OSD-771**: the API description gives a
  post-landing live-animal-return recovery before the terminal procedure, so the
  final tissue field records a clearly labelled euthanasia/sacrifice proxy
  (`30`, `LAR=2`, or `LAR=24`).  On-orbit subgroups are stated in `scope`.
- **OSD-239, OSD-242, OSD-576, OSD-665/666/770**: explicit API description
  collection/dissection statements are retained verbatim in the tissue evidence
  fields.  For RR-23, the final post-landing value is `FLT=1`; HGC/VGC control
  dates are not treated as landing intervals.
- **OSD-253/254**: API dissection dates are tied to mixed or non-comparable
  cohorts, so raw values are retained but no single post-landing number is
  assigned.
- Every evidence component in both the strict and tissue evidence columns is
  checked as a verbatim substring of that accession's fresh API description.

## Final appended columns

- `tissue_collection_after_landing_days` — blank, `<1`, a range, or a
  flight-qualified value such as `Space Flight=1` or `LAR=24`.
- `tissue_collection_timing_category` — the evidence class, including
  `post_landing_api_dissection_date_derived`,
  `post_landing_description_dissection_date_derived`,
  `post_landing_euthanasia_proxy`, `on_orbit`, and unresolved/mixed classes.
- `tissue_collection_scope` — cohort and tissue scope needed to interpret a
  group-qualified value.

The complete source field/value, exact evidence sentences, confidence, notes,
and the strict `euthanasia_*` fields are in the canonical accession audit.

## Reproduce the final files

Run from the `publication_release/` root:

```bash
python3 data/publication_input/mouse_metadata/euthanasia_audit_20260901/build_euthanasia_audit_reviewed_20260902.py
python3 data/publication_input/mouse_metadata/euthanasia_audit_20260901/build_mouse_tables_with_euthanasia_reviewed_20260902.py
```

The builders read the saved API JSON snapshots and the frozen source tables;
they verify the saved JSON hashes and require `fetch_status=ok` in the manifest.
The live fetch script is only needed when a new API snapshot is intentionally
requested.
