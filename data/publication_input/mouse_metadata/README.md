# NASA OSDR mouse metadata input for Fig. 6

These tables were rebuilt on 2026-08-30 directly from the current NASA OSDR
dataset API and each accession's ISA metadata ZIP. The old local sample-metadata
cache is not used by the figure.

- `00_analysis_sample_scope.csv` freezes the 761 expression sample columns that
  contributed to the selected seven-pathway analysis. It defines analysis scope,
  not biological metadata.
- `01_osdr_isa_sample_source_core_all.csv` contains the core source/sample fields
  from all 1,588 ISA rows across the 48 accessions.
- `02_osdr_download_manifest.csv` records the OSDR page, API and ZIP URLs, the
  selected ZIP filename, download time, byte size and SHA-256 for every accession.
- `03_mouse_level_metadata.csv` is the one-to-one join of the 761 analysis sample
  columns to the freshly downloaded ISA rows and is the direct input to Fig. 6.
- `04_sample_join_audit.csv` reports join coverage, missingness and conflicts per
  accession.
- `05_selected_source_value_audit.csv` lists the exact source age, sex and mission
  values retained for the selected samples.
- `06_metadata_conflict_audit.csv` preserves source conflicts instead of hiding
  them. In particular, OSD-162's legacy ISA rows say `SpaceX-3`, whereas the
  current dataset API and OSDR study page identify the mission as `SpaceX-8`.

`R/data/01_refresh_osdr_mouse_metadata.R` repeats the download and rebuild. By
default it downloads new ZIPs into the ignored `data/.cache/osdr_isa/` directory.
Use `--reuse-cache` only for an offline rerun after verifying the cached hashes.

The plot counts unique `accession::Source Name` values, so technical replicates
do not inflate point size and short source labels from unrelated datasets cannot
collide. Both selected Flight and Ground samples are included and separately
counted in `results/tables/Fig_6_mouse_metadata_grouped.csv`.

Age labels are not inferred. Exact values remain exact (`32 {week}` becomes
`32 wk`); ranges remain ranges only when the downloaded NASA ISA cell itself is
a range (for example `16-17 {week}`). `age_order_weeks` is a display-order value,
not a replacement for the source age.
