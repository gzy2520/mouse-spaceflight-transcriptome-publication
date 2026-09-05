"""Append API-reviewed post-landing tissue timing to both 761-row mouse tables.

The publication-facing metric is when tissue was obtained after landing.  This
script reads each frozen source table, appends three traceable tissue-timing
columns (joined on ``accession``), and writes new files into this audit directory.
The source tables are never modified.  The strict euthanasia-only fields remain
available in the accession-level audit as provenance, but are not copied into the
publication-facing augmented tables.

Appended columns (identical in both outputs):
  tissue_collection_after_landing_days : blank, ``<1``, a day/range, or
                                         group-qualified values such as
                                         ``Space Flight=1`` or ``LAR=24``
  tissue_collection_timing_category    : how the value was established
  tissue_collection_scope              : cohort/tissue scope of the value

Exact evidence and raw API fields stay in
``post_landing_tissue_collection_audit_20260902.csv``; they are deliberately not
repeated 761 times here.  The accession column is the join key back to that audit.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
MM = HERE.parent                      # .../mouse_metadata
RELEASE_ROOT = HERE.parents[3]        # .../publication_release

AUDIT_CSV = HERE / "post_landing_tissue_collection_audit_20260902.csv"

# (source table, new output name)
TARGETS = [
    (MM / "03_mouse_level_metadata.csv",
     HERE / "09_mouse_level_metadata_with_tissue_collection_20260902.csv"),
    (RELEASE_ROOT / "results" / "tables"
     / "Fig_6_mouse_sample_metadata_complete_audit.csv",
     HERE / "10_fig6_mouse_sample_metadata_audit_with_tissue_collection_20260902.csv"),
]

NEW_COLS = [
    "tissue_collection_after_landing_days",
    "tissue_collection_timing_category",
    "tissue_collection_scope",
]


def load_audit() -> dict[str, dict]:
    with AUDIT_CSV.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    accessions = [r["accession"] for r in rows]
    duplicates = sorted({acc for acc in accessions if accessions.count(acc) > 1})
    if duplicates:
        raise SystemExit(f"audit has duplicate accessions: {duplicates}")
    if not rows:
        raise SystemExit(f"audit is empty: {AUDIT_CSV}")
    required = {
        "tissue_collection_after_landing_days",
        "tissue_collection_timing_category",
        "tissue_collection_scope",
    }
    missing = sorted(required - set(rows[0]))
    if missing:
        raise SystemExit(f"audit missing tissue columns: {missing}")
    return {r["accession"]: r for r in rows}


def process(src: Path, dst: Path, audit: dict[str, dict]) -> None:
    with src.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if "accession" not in fieldnames:
        raise SystemExit(f"{src.name}: no 'accession' column to join on")
    if any(c in fieldnames for c in NEW_COLS):
        raise SystemExit(f"{src.name}: already has a target column; refusing")

    unmapped: Counter = Counter()
    source_accessions = {row["accession"] for row in rows}
    missing_audit = sorted(source_accessions - set(audit))
    if missing_audit:
        raise SystemExit(f"{src.name}: audit missing accessions: {missing_audit}")
    for row in rows:
        acc = row["accession"]
        a = audit.get(acc)
        if a is None:
            unmapped[acc] += 1
            row.update({c: "" for c in NEW_COLS})
        else:
            for col in NEW_COLS:
                row[col] = a[col]

    with dst.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames + NEW_COLS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"{src.name} -> {dst.name}")
    print(f"   rows={len(rows)}  out_cols={len(fieldnames) + len(NEW_COLS)}")
    if unmapped:
        raise SystemExit(f"{src.name}: unexpected unmapped accessions: {dict(unmapped)}")
    # distribution in this table
    dist = Counter(r["tissue_collection_timing_category"] for r in rows)
    print(f"   category distribution: {dict(dist)}")


def main() -> int:
    audit = load_audit()
    print(f"audit accessions: {len(audit)}")
    for src, dst in TARGETS:
        if not src.exists():
            raise SystemExit(f"missing source table: {src}")
        process(src, dst, audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
