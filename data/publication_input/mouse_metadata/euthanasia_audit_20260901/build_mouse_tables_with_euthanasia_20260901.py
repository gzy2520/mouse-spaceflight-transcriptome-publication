"""Append the post-landing euthanasia columns to BOTH 761-row mouse tables.

The user asked to add the post-landing euthanasia timing to the mouse info tables
WITHOUT overwriting the frozen / approved originals.  This script therefore reads
each source table, appends three traceable columns (joined on ``accession``), and
writes a NEW file into this audit directory.  The source files are never modified.

Appended columns (identical in both outputs):
  euthanasia_after_landing_days  : "" or "<1" / "1" / "1-2" / "2" / "24" / "30"
  euthanasia_timing_category     : post_landing | post_landing_date_derived
                                   | on_orbit | not_stated
  euthanasia_scope               : which cohort the value applies to (see audit)

The exact evidence sentence(s) for each accession stay in
``euthanasia_after_landing_audit_20260901_reviewed.csv`` and are deliberately not
repeated 761 times here; the accession column is the join key back to that audit.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
MM = HERE.parent                      # .../mouse_metadata
RELEASE_ROOT = HERE.parents[3]        # .../publication_release

AUDIT_CSV = HERE / "euthanasia_after_landing_audit_20260901_reviewed.csv"

# (source table, new output name)
TARGETS = [
    (MM / "03_mouse_level_metadata.csv",
     HERE / "09_mouse_level_metadata_with_euthanasia_20260901_reviewed.csv"),
    (RELEASE_ROOT / "results" / "tables"
     / "Fig_6_mouse_sample_metadata_complete_audit.csv",
     HERE / "10_fig6_mouse_sample_metadata_audit_with_euthanasia_20260901_reviewed.csv"),
]

NEW_COLS = ["euthanasia_after_landing_days", "euthanasia_timing_category",
            "euthanasia_scope"]


def load_audit() -> dict[str, dict]:
    with AUDIT_CSV.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    accessions = [r["accession"] for r in rows]
    duplicates = sorted({acc for acc in accessions if accessions.count(acc) > 1})
    if duplicates:
        raise SystemExit(f"audit has duplicate accessions: {duplicates}")
    if not rows:
        raise SystemExit(f"audit is empty: {AUDIT_CSV}")
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
            row["euthanasia_after_landing_days"] = a["euthanasia_after_landing_days"]
            row["euthanasia_timing_category"] = a["category"]
            row["euthanasia_scope"] = a["euthanasia_group_scope"]

    with dst.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames + NEW_COLS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"{src.name} -> {dst.name}")
    print(f"   rows={len(rows)}  out_cols={len(fieldnames) + len(NEW_COLS)}")
    if unmapped:
        raise SystemExit(f"{src.name}: unexpected unmapped accessions: {dict(unmapped)}")
    # distribution in this table
    dist = Counter(r["euthanasia_timing_category"] for r in rows)
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
