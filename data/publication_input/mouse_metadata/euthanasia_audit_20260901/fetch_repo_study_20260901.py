"""Re-fetch the OSDR repository-study JSON for the 48 analysis accessions.

Purpose
-------
Produce a fresh, time-stamped, hash-audited copy of each accession's repository
study record (which carries the free-text ``description`` field) so the
post-landing euthanasia-timing audit can cite a current API snapshot rather than
the stale June-2026 workbook cache.

The script is transport-only in spirit: it reuses the shared ``osdr_api`` client
for the GET / retry / decode, then persists (1) one JSON file per accession and
(2) two audit CSVs. It does no scientific interpretation.

Run from the ``publication_release/`` repository root:

    PYTHONPATH=<main-repo>/07_scripts python3 \
        data/publication_input/mouse_metadata/euthanasia_audit_20260901/fetch_repo_study_20260901.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from shared_rigor.osdr_api import OSDRApiClient, OSDRApiError

HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parents[3]  # .../publication_release
MANIFEST = (
    RELEASE_ROOT
    / "data"
    / "publication_input"
    / "mouse_metadata"
    / "02_osdr_download_manifest.csv"
)
JUNE_CACHE_DIR = (
    RELEASE_ROOT
    / ".."
    / "01_osdr_database"
    / "osdr_catalog_20260601"
    / "resource_workbook_v2"
    / "repo_study_cache"
)
FRESH_DIR = HERE / "repo_study_fresh_20260901"
FETCH_MANIFEST_CSV = HERE / "fetch_manifest_20260901.csv"
DESCRIPTIONS_CSV = HERE / "study_descriptions_20260901.csv"


def _accessions() -> list[str]:
    rows = list(csv.DictReader(MANIFEST.open(newline="", encoding="utf-8")))
    accs = [r["accession"].strip() for r in rows]
    if len(accs) != len(set(accs)):
        raise SystemExit("Duplicate accessions in the download manifest")
    return accs


def _sha256(payload: dict) -> str:
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def main() -> int:
    FRESH_DIR.mkdir(parents=True, exist_ok=True)
    client = OSDRApiClient()
    accs = _accessions()
    print(f"Fetching {len(accs)} accessions from the live OSDR API ...")

    manifest_rows: list[dict] = []
    desc_rows: list[dict] = []
    failures: list[str] = []

    for i, acc in enumerate(accs, start=1):
        fetched_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        url = client.endpoints.repo_study.format(accession=acc)
        june_path = JUNE_CACHE_DIR / f"{acc}.json"
        june_sha = ""
        june_desc = ""
        if june_path.exists():
            june = json.loads(june_path.read_text(encoding="utf-8"))
            june_sha = _sha256(june)
            june_desc = june.get("description", "")

        status = "ok"
        payload: dict = {}
        try:
            payload = client.fetch_repo_study(acc, cache=False)
        except OSDRApiError as exc:
            status = f"fetch_failed: {exc}"
            failures.append(acc)
            payload = {}
            # Fall back to the June cache so a single network hiccup does not
            # block the audit; the status column records that this row is stale.
            if june_path.exists():
                payload = json.loads(june_path.read_text(encoding="utf-8"))
                status = "fetch_failed_used_june_cache"

        fresh_path = FRESH_DIR / f"{acc}.json"
        if payload:
            fresh_path.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )

        fresh_sha = _sha256(payload) if payload else ""
        desc = payload.get("description", "") if payload else ""
        manifest_rows.append(
            {
                "accession": acc,
                "repo_study_url": url,
                "fetched_at_utc": fetched_at,
                "fetch_status": status,
                "fresh_sha256": fresh_sha,
                "june_cache_sha256": june_sha,
                "description_matches_june": (
                    "" if (not payload or not june_path.exists())
                    else str(desc == june_desc)
                ),
                "study_version": payload.get("version", "") if payload else "",
                "modified_date": payload.get("modifiedDate", "") if payload else "",
            }
        )
        desc_rows.append(
            {
                "accession": acc,
                "title": payload.get("title", "") if payload else "",
                "flight_program": payload.get("flightProgram", "") if payload else "",
                "mission_name": payload.get("missionName", "") if payload else "",
                "doi": payload.get("doi", "") if payload else "",
                "description": desc,
            }
        )
        print(f"  [{i:2d}/{len(accs)}] {acc}  status={status}  desc_len={len(desc)}")
        time.sleep(0.3)

    with FETCH_MANIFEST_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)

    with DESCRIPTIONS_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(desc_rows[0].keys()))
        writer.writeheader()
        writer.writerows(desc_rows)

    print(f"\nWrote {FETCH_MANIFEST_CSV.name} and {DESCRIPTIONS_CSV.name}")
    if failures:
        print(f"WARNING: {len(failures)} fetch(es) fell back to the June cache: "
              f"{', '.join(failures)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
