"""Build the post-landing euthanasia-timing audit table (one row per accession).

For each of the 48 analysis accessions, this records:
  * the category of the description's euthanasia statement,
  * the post-landing euthanasia duration (where determinate),
    * exact evidence sentence component(s) taken verbatim from the freshly
      fetched ``description`` field, and
  * a short human note.

Rigor control: every evidence sentence component is verified to be a verbatim
substring of the corresponding accession's description.  The joined
``evidence_sentence`` column is a display field; ``evidence_sentence_1`` and
``evidence_sentence_2`` preserve the exact sentence boundaries for auditability.
Any anchor that fails to locate a sentence is reported, so a transcription slip
cannot silently enter the audit.
Only the ``description`` field is used (per the task); no other field is read.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DESCRIPTIONS_CSV = HERE / "study_descriptions_20260901.csv"
# Keep the first-pass file untouched.  The reviewed file is the publication
# candidate produced by this script after the independent semantic audit.
OUT_CSV = HERE / "euthanasia_after_landing_audit_20260901_reviewed.csv"

# ---------------------------------------------------------------------------
# Curated per-accession judgments.
#
#   anchors : substrings that locate the evidence sentence(s) in the description
#   category: post_landing | post_landing_date_derived | on_orbit | not_stated
#   days    : post-landing duration string ("" when not determinate; ranges are allowed)
#   scope   : which cohort the value applies to
#   conf    : explicit | date_derived | conditional_date_derived |
#             inferred_on_orbit | proxy_tissue_collection | not_stated
#   notes   : context for a human reviewer
# ---------------------------------------------------------------------------
def _rr23(scope: str = "flight (HGC/VGC are controls)") -> dict:
    return dict(
        anchors=["returned to Earth alive (Jan", "euthanized and dissected on Jan 14th"],
        category="post_landing_date_derived", days="1", scope=scope,
        conf="date_derived",
        notes="RR-23: flight mice returned alive Jan 13 2021; flight animals euthanized "
              "Jan 14 2021 => ~1 d after landing (HGC Jan 17, VGC Jan 20 are ground controls).")


def _notstated(note: str) -> dict:
    return dict(anchors=[], category="not_stated", days="", scope="n/a",
                conf="not_stated", notes=note)


J = {
    # --- Explicit post-landing day counts ---------------------------------
    "OSD-239": dict(
        anchors=["less than 1 day after splashdown"],
        category="post_landing", days="<1", scope="flight", conf="explicit",
        notes="Returned live; euthanized and dissected <1 day after splashdown."),
    "OSD-240": dict(
        anchors=["recover for 30 days in standard habitats before euthanasia"],
        category="post_landing", days="30", scope="flight (LAR)", conf="explicit",
        notes="Live animal return; 30-day recovery before euthanasia (ketamine/xylazine)."),
    "OSD-241": dict(
        anchors=["recover for 30 days in standard habitats before euthanasia"],
        category="post_landing", days="30", scope="flight (LAR)", conf="explicit",
        notes="Live animal return; 30-day recovery before euthanasia (ketamine/xylazine)."),
    "OSD-714": dict(
        anchors=["skeletal muscle was collected two days after landing"],
        category="not_stated", days="", scope="flight (tissue collected 2 d; euthanasia not stated)",
        conf="proxy_tissue_collection",
        notes="Description states that skeletal muscle was collected 2 days after landing, "
              "but does not state when the mice were euthanized; do not use 2 d as a "
              "euthanasia interval."),
    # --- RRRM live-animal-return: explicit recovery, on-orbit subgroup -----
    "OSD-379": dict(
        anchors=["allowed to recover for 2 days (Live Animal Return, LAR) before sacrifice"],
        category="post_landing", days="2", scope="flight_LAR only (ISS-T subgroup on-orbit)",
        conf="explicit",
        notes="RRRM-1: LAR animals recovered 2 days before sacrifice; the ISS-T subgroup "
              "was sacrificed on-orbit after 22-23 days."),
    "OSD-511": dict(
        anchors=["allowed to recover for 2 days (Live Animal Return, LAR) before sacrifice"],
        category="post_landing", days="2", scope="flight_LAR only (ISS-T subgroup on-orbit)",
        conf="explicit",
        notes="RRRM-1 (mammary): LAR recovered 2 days; ISS-T subgroup sacrificed on-orbit."),
    "OSD-580": dict(
        anchors=["allowed to recover for 24 days (Live Animal Return, LAR) before sacrifice"],
        category="post_landing", days="24", scope="flight_LAR only (ISS-T subgroup on-orbit)",
        conf="explicit",
        notes="RRRM-2 (heart): LAR recovered 24 days before sacrifice; ISS-T subgroup "
              "sacrificed on-orbit after 55-58 days."),
    "OSD-771": dict(
        anchors=["allowed to recover for 24 days (Live Animal Return, LAR) before sacrifice"],
        category="post_landing", days="24", scope="flight_LAR only (ISS-T subgroup on-orbit)",
        conf="explicit",
        notes="RRRM-2 (kidney): LAR recovered 24 days before sacrifice; ISS-T on-orbit."),
    # --- Date-derived post-landing timing ---------------------------------
    "OSD-242": dict(
        anchors=["for testing, euthanasia and dissection on 9/18/2018"],
        category="post_landing_date_derived", days="1-2", scope="flight",
        conf="conditional_date_derived",
        notes="Description gives euthanasia on 9/18/2018, while the same fresh API record "
              "identifies the SpaceX-12 mission as ending 9/17/2017; the year is a source "
              "typo and the description alone cannot determine an exact integer interval. "
              "Record as conditional ~1-2 d, not a confirmed 2 d."),
    "OSD-420": dict(
        anchors=["for testing, euthanasia, and dissection on 9/18/2018"],
        category="post_landing_date_derived", days="1-2", scope="flight",
        conf="conditional_date_derived",
        notes="Description gives euthanasia on 9/18/2018, while the same fresh API record "
              "identifies the SpaceX-12 mission as ending 9/17/2017; the year is a source "
              "typo and the description alone cannot determine an exact integer interval. "
              "Record as conditional ~1-2 d, not a confirmed 2 d."),
    "OSD-421": dict(
        anchors=["for testing, euthanasia, and dissection on 9/18/2018"],
        category="post_landing_date_derived", days="1-2", scope="flight",
        conf="conditional_date_derived",
        notes="Description gives euthanasia on 9/18/2018, while the same fresh API record "
              "identifies the SpaceX-12 mission as ending 9/17/2017; the year is a source "
              "typo and the description alone cannot determine an exact integer interval. "
              "Record as conditional ~1-2 d, not a confirmed 2 d."),
    # --- RR-23: returned alive Jan 13 2021, flight euthanized Jan 14 2021 --
    "OSD-463": _rr23(),
    "OSD-464": _rr23(),
    "OSD-506": _rr23(),
    "OSD-512": _rr23(),
    "OSD-513": _rr23(),
    "OSD-515": _rr23(),
    "OSD-525": _rr23(),
    "OSD-576": _rr23(),
    "OSD-665": _rr23(),
    "OSD-666": _rr23(),
    "OSD-770": _rr23(scope="flight (GC/VIV are controls)"),
    # --- On-orbit euthanasia (post-landing time not applicable) ----------
    "OSD-419": dict(
        anchors=["euthanized and preserved intact on-orbit for subsequent dissection on the ground"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-1: flight mice euthanized on-orbit after 37 d; carcasses stowed in MELFI "
              "and returned. Post-landing euthanasia time not applicable."),
    "OSD-47": dict(
        anchors=["dissected on-orbit (21 or 22 days after launch)"],
        category="on_orbit", days="", scope="all flight (dissected on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-1 CASIS: flight group dissected on-orbit (21-22 d after launch). "
              "Post-landing time not applicable."),
    "OSD-253": dict(
        anchors=["euthanized on the ISS with Ketamine/Xylazine/Acepromazine"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-7: half euthanized on the ISS after 25 d, the rest after 75-76 d; "
              "carcasses returned on SpX-15/16. Post-landing time not applicable."),
    "OSD-254": dict(
        anchors=["euthanized on the ISS with Ketamine/Xylazine/Acepromazine"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-7 (dorsal skin): same on-orbit schedule as OSD-253. "
              "Post-landing time not applicable."),
    "OSD-462": dict(
        anchors=["28-29 days in microgravity, the Flight mice were euthanized"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-10: flight mice euthanized on-orbit after 28-29 d; carcasses cryopreserved "
              "and returned. Post-landing time not applicable."),
    "OSD-667": dict(
        anchors=["28-29 days in microgravity, the Flight mice were euthanized"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-10 (colon): same on-orbit schedule as OSD-462. "
              "Post-landing time not applicable."),
    "OSD-899": dict(
        anchors=["28-29 days in microgravity, the Flight mice were euthanized"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-10 (cecum): same on-orbit schedule as OSD-462. "
              "Post-landing time not applicable."),
    "OSD-900": dict(
        anchors=["28-29 days in microgravity, the Flight mice were euthanized"],
        category="on_orbit", days="", scope="all flight (euthanized on the ISS)",
        conf="inferred_on_orbit",
        notes="RR-10 (lung): same on-orbit schedule as OSD-462. "
              "Post-landing time not applicable."),
    # --- Not stated in the description -----------------------------------
    "OSD-99": _notstated("RR-1 archived-tissue re-analysis; description says tissue was "
                          "'stored at least a year at -80C after return to Earth' but gives no "
                          "post-landing euthanasia timing."),
    "OSD-100": _notstated("RR-1 archived-tissue re-analysis; tissue 'stored at least a year at "
                           "-80C after return to Earth'; no post-landing euthanasia timing."),
    "OSD-101": _notstated("RR-1 archived-tissue re-analysis; tissue 'stored at least a year at "
                           "-80C after return to Earth'; no post-landing euthanasia timing."),
    "OSD-104": _notstated("RR-1 archived-tissue re-analysis; tissue 'stored at least a year at "
                           "-80C after return to Earth'; no post-landing euthanasia timing."),
    "OSD-137": _notstated("Description states animals 'were euthanized' but gives no number of "
                           "days after landing."),
    "OSD-161": _notstated("Description states animals 'were euthanized' but gives no number of "
                           "days after landing."),
    "OSD-162": _notstated("Description states animals 'were euthanized' but gives no number of "
                           "days after landing."),
    "OSD-163": _notstated("Description states animals 'were euthanized' but gives no number of "
                           "days after landing."),
    "OSD-194": _notstated("Description states animals 'were euthanized' but gives no number of "
                           "days after landing."),
    "OSD-255": _notstated("Description states the study purpose only; no euthanasia timing."),
    "OSD-289": _notstated("Description states the study purpose only; no euthanasia timing."),
    "OSD-397": _notstated("RR-1 cohort (mice 'kept inside ISS middeck for 37 days'); this "
                           "accession's description gives no post-landing euthanasia timing."),
    "OSD-48": _notstated("Secondary dataset-mining analysis; description gives no euthanasia "
                          "timing."),
    "OSD-599": _notstated("Description says mice were 'flown on the ISS for 30 days' but gives no "
                           "post-landing euthanasia timing."),
    "OSD-686": _notstated("Description says mice were 'on board the International Space Station "
                           "for 30 days' but gives no post-landing euthanasia timing."),
    "OSD-690": _notstated("Description says mice 'travelled in space for 31 days' but gives no "
                           "post-landing euthanasia timing."),
    "OSD-758": _notstated("Description says tissue processed 'after returning them to Earth "
                           "alive' but gives no number of days after landing."),
    "OSD-759": _notstated("Description says tissue processed 'after returning them to Earth "
                           "alive' but gives no number of days after landing."),
}


def _sentences(text: str) -> list[str]:
    parts = re.split(r"(?<=[.!?])\s+", text.strip())
    return [p for p in parts if p]


def main() -> int:
    rows_desc = {r["accession"]: r for r in
                 csv.DictReader(DESCRIPTIONS_CSV.open(newline="", encoding="utf-8"))}
    if set(J) != set(rows_desc):
        missing = set(J) - set(rows_desc)
        extra = set(rows_desc) - set(J)
        raise SystemExit(f"Accession mismatch: missing={missing} extra={extra}")

    out_rows: list[dict] = []
    problems: list[str] = []

    for acc in sorted(rows_desc):
        spec = J[acc]
        desc = rows_desc[acc]["description"]
        sents = _sentences(desc)
        matched: list[str] = []
        for anchor in spec["anchors"]:
            hits = [s for s in sents if anchor in s]
            if not hits:
                problems.append(f"{acc}: anchor not found -> {anchor!r}")
            for h in hits:
                if h not in matched:
                    matched.append(h)
        # verbatim check: every matched sentence must be a substring of desc
        for m in matched:
            if m not in desc:
                problems.append(f"{acc}: matched sentence not verbatim -> {m!r}")
        evidence = " ... ".join(matched) if matched else ""

        out_rows.append({
            "accession": acc,
            "title": rows_desc[acc]["title"],
            "mission_name": rows_desc[acc]["mission_name"],
            "category": spec["category"],
            "confidence": spec["conf"],
            "euthanasia_after_landing_days": spec["days"],
            "euthanasia_group_scope": spec["scope"],
            "evidence_sentence": evidence,
            "evidence_sentence_1": matched[0] if len(matched) >= 1 else "",
            "evidence_sentence_2": matched[1] if len(matched) >= 2 else "",
            "evidence_sentence_count": str(len(matched)),
            "source_field": "description",
            "notes": spec["notes"],
        })

    fieldnames = ["accession", "title", "mission_name", "category", "confidence",
                  "euthanasia_after_landing_days", "euthanasia_group_scope",
                  "evidence_sentence", "evidence_sentence_1", "evidence_sentence_2",
                  "evidence_sentence_count", "source_field", "notes"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {OUT_CSV.name} with {len(out_rows)} rows.")
    # category summary
    from collections import Counter
    c = Counter(r["category"] for r in out_rows)
    print("category counts:", dict(c))
    if problems:
        print("\n!!! VERIFICATION PROBLEMS:")
        for p in problems:
            print("   ", p)
        return 1
    print("Verification OK: each evidence sentence component is a verbatim substring.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
