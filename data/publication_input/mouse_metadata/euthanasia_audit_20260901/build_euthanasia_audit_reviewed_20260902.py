"""Build the reviewed post-landing timing audit (one row per accession).

For each of the 48 analysis accessions, this records:
  * the category of the description's euthanasia statement,
  * the post-landing euthanasia duration (where determinate),
  * an independent tissue-collection timing summary derived from the fresh API
    payload's ``samples.table`` dissection/euthanasia-date fields or explicit
    description text,
  * exact evidence sentence component(s) taken verbatim from the API
    ``description`` field, and
  * a short human note.

Rigor control: the fresh JSON payload is the primary source.  The API
``description`` is used for sentence-level evidence, while API
``samples.table`` is used for collection/dissection dates when available, with
a clearly labelled euthanasia-date proxy when no dissection field exists.  Each
evidence sentence component is verified to be a verbatim substring of the API
description.  The joined ``evidence_sentence`` columns are display fields;
their numbered components preserve exact sentence boundaries for auditability.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import date, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
DESCRIPTIONS_CSV = HERE / "study_descriptions_20260901.csv"
FETCH_MANIFEST_CSV = HERE / "fetch_manifest_20260901.csv"
FRESH_DIR = HERE / "repo_study_fresh_20260901"
# Keep the first-pass file untouched.  The reviewed file is the publication
# candidate produced by this script after the independent semantic audit.
OUT_CSV = HERE / "post_landing_tissue_collection_audit_20260902.csv"

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
        anchors=[],
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


# ---------------------------------------------------------------------------
# Independent tissue-collection timing audit.
#
# This is deliberately separate from J: the publication-facing question is
# when tissue was obtained, while J preserves the stricter euthanasia-only
# judgment.  API samples.table dissection dates are preferred; description
# text is used only when it explicitly states a tissue/dissection time.
_PLACEHOLDER_VALUES = {"", "not applicable", "not available", "n/a", "na"}
_SIMPLE_DATE_FORMATS = (
    "%d-%b-%Y", "%d-%b-%y", "%m/%d/%Y", "%m/%d/%y",
    "%m-%b-%Y", "%b-%d-%Y", "%b-%d-%y",
)


def _parse_simple_date(value: str) -> date | None:
    value = value.strip()
    for fmt in _SIMPLE_DATE_FORMATS:
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            continue
    return None


def _single_mission_end(payload: dict) -> date | None:
    raw = str(payload.get("missionEnd") or "").strip()
    if not raw or "," in raw:
        return None
    return _parse_simple_date(raw)


def _fmt_day_range(days: list[int]) -> str:
    vals = sorted(set(days))
    if not vals:
        return ""
    return str(vals[0]) if len(vals) == 1 else f"{vals[0]}-{vals[-1]}"


def _is_flight_group(group: str) -> bool:
    """Return whether a factor label represents the flight cohort.

    Ground, basal, vivarium, and cohort-control labels do not have a landing
    event.  Their raw dates stay in ``source_value`` but must not be turned into
    a post-landing interval by subtracting the flight mission end date.
    """
    normalized = re.sub(r"[^a-z]+", " ", group.lower()).strip()
    return normalized in {"flight", "space flight", "flt"}


def _api_tissue_record(payload: dict) -> dict:
    """Summarize explicit dissection-date fields in one API payload.

    The returned day strings are group-qualified when the API provides a
    ``Factor Value[Spaceflight]``.  They are differences from the API
    ``missionEnd`` date, not inferred euthanasia durations.
    """
    rows = payload.get("samples", {}).get("table", [])
    keys = sorted({
        key for row in rows for key in row
        if "dissection" in key.lower()
        and "date" in key.lower()
        and "condition" not in key.lower()
    })
    date_kind = "dissection"
    if not keys:
        # Some API records expose the terminal tissue procedure only as an
        # euthanasia date.  Keep this as an explicitly labelled proxy rather
        # than silently presenting it as a dissection timestamp.
        keys = sorted({
            key for row in rows for key in row
            if "euthanasia" in key.lower()
            and "date" in key.lower()
            and "condition" not in key.lower()
        })
        date_kind = "euthanasia"
    by_field_group: dict[str, dict[str, Counter[str]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    for row in rows:
        group = str(row.get("Factor Value[Spaceflight]", "")).strip() or "all API samples"
        for key in keys:
            value = str(row.get(key, "")).strip()
            if value.lower() in _PLACEHOLDER_VALUES:
                continue
            by_field_group[key][group][value] += 1

    if not by_field_group:
        return dict(
            days="", category="not_stated", scope="n/a",
            confidence="not_stated", source_field="", source_value="",
            anchors=[], notes=(
                "Fresh API samples.table has no explicit dissection-date field, "
                "and no separate tissue-collection timing was stated."
            ),
        )

    source_parts: list[str] = []
    scope_groups: set[str] = set()
    for key in keys:
        for group in sorted(by_field_group[key]):
            scope_groups.add(group)
            counts = by_field_group[key][group]
            values = "; ".join(
                f"{value} (n={counts[value]})" for value in sorted(counts)
            )
            source_parts.append(f"{key} [{group}]: {values}")

    mission_end = _single_mission_end(payload)
    days_by_group: dict[str, list[int]] = defaultdict(list)
    flight_groups = {
        group for key in keys for group in by_field_group[key]
        if _is_flight_group(group)
    }
    # A landing interval is defined only for a flight group.  If the API field
    # contains control-only records, retain those raw values but do not invent
    # a landing anchor for them.
    unresolved = mission_end is None or not flight_groups
    for key in keys:
        for group, counts in by_field_group[key].items():
            if not _is_flight_group(group):
                continue
            for value, count in counts.items():
                parsed = _parse_simple_date(value)
                if parsed is None or mission_end is None:
                    unresolved = True
                    continue
                delta = (parsed - mission_end).days
                if delta < 0:
                    unresolved = True
                    continue
                days_by_group[group].extend([delta] * count)

    days_parts = [
        f"{group}={_fmt_day_range(days_by_group[group])}"
        for group in sorted(days_by_group)
        if days_by_group[group]
    ]
    days = "; ".join(days_parts)
    source_field = "API samples.table: " + " | ".join(keys)
    scope = "; ".join(sorted(scope_groups))
    if days and not unresolved:
        if date_kind == "dissection":
            category = "post_landing_api_dissection_date_derived"
            confidence = "api_sample_dissection_date"
            notes = (
                f"API samples.table provides dissection date(s) for the flight "
                f"group; day differences use the API missionEnd "
                f"{payload.get('missionEnd')!r}. Control-group dates are retained "
                "as raw provenance but have no landing interval. This is a tissue "
                "collection/dissection timing and is not a separately inferred "
                "euthanasia interval."
            )
        else:
            category = "post_landing_api_euthanasia_date_proxy"
            confidence = "api_sample_euthanasia_date_proxy"
            notes = (
                f"API samples.table provides euthanasia date(s) for the flight "
                f"group; day differences use the API missionEnd "
                f"{payload.get('missionEnd')!r}. Report this as a terminal "
                "tissue-acquisition proxy because no dissection date field is "
                "available; control-group dates have no landing interval."
            )
    else:
        category = "api_dissection_date_present_interval_unresolved"
        confidence = "api_sample_dissection_date_unresolved"
        notes = (
            "API samples.table provides dissection date text, but the available "
            "missionEnd/date formats or mixed cohorts do not support one "
            "unambiguous post-landing day value; raw API values are retained."
        )
        # Do not expose a partial numeric summary when any group/date is not
        # safely comparable to one mission-end date.
        days = ""
    return dict(
        days=days, category=category, scope=scope,
        confidence=confidence, source_field=source_field,
        source_value=" | ".join(source_parts), anchors=[], notes=notes,
    )


def _tissue_manual(
    *, days: str, category: str, scope: str, confidence: str,
    source_value: str, anchors: list[str], notes: str,
) -> dict:
    return dict(
        days=days, category=category, scope=scope, confidence=confidence,
        source_field="API description", source_value=source_value,
        anchors=anchors, notes=notes,
    )


def _rr23_tissue(days: str = "FLT=1") -> dict:
    return _tissue_manual(
        days=days,
        category="post_landing_description_dissection_date_derived",
        scope="Flight (post-landing); HGC/VGC controls (no landing interval)",
        confidence="api_description_dissection_date",
        source_value="API description gives return on Jan 13 and dissection on Jan 14/17/20",
        anchors=["euthanized and dissected on Jan 14th"],
        notes=(
            "The API description explicitly says the animals were dissected on "
            "Jan 14, 17 or 20, respectively. The tissue field reports the "
            "flight dissection timing (FLT=1 d after the Jan 13 return); HGC and "
            "VGC are controls without a landing interval. It does not make an "
            "additional euthanasia assumption."
        ),
    )


TISSUE_MANUAL = {
    "OSD-47": _tissue_manual(
        days="", category="on_orbit", scope="flight (dissected on ISS)",
        confidence="api_description_on_orbit",
        source_value="API description: dissected on-orbit (21 or 22 days after launch)",
        anchors=["dissected on-orbit (21 or 22 days after launch)"],
        notes="Tissue was obtained on-orbit; a post-landing interval is not applicable.",
    ),
    "OSD-239": _tissue_manual(
        days="<1", category="post_landing_explicit", scope="flight",
        confidence="api_description_dissection_date",
        source_value="API description explicitly states dissection less than 1 day after splashdown",
        anchors=[
            "Mice were returned live and euthanized and dissected less than 1 day after splashdown.",
            "Femoral skin was dissected 30 minutes after euthanasia and snap frozen in liquid nitrogen.",
        ],
        notes=(
            "Uses the API description's explicit dissection timing; the value is "
            "reported as tissue acquisition timing, not inferred from the euthanasia date."
        ),
    ),
    "OSD-240": _tissue_manual(
        days="30", category="post_landing_euthanasia_proxy", scope="flight LAR",
        confidence="api_sample_euthanasia_date_proxy",
        source_value=(
            "API description: animals were allowed to recover for 30 days in "
            "standard habitats before euthanasia"
        ),
        anchors=["allowed to recover for 30 days in standard habitats before euthanasia"],
        notes=(
            "The API description gives a 30-day live-animal-return recovery before "
            "euthanasia. Because tissue was obtained from these animals at the "
            "terminal procedure, report 30 d as a tissue-acquisition proxy; the "
            "description does not give a separate dissection timestamp."
        ),
    ),
    "OSD-241": _tissue_manual(
        days="30", category="post_landing_euthanasia_proxy", scope="flight LAR",
        confidence="api_sample_euthanasia_date_proxy",
        source_value=(
            "API description: animals were allowed to recover for 30 days in "
            "standard habitats before euthanasia"
        ),
        anchors=["allowed to recover for 30 days in standard habitats before euthanasia"],
        notes=(
            "The API description gives a 30-day live-animal-return recovery before "
            "euthanasia. Because tissue was obtained from these animals at the "
            "terminal procedure, report 30 d as a tissue-acquisition proxy; the "
            "description does not give a separate dissection timestamp."
        ),
    ),
    "OSD-242": _tissue_manual(
        days="1-2", category="post_landing_date_conditional", scope="flight",
        confidence="conditional_api_description_dissection_date",
        source_value="API description: euthanasia and dissection on 9/18/2018; API missionEnd=9/17/17",
        anchors=["testing, euthanasia and dissection on 9/18/2018"],
        notes=(
            "The API description gives the flight dissection date, but its year "
            "conflicts with the API mission record. Keep the conservative 1-2 d "
            "conditional range rather than a fixed integer."
        ),
    ),
    "OSD-379": _tissue_manual(
        days="LAR=2", category="mixed_post_landing_and_on_orbit",
        scope="LAR flight (analyzed); ISS-T flight on-orbit",
        confidence="api_description_euthanasia_proxy",
        source_value=(
            "API description: LAR animals were returned live, allowed to recover "
            "for 2 days, then sacrificed; ISS-T animals were sacrificed on-orbit"
        ),
        anchors=["allowed to recover for 2 days (Live Animal Return, LAR) before sacrifice"],
        notes=(
            "The analyzed LAR flight samples have a 2-day post-landing recovery "
            "before the terminal tissue procedure. The ISS-T subgroup was handled "
            "on-orbit, so it has no post-landing interval; the API sample table also "
            "contains mixed basal-control dissection text that is retained below "
            "as raw API provenance."
        ),
    ),
    "OSD-665": _rr23_tissue(),
    "OSD-666": _rr23_tissue(),
    "OSD-576": _rr23_tissue(),
    "OSD-580": _tissue_manual(
        days="LAR=24", category="mixed_post_landing_and_on_orbit",
        scope="LAR flight (analyzed); ISS-T flight on-orbit",
        confidence="api_description_euthanasia_proxy",
        source_value=(
            "API description: LAR animals were returned live, allowed to recover "
            "for 24 days, then sacrificed; ISS-T animals were sacrificed on-orbit"
        ),
        anchors=["allowed to recover for 24 days (Live Animal Return, LAR) before sacrifice"],
        notes=(
            "The analyzed flight samples are LAR samples and therefore have a "
            "24-day post-landing recovery before the terminal tissue procedure. "
            "The ISS-T subgroup was handled on-orbit and has no post-landing interval."
        ),
    ),
    "OSD-771": _tissue_manual(
        days="LAR=24", category="mixed_post_landing_and_on_orbit",
        scope="LAR flight (analyzed); ISS-T flight on-orbit",
        confidence="api_description_euthanasia_proxy",
        source_value=(
            "API description: LAR animals were returned live, allowed to recover "
            "for 24 days, then sacrificed; ISS-T animals were sacrificed on-orbit"
        ),
        anchors=["allowed to recover for 24 days (Live Animal Return, LAR) before sacrifice"],
        notes=(
            "The analyzed flight samples are LAR samples and therefore have a "
            "24-day post-landing recovery before the terminal tissue procedure. "
            "The ISS-T subgroup was handled on-orbit and has no post-landing interval."
        ),
    ),
    "OSD-714": _tissue_manual(
        days="2", category="post_landing_explicit", scope="flight (skeletal muscle)",
        confidence="api_description_explicit",
        source_value="API description: skeletal muscle was collected two days after landing",
        anchors=["All mice survived and returned to Earth, and skeletal muscle was collected two days after landing."],
        notes=(
            "This is the tissue-collection timing requested for the publication "
            "summary. The API description does not separately state euthanasia timing."
        ),
    ),
    "OSD-770": _rr23_tissue(),
}


TISSUE_UNRESOLVED_API_DATES = {
    "OSD-253": "API dates apply to Ground Control Rerun records and do not uniquely describe the analyzed flight cohort.",
    "OSD-254": "API dissection dates span multiple years/cohorts; no single post-landing interval is safe to collapse.",
    "OSD-379": "API dissection-date field contains mixed textual ranges (ISS-T/LAR cohorts); no single day interval is safe to collapse.",
}


def _tissue_record(acc: str, payload: dict) -> dict:
    record = _api_tissue_record(payload)
    if acc in TISSUE_MANUAL:
        manual = dict(TISSUE_MANUAL[acc])
        # Retain raw API samples.table values when the description supplies a
        # more useful cohort-qualified timing (for example RRRM-1 OSD-379).
        if record.get("source_field") and record.get("source_value"):
            manual["source_field"] = (
                f"{record['source_field']} | {manual['source_field']}"
            )
            manual["source_value"] = (
                f"{record['source_value']} | {manual['source_value']}"
            )
        record.update(manual)
    if acc in TISSUE_UNRESOLVED_API_DATES and acc not in TISSUE_MANUAL:
        record["days"] = ""
        record["category"] = "api_dissection_date_present_interval_unresolved"
        record["confidence"] = "api_sample_dissection_date_unresolved"
        record["notes"] = TISSUE_UNRESOLVED_API_DATES[acc]
    if acc == "OSD-419":
        record["category"] = "mixed_on_orbit_and_post_landing"
        record["notes"] = (
            "API samples.table provides final dissection dates for the analyzed "
            "Flight/Ground records; the description also notes that two flight "
            "animals were dissected on the ISS. Read the scope and raw API values "
            "rather than treating every sample as one collection event."
        )
    return record


def _sentences(text: str) -> list[str]:
    parts = re.split(r"(?<=[.!?])\s+", text.strip())
    return [p for p in parts if p]


def _load_api_payloads() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for path in sorted(FRESH_DIR.glob("OSD-*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        acc = str(payload.get("accession") or path.stem).strip()
        if acc in payloads:
            raise SystemExit(f"Duplicate API payload accession: {acc}")
        payloads[acc] = payload
    if not payloads:
        raise SystemExit(f"No saved API payloads found in {FRESH_DIR}")
    if FETCH_MANIFEST_CSV.exists():
        manifest = {
            r["accession"]: r
            for r in csv.DictReader(
                FETCH_MANIFEST_CSV.open(newline="", encoding="utf-8")
            )
        }
        if set(manifest) != set(payloads):
            raise SystemExit("Saved API payloads and fetch manifest do not match")
        for acc, payload in payloads.items():
            if manifest[acc].get("fetch_status") != "ok":
                raise SystemExit(
                    f"API payload for {acc} was not fetched fresh: "
                    f"{manifest[acc].get('fetch_status')!r}"
                )
            expected = manifest[acc].get("fresh_sha256", "")
            actual = hashlib.sha256(
                json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
            ).hexdigest()
            if expected and expected != actual:
                raise SystemExit(f"API payload hash mismatch for {acc}")
    # The CSV is a useful cross-check, but never the primary source for the
    # reviewed judgments.
    if DESCRIPTIONS_CSV.exists():
        csv_desc = {
            r["accession"]: r["description"]
            for r in csv.DictReader(DESCRIPTIONS_CSV.open(newline="", encoding="utf-8"))
        }
        for acc, payload in payloads.items():
            if acc in csv_desc and payload.get("description", "") != csv_desc[acc]:
                raise SystemExit(f"API/description CSV mismatch for {acc}")
    return payloads


def main() -> int:
    payloads = _load_api_payloads()
    rows_desc = {
        acc: {
            "accession": acc,
            "title": payload.get("title", ""),
            "mission_name": payload.get("missionName", ""),
            "description": payload.get("description", ""),
        }
        for acc, payload in payloads.items()
    }
    if set(J) != set(rows_desc):
        missing = set(J) - set(rows_desc)
        extra = set(rows_desc) - set(J)
        raise SystemExit(f"Accession mismatch: missing={missing} extra={extra}")

    out_rows: list[dict] = []
    problems: list[str] = []

    for acc in sorted(rows_desc):
        spec = J[acc]
        tissue = _tissue_record(acc, payloads[acc])
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

        tissue_matched: list[str] = []
        for anchor in tissue["anchors"]:
            hits = [s for s in sents if anchor in s]
            if not hits:
                problems.append(f"{acc}: tissue anchor not found -> {anchor!r}")
            for hit in hits:
                if hit not in tissue_matched:
                    tissue_matched.append(hit)
        for sentence in tissue_matched:
            if sentence not in desc:
                problems.append(
                    f"{acc}: tissue sentence not verbatim -> {sentence!r}"
                )
        tissue_evidence = " ... ".join(tissue_matched) if tissue_matched else ""

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
            "tissue_collection_after_landing_days": tissue["days"],
            "tissue_collection_timing_category": tissue["category"],
            "tissue_collection_scope": tissue["scope"],
            "tissue_collection_confidence": tissue["confidence"],
            "tissue_collection_source_field": tissue["source_field"],
            "tissue_collection_source_value": tissue["source_value"],
            "tissue_collection_evidence_sentence": tissue_evidence,
            "tissue_collection_evidence_sentence_1": (
                tissue_matched[0] if len(tissue_matched) >= 1 else ""
            ),
            "tissue_collection_evidence_sentence_2": (
                tissue_matched[1] if len(tissue_matched) >= 2 else ""
            ),
            "tissue_collection_evidence_sentence_count": str(len(tissue_matched)),
            "tissue_collection_notes": tissue["notes"],
        })

    fieldnames = ["accession", "title", "mission_name", "category", "confidence",
                  "euthanasia_after_landing_days", "euthanasia_group_scope",
                  "evidence_sentence", "evidence_sentence_1", "evidence_sentence_2",
                  "evidence_sentence_count", "source_field", "notes",
                  "tissue_collection_after_landing_days",
                  "tissue_collection_timing_category", "tissue_collection_scope",
                  "tissue_collection_confidence", "tissue_collection_source_field",
                  "tissue_collection_source_value",
                  "tissue_collection_evidence_sentence",
                  "tissue_collection_evidence_sentence_1",
                  "tissue_collection_evidence_sentence_2",
                  "tissue_collection_evidence_sentence_count",
                  "tissue_collection_notes"]
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
