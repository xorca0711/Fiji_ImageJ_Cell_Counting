#!/usr/bin/env python3
"""
aggregate_to_mouse.py
=====================================================================
Roll the per-region `run_summary.csv` produced by IF_Quant_Pipeline.groovy
up to the BIOLOGICAL replicate (mouse) level, then summarise per group.

Why this exists
---------------
The image pipeline emits one row per tissue/lesion region. Several sections
(and several regions per section) come from the SAME animal. Statistics for
this study must use n = MICE, not n = sections/regions. Averaging region rows
directly would pseudo-replicate and inflate significance.

Pooling is area-weighted (the statistically correct pooling for a section):
  * KRT5 pod area fraction   = sum(pod_area) / sum(tissue_area)      per mouse
  * marker density (/mm2)    = sum(morphology-authoritative pos_count)
                               / sum(tissue_area_mm2) per mouse
  * mean pod size            = sum(pod_area) / sum(n_pods)           per mouse
Counts and areas are summed; thresholds are averaged (QC only).

Outputs (written next to the input, or to --outdir):
  mouse_level_summary.csv  -- one row per (mouse_id, genotype, condition, panel)
  group_level_summary.csv  -- mean / sd / sem / n_mice per (genotype, condition,
                              panel, metric); n_mice is the real n for stats.

Usage
-----
  python3 aggregate_to_mouse.py /path/to/analysis_output/run_summary.csv
  python3 aggregate_to_mouse.py run_summary.csv --outdir ./stats

No third-party dependencies (standard library only).
=====================================================================
"""
import argparse
import csv
import math
import os
import sys
from collections import defaultdict

# Columns that identify a row rather than measure something.
KEY_COLS = ["mouse_id", "genotype", "condition", "panel"]
ROW_ID_COLS = ["image", "region", "section_id"]


def _num(v):
    """Parse a cell to float; blanks / non-numeric -> None."""
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.upper() == "NA":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def read_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            sys.exit(f"ERROR: {path} is empty or has no header.")
        rows = [r for r in reader if any((v or "").strip() for v in r.values())]
    return reader.fieldnames, rows


def validate_rows(header, rows):
    """Fail before aggregation when biological identity or row identity is unsafe."""
    missing_columns = [c for c in KEY_COLS + ROW_ID_COLS if c not in header]
    if missing_columns:
        sys.exit("ERROR: run_summary.csv is missing required columns: " +
                 ", ".join(missing_columns))

    invalid_mouse_rows = [
        r for r in rows
        if (r.get("mouse_id") or "").strip().upper() in {"", "NA", "N/A", "UNKNOWN"}
    ]
    if invalid_mouse_rows:
        examples = ", ".join((r.get("image") or "<unknown>") for r in invalid_mouse_rows[:5])
        sys.exit(
            f"ERROR: {len(invalid_mouse_rows)} row(s) lack a valid mouse_id ({examples}). "
            "Do not aggregate unknown animals into one pseudo-mouse; fix samplesheet.csv first."
        )

    identities = defaultdict(set)
    for r in rows:
        identities[(r.get("mouse_id") or "").strip()].add(
            ((r.get("genotype") or "NA").strip(), (r.get("condition") or "NA").strip())
        )
    conflicts = {mouse: values for mouse, values in identities.items() if len(values) > 1}
    if conflicts:
        details = "; ".join(f"{mouse}: {sorted(values)}" for mouse, values in list(conflicts.items())[:5])
        sys.exit("ERROR: a mouse_id maps to multiple genotype/condition identities: " + details)

    identity_column = "output_key" if "output_key" in header else "image"
    seen = set()
    duplicates = []
    for r in rows:
        key_columns = [identity_column, "region", "section_id", "panel"]
        key = tuple((r.get(c) or "").strip() for c in key_columns)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
    if duplicates:
        sys.exit(
            f"ERROR: {len(duplicates)} duplicate output/image-region-section-panel row(s) detected. "
            "Combine only one run_summary row per analyzed region; rerun old ambiguous summaries "
            "with a pipeline version that exports output_key."
        )


def classify_columns(header):
    """Group measurement columns by how they must be pooled."""
    # The plain <marker>_pos_count field is morphology-authoritative. Keep the
    # explicitly named raw-mean and state-audit fields in separate categories so
    # they cannot silently become the statistical endpoint.
    pos_count = [
        c for c in header
        if c.endswith("_pos_count")
        and not c.endswith("_morphology_pos_count")
        and not c.endswith("_raw_mean_pos_count")
        and not c.endswith("_true_pos_count")
    ]
    raw_mean_pos_count = [c for c in header if c.endswith("_raw_mean_pos_count")]
    morphology_pos_count = [c for c in header if c.endswith("_morphology_pos_count")]
    morphology_negative_count = [c for c in header if c.endswith("_morphology_negative_count")]
    marker_indeterminate_count = [
        c for c in header if c.endswith("_indeterminate_count")
        and not c.startswith("class_")
    ]
    morphology_evaluable_count = [c for c in header if c.endswith("_morphology_evaluable_count")]
    marker_audit_count = [
        c for c in header
        if c.endswith((
            "_raw_positive_final_negative_count",
            "_raw_negative_final_positive_count",
            "_intensity_morphology_discordant_count",
            "_review_burden_proxy_count",
        ))
    ]
    nucleus_qc_count = [
        c for c in header
        if c in {
            "n_rejected_nucleus_candidates", "n_rejected_below_min_area",
            "n_rejected_at_image_edge", "n_rejected_by_particle_filter",
            "n_nucleus_candidates_total",
        }
    ]
    positive_area = [c for c in header if c.endswith("_positive_area_um2")]
    n_components = [c for c in header if c.endswith("_n_components")]
    # total pod area only -- exclude the derived per-region MEAN pod size, which
    # also ends in "_pod_area_um2" and would otherwise be summed as an area.
    pod_area = [c for c in header
                if c.endswith("_pod_area_um2") and not c.endswith("_mean_pod_area_um2")]
    n_pods = [c for c in header if c.endswith("_n_pods")]
    class_count = [
        c for c in header if c.startswith("class_") and c.endswith("_count")
        and not c.endswith("_evaluable_count")
        and not c.endswith("_indeterminate_count")
    ]
    class_evaluable_count = [c for c in header if c.startswith("class_") and c.endswith("_evaluable_count")]
    class_indeterminate_count = [c for c in header if c.startswith("class_") and c.endswith("_indeterminate_count")]
    # everything else numeric-ish that we simply sum or average
    state_counts = (set(raw_mean_pos_count) | set(morphology_pos_count) |
                    set(morphology_negative_count) | set(marker_indeterminate_count) |
                    set(morphology_evaluable_count) | set(class_evaluable_count) |
                    set(class_indeterminate_count) | set(marker_audit_count))
    sum_cols = (set(["region_area_um2", "n_nuclei"]) | set(pos_count) |
                set(pod_area) | set(n_pods) | set(class_count) | state_counts |
                set(nucleus_qc_count) | set(positive_area) | set(n_components))
    # derived columns we recompute (do NOT sum): fractions, densities, mean pod size, thresholds
    return {
        "pos_count": pos_count,
        "raw_mean_pos_count": raw_mean_pos_count,
        "state_counts": sorted(state_counts),
        "marker_audit_count": marker_audit_count,
        "nucleus_qc_count": nucleus_qc_count,
        "positive_area": positive_area,
        "n_components": n_components,
        "pod_area": pod_area,
        "n_pods": n_pods,
        "class_count": class_count,
        "sum_cols": sum_cols,
    }


def marker_of(col, suffix):
    return col[: -len(suffix)]


def aggregate_mice(header, rows):
    cats = classify_columns(header)
    groups = defaultdict(list)
    for r in rows:
        key = tuple(r.get(k, "NA") for k in KEY_COLS)
        groups[key].append(r)

    out_rows = []
    for key, grp in sorted(groups.items()):
        mouse_id, genotype, condition, panel = key
        rec = {"mouse_id": mouse_id, "genotype": genotype,
               "condition": condition, "panel": panel}
        rec["n_regions"] = len(grp)
        rec["n_sections"] = len({(g.get("section_id") or "NA") for g in grp})

        # --- sums ---
        sums = {}
        for c in cats["sum_cols"]:
            vals = [_num(g.get(c)) for g in grp]
            vals = [v for v in vals if v is not None]
            sums[c] = sum(vals) if vals else 0.0

        total_area_um2 = sums.get("region_area_um2", 0.0)
        total_area_mm2 = total_area_um2 / 1e6
        rec["total_tissue_area_um2"] = total_area_um2
        rec["total_nuclei"] = sums.get("n_nuclei", 0.0)

        # --- nucleus-candidate QC totals and pooled fractions ---
        for c in cats["nucleus_qc_count"]:
            rec[f"{c}_total"] = sums[c]
        candidate_total = sums.get(
            "n_nucleus_candidates_total",
            sums.get("n_nuclei", 0.0) + sums.get("n_rejected_nucleus_candidates", 0.0),
        )
        rejected_total = sums.get("n_rejected_nucleus_candidates", 0.0)
        rec["nucleus_candidate_acceptance_fraction"] = (
            sums.get("n_nuclei", 0.0) / candidate_total if candidate_total > 0 else 0.0
        )
        rec["nucleus_candidate_rejection_fraction"] = (
            rejected_total / candidate_total if candidate_total > 0 else 0.0
        )
        for source, target in (
            ("n_rejected_below_min_area", "rejected_below_min_fraction_of_rejected"),
            ("n_rejected_at_image_edge", "rejected_edge_fraction_of_rejected"),
            ("n_rejected_by_particle_filter", "rejected_particle_filter_fraction_of_rejected"),
        ):
            rec[target] = sums.get(source, 0.0) / rejected_total if rejected_total > 0 else 0.0

        # --- marker positive counts + pooled density ---
        for c in cats["pos_count"]:
            m = marker_of(c, "_pos_count")
            rec[f"{m}_pos_count_total"] = sums[c]
            rec[f"{m}_density_per_mm2"] = (sums[c] / total_area_mm2) if total_area_mm2 > 0 else 0.0

        # --- explicit audit/state totals; never substitute for the endpoint ---
        for c in cats["raw_mean_pos_count"]:
            m = marker_of(c, "_raw_mean_pos_count")
            rec[f"{m}_raw_mean_pos_count_total"] = sums[c]
            rec[f"{m}_raw_mean_density_per_mm2"] = (
                sums[c] / total_area_mm2 if total_area_mm2 > 0 else 0.0
            )
        for c in cats["state_counts"]:
            if c in cats["raw_mean_pos_count"]:
                continue
            rec[f"{c}_total"] = sums[c]

        # Recompute morphology/QC fractions from pooled counts. Region-level
        # percentages must never be averaged because region sizes differ.
        for c in [x for x in header if x.endswith("_morphology_evaluable_count")]:
            marker = marker_of(c, "_morphology_evaluable_count")
            evaluable = sums.get(c, 0.0)
            included = sums.get("n_nuclei", 0.0)
            positive = sums.get(f"{marker}_morphology_pos_count", 0.0)
            negative = sums.get(f"{marker}_morphology_negative_count", 0.0)
            indeterminate = sums.get(f"{marker}_indeterminate_count", 0.0)
            discordant = sums.get(f"{marker}_intensity_morphology_discordant_count", 0.0)
            review = sums.get(f"{marker}_review_burden_proxy_count", indeterminate + discordant)
            rec[f"{marker}_morphology_positive_fraction_of_evaluable"] = positive / evaluable if evaluable > 0 else 0.0
            rec[f"{marker}_morphology_negative_fraction_of_evaluable"] = negative / evaluable if evaluable > 0 else 0.0
            rec[f"{marker}_indeterminate_fraction_of_included"] = indeterminate / included if included > 0 else 0.0
            rec[f"{marker}_intensity_morphology_discordant_fraction_of_evaluable"] = discordant / evaluable if evaluable > 0 else 0.0
            rec[f"{marker}_review_burden_proxy_fraction_of_included"] = review / included if included > 0 else 0.0

        # --- generic regional area endpoints (AcTub, membranes, reporter, ECM) ---
        for c in cats["positive_area"]:
            marker = marker_of(c, "_positive_area_um2")
            area = sums[c]
            components = sums.get(f"{marker}_n_components", 0.0)
            rec[f"{marker}_positive_area_um2_total"] = area
            rec[f"{marker}_positive_area_fraction"] = area / total_area_um2 if total_area_um2 > 0 else 0.0
            rec[f"{marker}_n_components_total"] = components
            rec[f"{marker}_mean_component_area_um2"] = area / components if components > 0 else 0.0

        # --- pod area, fraction, count, mean size (per area-marker, e.g. KRT5) ---
        for c in cats["pod_area"]:
            m = marker_of(c, "_pod_area_um2")
            pod_area = sums[c]
            npods = sums.get(f"{m}_n_pods", 0.0)
            rec[f"{m}_pod_area_um2_total"] = pod_area
            rec[f"{m}_pod_area_frac"] = (pod_area / total_area_um2) if total_area_um2 > 0 else 0.0
            rec[f"{m}_n_pods_total"] = npods
            rec[f"{m}_mean_pod_area_um2"] = (pod_area / npods) if npods > 0 else 0.0

        # --- classification counts + pooled density ---
        for c in cats["class_count"]:
            base = c[: -len("_count")]  # e.g. class_KRT5+_AGER-
            rec[f"{base}_count_total"] = sums[c]
            rec[f"{base}_density_per_mm2"] = (sums[c] / total_area_mm2) if total_area_mm2 > 0 else 0.0

        out_rows.append(rec)
    return out_rows


def _stats(values):
    vals = [v for v in values if v is not None]
    n = len(vals)
    if n == 0:
        return 0, 0.0, 0.0, 0.0
    mean = sum(vals) / n
    if n > 1:
        var = sum((v - mean) ** 2 for v in vals) / (n - 1)
        sd = math.sqrt(var)
        sem = sd / math.sqrt(n)
    else:
        sd = 0.0
        sem = 0.0
    return n, mean, sd, sem


def group_stats(mouse_rows):
    """mean/sd/sem/n across mice, per (genotype, condition, panel, metric)."""
    metric_cols = []
    seen = set()
    skip = set(KEY_COLS) | {"n_regions", "n_sections"}
    for r in mouse_rows:
        for c in r:
            if c in skip or c in seen:
                continue
            if isinstance(r[c], (int, float)):
                metric_cols.append(c)
                seen.add(c)

    groups = defaultdict(list)
    for r in mouse_rows:
        groups[(r["genotype"], r["condition"], r["panel"])].append(r)

    out = []
    for (geno, cond, panel), grp in sorted(groups.items()):
        for metric in metric_cols:
            vals = [r[metric] for r in grp if isinstance(r.get(metric), (int, float))]
            if not vals:
                continue
            n, mean, sd, sem = _stats(vals)
            out.append({"genotype": geno, "condition": cond, "panel": panel,
                        "metric": metric, "n_mice": n,
                        "mean": mean, "sd": sd, "sem": sem})
    return out


def write_csv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    cols = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(description="Aggregate per-region run_summary.csv to mouse and group level.")
    ap.add_argument("run_summary", help="path to run_summary.csv from the Fiji pipeline")
    ap.add_argument("--outdir", default=None, help="output folder (default: alongside input)")
    args = ap.parse_args()

    if not os.path.isfile(args.run_summary):
        sys.exit(f"ERROR: not found: {args.run_summary}")
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.run_summary))
    os.makedirs(outdir, exist_ok=True)

    header, rows = read_rows(args.run_summary)
    if not rows:
        sys.exit("ERROR: no data rows in run_summary.csv")
    validate_rows(header, rows)

    mouse_rows = aggregate_mice(header, rows)
    grp_rows = group_stats(mouse_rows)

    mouse_path = os.path.join(outdir, "mouse_level_summary.csv")
    group_path = os.path.join(outdir, "group_level_summary.csv")
    write_csv(mouse_path, mouse_rows)
    write_csv(group_path, grp_rows)

    n_mice = len({(r["mouse_id"], r["genotype"], r["condition"]) for r in mouse_rows})
    print(f"Read {len(rows)} region rows.")
    print(f"Wrote {len(mouse_rows)} mouse x panel rows -> {mouse_path}")
    print(f"Wrote {len(grp_rows)} group x metric rows -> {group_path}")
    print(f"Distinct animals: {n_mice}  (this is your statistical n, split by group)")
    print("Reminder: compare groups on the mouse-level metrics; n = mice.")


if __name__ == "__main__":
    main()
