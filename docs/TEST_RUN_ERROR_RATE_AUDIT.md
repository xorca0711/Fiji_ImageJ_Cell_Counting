# Test-run error-proxy audit

Audit date: 2026-07-22  
Primary data set: 45 canonical 20x fields  
Sensitivity data set: all 47 attempts, including two exactly repeated G002 pilot analyses

## What “error rate” means here

No manually annotated ground truth exists for these fields. Therefore, none of
the percentages below is a validated false-positive or false-negative rate.
They are deliberately separated into three QC proxies:

1. **Nucleus-candidate rejection:** DAPI-derived objects removed before cell
   counting.
2. **Indeterminate marker calls:** included nuclei for which a marker could not
   safely receive a positive or negative morphology call.
3. **Intensity–morphology disagreement:** evaluable cells for which the legacy
   mean-intensity decision differs from the morphology-authoritative decision.

The recommended longitudinal script metric is **review burden**, defined as
`indeterminate + intensity–morphology discordant`. These categories do not
overlap.

## Overall canonical percentages

### Nucleus segmentation

| Outcome | Count | Denominator | Percentage |
|---|---:|---:|---:|
| Accepted nuclei | 132,384 | 255,035 candidates | 51.91% |
| Rejected candidates | 122,651 | 255,035 candidates | 48.09% |
| Rejected below minimum area | 116,476 | 122,651 rejected | 94.97% |
| Rejected at image edge | 6,175 | 122,651 rejected | 5.03% |
| Rejected by particle filter | 0 | 122,651 rejected | 0.00% |

The pooled rejection rate is distorted by panel-R field G014, which has a
96.36% candidate-rejection rate. Excluding only that field gives 40.86%; the
median field-level rejection rate is 40.37%. This is why a single pooled
“segmentation error rate” must not be used to tune the minimum nucleus area.

### Assignable marker calls

Acetylated tubulin is excluded from this combined denominator because its
whole-field per-cell calls were intentionally disabled in the absence of an
airway compartment ROI.

| Outcome | Count | Denominator | Percentage |
|---|---:|---:|---:|
| Evaluable | 286,799 | 301,602 marker-cell opportunities | 95.09% |
| Indeterminate | 14,803 | 301,602 opportunities | 4.91% |
| Morphology-positive | 89,331 | 286,799 evaluable | 31.15% |
| Morphology-negative | 197,468 | 286,799 evaluable | 68.85% |
| Intensity–morphology concordant | 251,970 | 286,799 evaluable | 87.86% |
| Intensity–morphology discordant | 34,829 | 286,799 evaluable | 12.14% |
| Raw positive changed to final negative | 1,010 | 286,799 evaluable | 0.35% |
| Raw negative changed to final positive | 33,819 | 286,799 evaluable | 11.79% |
| Review burden | 49,632 | 301,602 opportunities | **16.46%** |

The disagreement is strongly directional: morphology usually rescues localized
signal whose per-object mean is below threshold. Therefore, the 12.14% is a
method-sensitivity rate, not evidence that morphology produces 12.14% errors.

All 286,799 evaluable calls used adaptive Otsu thresholds. Thus, 100% are
exploratory and 0% are confirmatory fixed-threshold calls.

## Marker-level breakdown

| Panel | Marker | Positive / negative / indeterminate | Discordance among evaluable | Review burden |
|---|---|---|---:|---:|
| E | CC10 | 41,331 / 49,372 / 4,847 | 19.25% | 23.34% |
| E | tdTOM | 24,968 / 65,735 / 4,847 | 11.04% | 15.55% |
| R | T1A | 6,462 / 28,669 / 1,703 | 6.43% | 10.76% |
| R | tdTOM | 13,836 / 21,295 / 1,703 | 11.72% | 15.81% |
| R | mRAGE | 2,734 / 32,397 / 1,703 | 2.79% | 7.28% |
| E | AcTub | 0 / 0 / 95,550 | N/A | N/A—regional-area endpoint |

## Gate outcomes

Across assignable marker-cell opportunities:

- insufficient spatial coverage: 198,029 occurrences (65.66%);
- fragmented spatial pattern: 136,338 occurrences (45.20%);
- shared perinuclear support: 14,803 occurrences (4.91%);
- invalid projection, wrong compartment, or enrichment failure: 0 in these
  particular panels.

Coverage and fragmentation reasons can co-occur and normally describe
morphology-negative cells. They are not rejection errors. Shared support is the
cause of the indeterminate calls in this data set.

## Consequences for script restructuring

1. Keep morphology as the final authority and retain mean intensity only as an
   audit comparator. A mean-intensity-only workflow would miss many localized
   positives, especially CC10.
2. Track nucleus rejection, marker indeterminacy, disagreement, and review
   burden separately. Do not collapse them into one percentage.
3. Treat panel-R G014 as a field-level DAPI QC outlier. Add image-level review
   before changing global nucleus-size settings.
4. Inspect rejected-object area distributions against manual annotations before
   lowering the minimum nucleus area. Nearly all rejected objects are below that
   cutoff, but the audit cannot determine which are debris versus real nuclei.
5. For acetylated tubulin, either supply an airway ROI before per-cell
   association or keep `cellCall=false` and report regional ciliary area. Its
   current 100% indeterminate status is expected, not an error.
6. Freeze control-derived intensity thresholds before confirmatory analysis;
   adaptive-Otsu calls should remain explicitly labeled exploratory.
7. Use the 16.46% canonical review burden as the baseline for regression tests.
   A lower number is desirable only if manually reviewed accuracy and biological
   morphology are preserved.

## Reproducible audit

Run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\TestRunErrorAudit.ps1
```

The audit scans all `run_summary.csv` and `*_cells.csv` files beneath
`test_runs`, removes the two `FinalPilot_*` duplicates from the canonical view,
and also reports all attempts for sensitivity analysis.

The main Fiji pipeline now exports candidate acceptance/rejection fractions,
marker positive/negative/indeterminate fractions, intensity–morphology
disagreement counts, and marker review-burden proxies directly in future
`run_summary.csv` files.
