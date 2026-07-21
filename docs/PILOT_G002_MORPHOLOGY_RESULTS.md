# G002 Morphology-Primary Fiji Pilot Results

Date: 2026-07-21

This report records one original 20x field from each requested staining panel,
processed with the morphology-first three-state hierarchy in
[`../IF_Quant_Pipeline.groovy`](../IF_Quant_Pipeline.groovy).

These are technical pilot results, not study endpoints. All marker and area
thresholds were image-specific Otsu thresholds, so every positive/negative result
is exploratory until thresholds are derived from controls and frozen.

## 1. CC10 / tdTomato / acetylated tubulin

Input field: `A01_G002_0001.oir`, panel E, 20x, one optical section.

Output folder: `test_runs/current/FinalPilot_CC10_AcTub_G002_morphology_primary_v2`

The field contains mixed airway-like and parenchymal structures. It was not
forced into a single airway compartment. This preserves AcTub cell calls as
indeterminate while retaining a regional ciliary-signal measurement.

| Endpoint | Result |
|---|---:|
| Analysis area | 405,000.0 um2 |
| Included nuclei | 2,909 |
| CC10 morphology positive | 1,431 (49.2%) |
| CC10 morphology negative | 1,340 (46.1%) |
| CC10 indeterminate | 138 (4.7%) |
| CC10 raw mean-intensity positive, audit only | 905 |
| tdTomato morphology positive | 659 (22.7%) |
| tdTomato morphology negative | 2,112 (72.6%) |
| tdTomato indeterminate | 138 (4.7%) |
| tdTomato raw mean-intensity positive, audit only | 324 |
| tdTomato filtered positive area | 35,371.4 um2 (8.73%) |
| AcTub morphology positive/negative | 0 / 0 |
| AcTub indeterminate | 2,909 (100%) |
| AcTub filtered regional area | 25,171.6 um2 (6.22%) |
| AcTub accepted components | 1,071 |
| AcTub mean accepted component area | 23.50 um2 |

The difference between morphology-positive and raw mean-positive CC10/tdTomato
counts is expected: a localized connected signal can satisfy the morphology gate
while its mean is diluted over the full support ring. This is why mean intensity
is no longer the call authority.

The AcTub regional result remains `exploratory_compartment_unassigned`. At 20x,
it must not be presented as 1,071 ciliated cells or individual cilia. An airway
ROI and control-derived cutoff are required for cell-associated interpretation.

## 2. T1alpha / tdTomato / mRAGE

Input field: `A01_G002_0001.oir`, panel R, 20x, one optical section.

Output folder: `test_runs/current/FinalPilot_T1A_mRAGE_G002_morphology_primary_v2`

Prior raw/QC review showed alveolar-parenchymal architecture without a conducting
airway. The run therefore used an explicit, provenance-recorded whole-field
`alveolar` assignment with compartment mode required.

| Endpoint | Result |
|---|---:|
| Analysis area | 405,000.0 um2 |
| Included nuclei | 2,108 |
| T1alpha morphology positive | 429 (20.4%) |
| T1alpha morphology negative | 1,571 (74.5%) |
| T1alpha indeterminate | 108 (5.1%) |
| T1alpha raw mean-intensity positive, audit only | 321 |
| T1alpha filtered membrane area | 29,792.5 um2 (7.36%) |
| T1alpha accepted membrane components | 747 |
| tdTomato morphology positive | 1,183 (56.1%) |
| tdTomato morphology negative | 817 (38.8%) |
| tdTomato indeterminate | 108 (5.1%) |
| tdTomato filtered positive area | 73,498.7 um2 (18.15%) |
| mRAGE morphology positive | 113 (5.4%) |
| mRAGE morphology negative | 1,887 (89.5%) |
| mRAGE indeterminate | 108 (5.1%) |
| mRAGE raw mean-intensity positive, audit only | 79 |
| mRAGE filtered membrane area | 12,203.2 um2 (3.01%) |
| mRAGE accepted membrane components | 245 |

The per-cell and membrane-area endpoints are complementary and must not be
treated as interchangeable. Cell calls ask whether a DAPI-associated local ring
has sufficient connected membrane support. Area endpoints describe the fraction
of the alveolar field occupied by accepted thresholded membrane components.

## 3. Adjustments applied after the first confirmation run

- Conventional `<marker>_pos_count` and density fields now use the final
  morphology call. Raw mean-intensity counts have explicit `raw_mean` names.
- Positive/negative/indeterminate state totals are separately exported and must
  sum to the included nucleus count.
- Components below the declared physical minimum are removed before both area
  measurement and mask export. The old implementation filtered component counts
  but still included tiny components in total area.
- T1alpha, mRAGE, AGER, and PDPN now receive independent filtered membrane-area
  masks in addition to per-cell calls.
- Every non-DAPI marker receives a call-QC PNG: positive green, negative cyan,
  indeterminate magenta, analysis ROI orange.
- `IFQ_WHOLE_FIELD_COMPARTMENT` provides an explicit, provenance-recorded
  override for a visually reviewed homogeneous field. It must not be used for a
  mixed field.
- Mouse aggregation explicitly keeps morphology endpoints primary and raw
  mean-intensity results audit-only.

## 4. Validation checks

- Groovy script compilation: passed.
- Both Fiji runs: completed successfully with one image selected each.
- For every marker, positive + negative + indeterminate equaled included nuclei.
- Three call-QC PNGs and positive/indeterminate label masks were produced per
  panel.
- `aggregate_to_mouse.py` syntax check: passed.
- Mouse-level aggregation on both one-image summaries: passed.

## 5. Required next step before biological reporting

Run negative, secondary-only/unstained, and known positive controls to choose
fixed marker cutoffs. Freeze cutoffs, spatial fractions, connectedness, support
width, minimum component area, compartment ROIs, and DAPI segmentation settings.
Then process the full blinded cohort and aggregate to mouse level.
