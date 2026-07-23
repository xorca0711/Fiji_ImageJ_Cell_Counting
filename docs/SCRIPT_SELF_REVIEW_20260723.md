# Script self-review (2026-07-23)

This review covers `IF_Quant_Pipeline.groovy`, `aggregate_to_mouse.py`, the
samplesheet contract, and the test-run audit helper. It distinguishes software
defects from biological/model limitations. A QC review fraction is not a
measured false-positive or false-negative rate unless calls are compared with
independent blinded ground truth.

## High-impact defects corrected

| Defect | Why it mattered | Correction |
|---|---|---|
| Tissue detection silently fell back to the whole field | Empty background or failed segmentation could be analyzed as tissue | Automatic detection now fails closed; whole-field analysis must be requested explicitly |
| Segmentation on a manual/non-whole-field ROI could shift object coordinates | The duplicated image processor could inherit and crop to the ROI | Segmentation now remains in full-image coordinates and is restricted by the region mask |
| Different files could share the same output name | Results could be overwritten or mixed | A stable `output_key` is generated; collisions receive a relative-path hash |
| Duplicate basenames in recursive input were ambiguous | Samplesheet metadata could be assigned to the wrong acquisition | `relative_path` is supported and required for duplicate basenames |
| Existing output directories were reused silently | Old and new results could be combined | Non-empty output is rejected by default; override requires `IFQ_ALLOW_NONEMPTY_OUTPUT=true` |
| Image/panel/channel failures could be skipped without a complete record | A batch could appear successful while losing samples | The manifest records `complete`, `partial`, or `failed`, per-image failures, and success/failure counts; headless failure exits non-zero |
| Multi-series Bio-Formats files silently used series 0 | The analyzed series might not be the intended specimen | Multi-series inputs now fail and require an explicit series policy |
| Split/projected ImageJ images were not reliably closed | Large 20x batches could exhaust memory | Raw, split, projected, mask, and QC images are closed in `finally` cleanup |
| Manual ROIs could be empty, invalid, duplicated, or overlapping | Cells could be omitted or double-counted | ROI files, bounds, names, and pairwise overlap are validated |
| Wrong anatomical compartment was counted as marker-negative | Compartment mismatch inflated biological negatives | Such calls are now `indeterminate`; they are excluded from evaluable positive/negative denominators |
| Boolean options accepted arbitrary strings as false | Configuration typos could silently change behavior | Runtime booleans now accept only `true` or `false` |
| `IFQ_MORPHOLOGY_PRIMARY=false` was recorded but did not implement a valid alternate decision path | Provenance could imply intensity-only classification even though the algorithm remained morphology-based | `false` is now rejected; morphology is always authoritative for final calls |
| Pixel calibration assumptions were implicit | Physical areas and ring size could be wrong | Positive micrometre calibration and approximately square pixels are required |
| CSV parsing used simple comma splitting | Quoted filenames or fields could corrupt metadata and audits | Samplesheet and audit parsing now support quoted CSV and UTF-8 BOM |
| Mouse aggregation averaged regional percentages | Unequal region sizes biased pooled estimates | Fractions are recomputed from summed numerators and denominators |
| Regional area endpoints were dropped during aggregation | Ciliary, reporter, membrane, and other area results were lost | Generic area totals, fractions, component counts, and component areas are retained |
| Missing mouse IDs became pseudo-mice | Invalid metadata could enter group summaries | Missing IDs, inconsistent mouse metadata, and duplicate region rows now fail validation |

## Verification performed

- Fiji headless smoke test on the established 20x G002 CC10/tdTOM/acetylated-tubulin acquisition completed with exit code 0.
- Counts matched the previous accepted pilot exactly: 2,909 included nuclei, 2,198 rejected nuclei, 1,431 CC10 morphology-positive cells, 659 tdTOM morphology-positive cells, and 2,909 AcTub calls held as indeterminate without an airway compartment annotation.
- A deliberate no-match run exited with code 1 and did not emit a misleading summary.
- Aggregation was checked on the 17-field test batch and on the new one-image smoke output.
- The C# audit helper re-read all 47 cell tables, including quoted CSV handling; 45 were canonical runs and two were retained historical attempts.

## Remaining limitations (not silently fixable)

### Requires blinded biological ground truth

1. **Pilot morphology cutoffs are not validated operating points.** Coverage,
   connected-component, enrichment, and ownership cutoffs must be frozen from
   control material and blinded annotations. The current review burden is a QC
   signal, not an error rate.
2. **Adaptive Otsu thresholds are exploratory.** Final cohort analysis should
   use acquisition-matched negative/positive controls and frozen thresholds.
3. **Disease and lineage identity require panels and context.** A single marker
   or morphology rule does not diagnose IPF, adenocarcinoma, or injury state.

### Known geometric/model constraints

1. The perinuclear ring is a nucleus-associated measurement proxy, not a full
   cell boundary. It is weakest for elongated, flat, multinucleated, or highly
   polarized cells.
2. Classic watershed has a minimum area gate but no validated maximum
   area/shape gate. Dense fields may retain fused nuclei or split irregular
   nuclei. Adding an arbitrary maximum would exchange one unmeasured error for
   another.
3. Ownership checks use neighboring-nucleus geometry and can be slow in very
   dense tissue. They do not replace membrane/cytoplasm segmentation.
4. Cells intersecting a manual ROI boundary follow the nucleus-based inclusion
   rule; partial cytoplasm near the boundary can still affect marker support.
5. Regional adaptive area masks use a field-derived threshold. A fixed
   control-derived area threshold is preferable when comparing cohorts.
6. Multi-series inputs intentionally fail; automated series selection or
   per-series expansion has not yet been implemented.
7. StarDist remains an optional, unvalidated path on this Windows ARM64 setup;
   the verified runtime uses classic watershed.
8. Isotropic ROI enlargement requires approximately square pixels. Anisotropic
   physical expansion is not implemented.
9. A generic `RoiSet.zip` can still be reused intentionally for multiple images
   in one directory. Prefer image-specific ROI filenames to avoid accidental
   reuse.

## Recommended validation before a final cohort run

1. Annotate stratified fields blindly, including weak signal, dense injury,
   airway, alveolar, tumor, fibrotic, and negative-control regions.
2. Evaluate nucleus detection and each marker call separately using confusion
   matrices and object-matching rules decided before tuning.
3. Freeze thresholds and morphology parameters on a training subset; report
   sensitivity, specificity, precision, recall, and indeterminate fraction on
   a held-out subset.
4. Lock acquisition settings, panel configuration, script commit, Fiji/plugin
   versions, and parameter JSON before running the study cohort.
5. Review QC overlays and manifest status for every batch; aggregate only
   complete, metadata-valid outputs.
