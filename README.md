# Fiji/ImageJ Lung Immunofluorescence Quantification

A morphology-first Fiji/ImageJ pipeline for lung immunofluorescence across
different image layouts and research questions. The validated built-in panels
retain the **IFN-γ knockout / PR8 influenza-injury** workflow and KRT5 pod
readout. An opt-in marker registry and study-owned JSON panel map extend the
same engine to acute injury/regeneration, IPF/fibrosis, stromal, vascular,
immune, and lung-adenocarcinoma lineage research without hard-coding each new
antibody combination.

- **`IF_Quant_Pipeline.groovy`** — the analysis pipeline (run inside Fiji).
- **`aggregate_to_mouse.py`** — rolls per-region results up to the **mouse**
  (biological replicate) level and produces a stats-ready group summary.
- **`samplesheet_template.csv`** — per-image metadata template.
- **[`WORKFLOW.md`](WORKFLOW.md)** — current end-to-end operational sequence.
- **[`docs/`](docs/README.md)** — marker morphology and validated pilot results.
- **[`config/`](config/README.md)** — universal marker registry and custom-panel
  examples.
- **[`legacy/`](legacy/README.md)** — non-authoritative historical archive.

---

## 1. Biological rationale

After severe influenza (PR8/H1N1) injury, p63⁺/KRT5⁺ distal-airway progenitors
(DASCs / LNEPs) expand and migrate into denuded alveolar zones, forming ectopic
**KRT5⁺ "pods"** — a marker of *dysplastic*, non-resolving repair that competes
with functional AT1/AT2 alveolar regeneration.

- **Reference result:** IFN-γ **receptor** KO → decrease in pod region.
- **This study's question:** does IFN-γ **ligand** KO (Het / KO on a PR8
  background) also reduce the pod region, while tracking the immune (CD4/CD8)
  and regenerative (AT1/AT2) response?

The pipeline therefore treats **KRT5⁺ pod area / fraction of tissue** as the
primary endpoint and expresses everything else relative to tissue area so that
genotypes can be compared fairly.

> ⚠️ Global IFN-γ **ligand** KO alters injury severity differently from
> epithelial **receptor** KO. Keep an independent viral-clearance control
> (NP stain or qPCR) outside the image analysis — this pipeline does not
> correct for differences in viral load.

---

## 2. Built-in panels and universal marker configuration

The table below preserves the original study panels. Choose a built-in panel
through the samplesheet or `IFQ_PANEL`. The engine itself is not limited to a
fixed number of markers: a panel can map any available acquisition channels,
including non-contiguous channel indices, as long as exactly one channel is the
nuclear segmentation channel.

| Panel | Channels (acquisition order)      | Purpose                                | Key classification      | Scheme 1 slides |
|:-----:|-----------------------------------|----------------------------------------|-------------------------|:---------------:|
| **A** | DAPI · KRT5 · AGER                | Pod area + AT1 (RAGE) boundary         | KRT5⁺/AGER⁻             | ×3              |
| **B** | DAPI · KRT5 · Pro-SPC             | Regeneration readout (AT2)             | KRT5⁻/Pro-SPC⁺          | ×2              |
| **C** | DAPI · KRT5 · CD8                 | Cytotoxic T-cell infiltrate            | CD8⁺, KRT5⁺/CD8⁺        | ×1              |
| **D** | DAPI · KRT5 · CD4                 | Helper T-cell infiltrate               | CD4⁺, KRT5⁺/CD4⁺        | —               |
| **P** | DAPI · KRT5 · PDPN                | AT1 alt (T1-α/podoplanin)              | **KRT5⁺/PDPN⁻**         | —               |
| **S** | DAPI · KRT5 · Sox2                | Airway/epithelial (optional)           | KRT5⁺/Sox2⁺, KRT5⁺/Sox2⁻ | —             |
| **S2**| DAPI · KRT5 · p63 · YAP           | **Scheme 2 (future)** mechanistic      | KRT5⁺/p63⁺, KRT5⁺/YAP⁺  | if possible     |

Marker roles the pipeline understands:

- **nuclear** (DAPI) — segmentation channel.
- **cyto** (KRT5, Pro-SPC, CC10, tdTomato) — measured in a perinuclear ring.
- **membrane** (AGER, PDPN, T1A, mRAGE, CD4, CD8) — measured in a local
  membrane-support ring.
- **nuc_marker** (p63, Sox2) — measured inside the nucleus.
- **nuc_ratio** (YAP) — nucleus vs. a true cytoplasmic ring → nuclear:cytoplasmic
  ratio (needs a single Z-plane; see caveats).
- **apical_cilia** (acetylated alpha-tubulin) — thresholded apical-cilia
  patches plus an explicitly approximate nucleus-proximity association.
- **regional_area** (COL1A1, CTHRC1, ACTA2, mucins) — independent positive-area
  mask without a per-nucleus positive/negative call.

For a new study, copy
[`config/custom_panels.example.json`](config/custom_panels.example.json), edit
the channel indices and marker set, and load it with `IFQ_PANEL_CONFIG`. Marker
aliases, localization, lineage/state notes, and default analytical roles are in
[`config/lung_marker_registry.json`](config/lung_marker_registry.json). The
registry is not a whitelist or diagnostic classifier; an unknown marker can be
used by declaring its role explicitly.

See
[`docs/UNIVERSAL_MARKER_CONFIGURATION.md`](docs/UNIVERSAL_MARKER_CONFIGURATION.md)
for the IPF, acute-injury, lung-adenocarcinoma, ROI-tag, and control hierarchy.

See [`docs/MARKER_MORPHOLOGY_GUIDE.md`](docs/MARKER_MORPHOLOGY_GUIDE.md) for the complete
marker-by-marker measurement, compartment-gating, optical-sectioning, control,
and QC rationale.

The concise end-to-end operational sequence is in [`WORKFLOW.md`](WORKFLOW.md).
Validated pilot results are indexed under [`docs/`](docs/README.md), while
superseded scripts, figures, and local run outputs are described in
[`legacy/`](legacy/README.md).

### Morphology-first call authority

The final marker call has three states: positive, negative, or indeterminate.
An intensity cutoff first defines candidate positive pixels. A final positive
then requires the marker-specific minimum spatial coverage, connected pattern,
localization/enrichment rule, unique ownership where applicable, and anatomical
compartment. A final negative is assigned only when the marker is evaluable.
Missing compartment, invalid YAP projection, empty support, or shared support is
indeterminate rather than negative.

`<marker>_pos` is the legacy mean-intensity audit result. Classifications and
counts use `<marker>_final_call`; they do not use the raw mean-intensity result.
Adaptive Otsu calls are explicitly exploratory. Use fixed, control-derived
`IFQ_<MARKER>_THRESHOLD` values for confirmatory analysis.

For panel E, marker morphology determines the analytical unit. CC10/SCGB1A1 is
measured in perinuclear cytoplasm as a secretory-cell phenotype, tdTomato is
measured in perinuclear cytoplasm plus a reporter-positive area mask, and
acetylated alpha-tubulin is measured as ciliary patches plus a 6-µm
nucleus-adjacent support zone. At 20x the latter is an association with nearby
cilia, not proof that an individual nucleus owns a resolved axoneme.

Do not edit built-in `PANELS` for a new acquisition. Put the new `idx:` mapping
in a study-owned JSON file so the channel map is versioned independently from
the analysis engine.

---

## 3. Feature coverage

Every requested feature maps to the pipeline:

| Requested feature | Where |
|---|---|
| Import original confocal files via Bio-Formats | `bfOpen()` — metadata + calibration preserved, no autoscale |
| Separate channels, preserve Z-stack / calibration | `ChannelSplitter.split` + `projectChannel` (calibration re-applied) |
| Define lung-tissue / lesion ROIs | `resolveTissueRois()` — manual `RoiSet.zip`/`.roi` if present, else auto from DAPI |
| Segment nuclei and cells consistently | `segmentNuclei()` — StarDist (preferred) or classic watershed fallback; perinuclear ring = "cell" |
| KRT5⁺ pod **area** and cell counts | independent threshold mask → `positiveAreaInRoi` + Analyze Particles; per-cell KRT5⁺ calls |
| AGER / PDPN and Pro-SPC populations | ring/membrane measurement + per-region adaptive threshold |
| Count CD4⁺ and CD8⁺ cells | panels C/D → per-cell positivity + density per mm² |
| Double-positive/negative (e.g. KRT5⁺/PDPN⁻) | `classify` rules per panel |
| Export per-cell, per-image, masks, QC overlays | `*__cells.csv`, `run_summary.csv`, label/pod masks (TIFF), QC PNGs |
| Record every threshold, filter, plugin version, parameter | `*__params.json` + `run_manifest.json`; resolved thresholds in `run_summary.csv` |

---

## 4. Requirements

- **[Fiji](https://fiji.sc/)** (ImageJ distribution). Bio-Formats is bundled.
- **StarDist + CSBDeep** update sites (recommended, for robust nuclei):
  `Help ▸ Update… ▸ Manage update sites` → tick **CSBDeep** and **StarDist** →
  apply → restart. If you cannot install them, set `SEGMENTER = "classic"` and
  the pipeline still runs with a watershed fallback.
- **Python 3** (standard library only) for `aggregate_to_mouse.py`.

---

## 5. Quick start

1. **Open the script in Fiji:** `File ▸ New ▸ Script…`, set **Language ▸ Groovy**,
   open `IF_Quant_Pipeline.groovy`.
2. **Edit the config block** (top of the script):
   - `INPUT_DIR` — folder of confocal files (`.czi/.lif/.nd2/.oib/.tif…`).
   - `OUTPUT_DIR` — where results are written.
   - `PANEL` — default panel if a file has no samplesheet/filename hint.
   - For a new acquisition, set `IFQ_PANEL_CONFIG` to a study-owned JSON file
     and confirm its channel `idx` map against Bio-Formats metadata.
3. **(Recommended) Add metadata.** Copy `samplesheet_template.csv` into
   `INPUT_DIR`, rename to **`samplesheet.csv`**, and fill one row per image.
   This is what carries `mouse_id` through every export — essential for correct
   statistics (see §8).
4. **(Optional) Draw lesion ROIs.** Save a `RoiSet.zip` (or `<image>.roi`) next
   to each image. Name ROIs `CTRL`/`KO` (or any label) to split a slide that
   carries two tissues. Without these the pipeline auto-detects tissue from DAPI.
5. **Run** (▶). Watch the Log window.
6. **Tune once, then batch.** Open a QC overlay, adjust per-marker
   `POS_SENSITIVITY` and pod settings (see §7), **freeze** the parameters, then
   run the whole set.
7. **Aggregate to mouse level:**
   ```bash
   python3 aggregate_to_mouse.py /path/to/analysis_output/run_summary.csv
   ```

For a headless run, use `SEGMENTER = "classic"`. On installations where the
Fiji launcher itself does not start (observed with one Windows ARM64 bundle),
invoke Fiji's bundled Java directly from PowerShell:

```powershell
$fiji = "X:\Fiji"
$java = Get-ChildItem "$fiji\java" -Recurse -Filter java.exe | Select-Object -First 1
& $java.FullName `
  '--add-opens=java.base/java.lang=ALL-UNNAMED' `
  "-javaagent:$fiji\jars\ij1-patcher-2.0.0.jar=init" `
  '-Djava.awt.headless=true' "-Dplugins.dir=$fiji" `
  -cp "$fiji\jars\*;$fiji\plugins\*" net.imagej.Main --headless `
  --run 'C:\path\to\IF_Quant_Pipeline.groovy'
```

The script exits automatically after all headless exports complete. StarDist's
ROI Manager output remains interactive-only; use the classic segmenter for
unattended jobs.

Paths and test selection can be supplied without editing the Groovy file:

```powershell
$env:IFQ_INPUT_DIR = 'G:\내 드라이브\260719-CW'
$env:IFQ_OUTPUT_DIR = "$PWD\test_runs\260719-CW_CC10_smoke"
$env:IFQ_PANEL = 'E'                 # M=4x DAPI/CC10/tdTOM; E/R=20x panels
$env:IFQ_MARKER_REGISTRY = "$PWD\config\lung_marker_registry.json"
# New studies only: $env:IFQ_PANEL_CONFIG = 'D:\study\panels.json'
$env:IFQ_RECURSIVE = 'true'
$env:IFQ_INCLUDE_REGEX = '.*CC10_488.*20x 2k_Cycle.*G001_0001\.oir$'
$env:IFQ_MAX_IMAGES = '1'            # 0 means all matching files
$env:IFQ_TISSUE_MODE = 'whole_field' # count disconnected DAPI objects across the field
$env:IFQ_COMPARTMENT_MODE = 'required' # final run: require alveoli/airway ROI names
$env:IFQ_T1A_THRESHOLD = '887.2'        # example only; freeze from controls
$env:IFQ_MRAGE_THRESHOLD = '503.1'      # example only; freeze from controls
$env:IFQ_T1A_MIN_RING_FRACTION = '0.30' # pilot value; validate before freezing
$env:IFQ_MRAGE_MIN_RING_FRACTION = '0.30'
$env:IFQ_DAPI_METHOD = 'local_phansalkar' # uneven-illumination alternative
$env:IFQ_DAPI_BACKGROUND_RADIUS_UM = '15'
$env:IFQ_DAPI_LOCAL_RADIUS_UM = '4'
$env:IFQ_DAPI_BLUR_SIGMA_PX = '1'
$env:IFQ_MIN_NUCLEUS_AREA_UM2 = '8' # pilot sensitivity; freeze after manual QC
$env:IFQ_ACTUB_SUPPORT_EXPAND_UM = '6'       # nucleus-to-apical-cilia support zone
$env:IFQ_ACTUB_MIN_SUPPORT_FRACTION = '0.02' # exploratory; freeze from controls
$env:IFQ_ACTUB_MIN_PATCH_AREA_UM2 = '0.5'    # 20x ciliary patch, not one axoneme
```

The 260719-CW filename convention is recognized automatically: mouse/date,
condition, section, and panel M/E/R are inferred. A `samplesheet.csv` remains
authoritative when supplied. Recursive runs add a stable suffix for duplicate
basenames so `Cycle` and `Cycle_01` data cannot overwrite each other.
Output folders and files use the concise pattern
`<mouse>_<condition>_<panel>_<section>`; the complete acquisition filename is
retained in `run_manifest.json` and each `__params.json` file.

Panel-R QC images display DAPI/T1A/tdTOM/mRAGE as blue/green/red/white.
Cyan is the only per-object outline and marks counted DAPI nuclei. Orange marks
the analysis ROI; green/red/white boundaries mark continuous fluorescence
regions for T1A/tdTOM/mRAGE. Rejected DAPI candidates remain available in the
separate `__rejected_nuclei_mask.tif` audit image.

For DAPI tuning, the default `local_phansalkar` path adds rolling-background
removal, gentle Gaussian smoothing, contrast normalization, local Phansalkar
thresholding, hole filling, and watershed; `global_otsu` preserves the original
fallback. Each run exports a display-balanced DAPI-only QC PNG and the actual
binary `__DAPI_candidate_mask.tif`. Display color balance does not alter the
mask; analytical preprocessing parameters are recorded in `__params.json`.

For morphology-aware final calls, draw compartment ROIs before looking at the
target-channel quantification and name each ROI `alveoli`, `airway`, or
`ambiguous` (suffixes such as `alveoli_01` are allowed). Put them in the usual
per-image ROI ZIP and set `IFQ_COMPARTMENT_MODE=required`. For T1A and mRAGE the
pipeline then exports three distinct decisions: intensity above a predeclared
cutoff, fraction of the perinuclear ring above that cutoff, and compartment
consistency. `<marker>_true_pos` is assigned only when all three pass; it is
left blank when morphology is unassigned. Thresholds and minimum ring fractions
must be selected from negative/positive controls and frozen before the study
groups are compared; the example values above are not validated cutoffs.

---

## 6. Input data expectations

- **Formats:** anything Bio-Formats reads (`.czi`, `.lif`, `.nd2`, `.oir`, `.oib`,
  `.oif`, `.ics`, `.tif/.tiff`). Adjust `FILE_GLOB` to widen/narrow.
- **Calibration:** must be embedded in the file (µm/pixel, Z-step). The pipeline
  reads and preserves it; all areas/distances are reported in µm/µm².
- **Metadata resolution order** for each file:
  1. `samplesheet.csv` row (matched by exact filename) — **preferred**;
  2. otherwise a filename token convention `mouseID_condition_panel_section.ext`
     (e.g. `m01_PR8_A_s1.czi`);
  3. otherwise the script defaults (`PANEL`, `genotype = NA`).

---

## 7. Configuration reference (most-used knobs)

| Parameter | Meaning |
|---|---|
| `SEGMENTER` | `"stardist"` (preferred) or `"classic"` watershed |
| `STARDIST_PROB` / `STARDIST_NMS` | detection probability / overlap thresholds |
| `PROJECTION` | `"max"` (default) · `"sum"` · `"avg"` · `"single"` |
| `SINGLE_PLANE` | plane index when `PROJECTION="single"` (−1 = middle) |
| `RING_EXPAND_UM` | perinuclear ring width → cytoplasm/membrane "cell" proxy |
| `POD_MIN_AREA_UM2` | minimum particle size to count as a KRT5⁺ pod |
| `POD_THRESH_METHOD` | auto-threshold method for the pod mask (`Otsu`, `Li`, …) |
| `POS_SENSITIVITY` | per-marker multiplier on the auto threshold (`>1` stricter, `<1` more permissive) |
| `IFQ_MORPHOLOGY_PRIMARY` | `true` by default; makes the three-state morphology call authoritative |
| `IFQ_WHOLE_FIELD_COMPARTMENT` | explicit compartment for a visually reviewed homogeneous field: `airway`, `alveolar`, `ambiguous`, or `unassigned` |
| `IFQ_<MARKER>_THRESHOLD` | fixed control-derived candidate-pixel cutoff for a marker |
| `IFQ_<MARKER>_MIN_POSITIVE_FRACTION` | minimum fraction of the role-specific support above cutoff |
| `IFQ_<MARKER>_MIN_LARGEST_COMPONENT_SHARE` | minimum connectedness of positive support |
| `IFQ_ACTUB_SUPPORT_EXPAND_UM` | radius beyond each nucleus used for AcTub proximity association |
| `IFQ_ACTUB_MIN_SUPPORT_FRACTION` | minimum fraction of that support zone above the AcTub threshold |
| `IFQ_ACTUB_MIN_PATCH_AREA_UM2` | minimum connected AcTub ciliary-patch area at 20x |
| `TISSUE_THRESH_METHOD` / `TISSUE_MIN_AREA_UM2` | auto tissue detection when no manual ROI |

**Positivity rule:** the resolved intensity cutoff is applied to pixels, and the
marker-specific morphology gates authorize the final call. Object mean intensity
is exported for audit but is not final-call authority. AcTub additionally
requires an airway ROI for a per-cell interpretation; its primary 20x outputs
remain regional ciliary area and connected patches. CC10 positivity denotes
protein phenotype and must not be used alone to infer club-cell lineage after
injury.

---

## 8. Outputs

Written under `OUTPUT_DIR`:

```
OUTPUT_DIR/
├── run_summary.csv                     # per-region summary (all images) — feed to aggregate_to_mouse.py
├── run_manifest.json                   # run-level provenance: versions + full config + per-image status
└── <image_stem>/
    ├── <stem>__cells.csv               # per-cell measurements (one row per nucleus/cell)
    ├── <stem>__params.json             # per-image parameters, calibration, channel map, thresholds
    ├── <stem>__<region>__QC.png        # QC overlay: tissue (white), nuclei (cyan), KRT5⁺ (magenta), pods (yellow)
    ├── <stem>__<region>__nuclei_mask.tif  # 16-bit label mask of nuclei
    ├── <stem>__<region>__<marker>_morphology_positive_nuclei_mask.tif
    ├── <stem>__<region>__<marker>_indeterminate_nuclei_mask.tif
    ├── <stem>__<region>__<marker>_CALL_QC.png  # final positive/negative/indeterminate overlay
    ├── <stem>__KRT5_pod_mask.tif       # binary pod mask (255 = KRT5⁺ pod)
    ├── <stem>__tdTOM_reporter_positive_mask.tif
    ├── <stem>__AcTub_ciliary_mask.tif  # area-filtered ciliary patches
    └── <stem>__T1A_membrane_positive_mask.tif  # analogous AGER/PDPN/mRAGE masks
```

`<stem>` is a concise channel signature so every exported file is
self-describing, for example
`C1-DAPI_C2-T1alpha-488_C3-tdTOM_C4-mRAGE-647__cells.csv`. The containing
directory carries the specimen ID, avoiding the long duplicated paths that
Windows may fail to open. The same signature is stored in `__params.json` and
`run_manifest.json`. For panels without a recorded fluorophore, the channel's
biological marker name is used.

**Key `run_summary.csv` columns** (marker columns vary by panel):

- `mouse_id, section_id, genotype, condition, panel, region` — identifiers.
- `region_area_um2, n_nuclei` — denominators.
- `<marker>_pos_count`, `<marker>_density_per_mm2` — morphology-authoritative
  positives; `<marker>_raw_mean_pos_count` preserves the legacy intensity audit.
- `<marker>_pos_threshold` and `<marker>_threshold_source` — resolved cutoff and
  whether it was fixed or adaptive.
- `<marker>_morphology_pos_count`, `<marker>_morphology_negative_count`,
  `<marker>_indeterminate_count`, `<marker>_morphology_evaluable_count` — the
  authoritative three-state decision summary.
- `KRT5_pod_area_um2`, `KRT5_pod_area_frac`, `KRT5_n_pods`,
  `KRT5_mean_pod_area_um2`, `KRT5_pod_threshold` — **primary endpoint**.
- `<marker>_positive_area_um2`, `<marker>_positive_area_frac`,
  `<marker>_n_components`, `<marker>_area_mode` — morphology-neutral regional
  fields; for panel E the modes are `reporter` (tdTomato) and `ciliary` (AcTub).
  Components below the declared physical minimum are excluded from the saved
  mask and from the area endpoint.
- `class_<rule>_count`, `class_<rule>_evaluable_count`, and
  `class_<rule>_indeterminate_count` — classifications consume only final calls.

After `aggregate_to_mouse.py` you also get:

- **`mouse_level_summary.csv`** — one row per (mouse × panel); pod fraction and
  densities are **area-weighted** across that animal's sections
  (`Σ pod_area / Σ tissue_area`, `Σ counts / Σ area`).
- **`group_level_summary.csv`** — mean / sd / **sem** / **n_mice** per
  (genotype × condition × panel × metric). `n_mice` is the real statistical n.

---

## 9. Statistics — n is counted by MICE, not sections

This is the single most important analysis rule for this study. Three technical
sections from one animal are still **n = 1**. `mouse_id` and `section_id` travel
through every export precisely so you can collapse to the animal before testing.

`aggregate_to_mouse.py` does that collapse correctly (area-weighted pooling) and
reports `n_mice`. Run your comparison (e.g. KO vs Het vs WT pod fraction) on the
**mouse-level** values. The script intentionally does **not** compute p-values —
export `mouse_level_summary.csv` to Prism/R and apply the appropriate test
(e.g. Mann-Whitney / Kruskal-Wallis for small n, or a mixed model if you keep
sections as repeated measures).

---

## 10. Threshold-tuning workflow

1. Run 2–3 representative images with default settings.
2. Open the `__QC.png` overlays. Check: are nuclei split correctly? Are KRT5⁺
   pods (yellow) matching what you see by eye? Are AGER/PDPN/CD calls sensible?
3. Adjust:
   - **Under-segmented nuclei** → lower `STARDIST_PROB` / raise `STARDIST_NMS`.
   - **Too many/few positive cells** → tune `POS_SENSITIVITY[marker]`.
   - **Pod mask too greedy/sparse** → change `POD_THRESH_METHOD`, `POD_BLUR_SIGMA_PX`,
     `POD_MIN_AREA_UM2`.
4. **Freeze** the parameters and batch the whole cohort with identical settings.
   Never tune per-image across the cohort — that biases group comparisons.

---

## 11. Caveats & interpretation

- **AGER/PDPN are thin membrane signals.** Per-cell mean is weak; KRT5⁺/AGER⁻
  (or /PDPN⁻) is most robust as an **area** relationship (pod area vs. AT1 area).
  Interpret the per-object AT1 call conservatively.
- **Projection matters.** MAX projection is fine for pod **area**. For **YAP**
  (Scheme 2) a MAX projection corrupts the nuclear:cytoplasmic ratio — set
  `PROJECTION="single"` and use one representative plane.
- **Auto thresholds are exploratory placeholders.** They adapt per image and
  define candidate pixels for morphology gates. Confirm them against controls,
  replace them with fixed marker cutoffs, and freeze every morphology parameter
  before reporting confirmatory results.
- **Ligand vs. receptor KO** and **viral load** — see §1.

---

## 12. Reproducibility / provenance

Every run records, so results can be reconstructed exactly:

- ImageJ, Bio-Formats, Java and OS versions (`run_manifest.json`).
- Full configuration (segmenter + parameters, projection, `blackBackground`,
  tissue method, sensitivities).
- Per-image calibration, channel→marker map, and the **resolved numeric
  thresholds** for pods and each marker.
- StarDist/CSBDeep: model + parameters are captured; if you need exact plugin
  build strings, note them from `Help ▸ Update` (they are not reliably
  queryable from a script).

---

## 13. Troubleshooting

| Symptom | Fix |
|---|---|
| `INPUT_DIR is not a folder` | set `INPUT_DIR`/`OUTPUT_DIR` at the top of the script |
| `found N channels, panel expects M → skipped` | fix the panel's `idx:`/marker map or the samplesheet `panel` |
| StarDist errors / not found | install CSBDeep + StarDist update sites, or set `SEGMENTER="classic"` |
| All areas look inverted / zero | pipeline forces `blackBackground=true`; if you edited masking, keep foreground = 255 |
| Densities look 10–100× off | check embedded calibration (µm/pixel) in the source files |
