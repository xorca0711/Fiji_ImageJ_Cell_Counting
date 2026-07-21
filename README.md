# Fiji/ImageJ confocal quantification — IFN-γ KO / PR8 influenza injury

A reproducible Fiji (ImageJ) pipeline for quantifying immunofluorescence
confocal sections from the **IFN-γ knockout / PR8 (H1N1) influenza lung-injury**
model. The headline readout is the **KRT5⁺ dysplastic "pod" area**, with
supporting alveolar (AT1/AT2), immune (CD4/CD8) and airway (Sox2) markers, plus
a future mechanistic panel (p63/YAP).

- **`IF_Quant_Pipeline.groovy`** — the analysis pipeline (run inside Fiji).
- **`aggregate_to_mouse.py`** — rolls per-region results up to the **mouse**
  (biological replicate) level and produces a stats-ready group summary.
- **`samplesheet_template.csv`** — per-image metadata template.

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

## 2. Panel / antibody design

Acquisition limit is **3 markers per slide = DAPI + 2 primaries**. Choose the
panel per slide via the samplesheet. Panel keys are single tokens so they also
survive filename parsing.

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
- **cyto** (KRT5, Pro-SPC, CD4, CD8) — measured in a perinuclear ring.
- **membrane** (AGER, PDPN) — measured in the ring; also interpretable as area.
- **nuc_marker** (p63, Sox2) — measured inside the nucleus.
- **nuc_ratio** (YAP) — nucleus vs. a true cytoplasmic ring → nuclear:cytoplasmic
  ratio (needs a single Z-plane; see caveats).

To change which acquisition channel is which marker, edit the `idx:` values in
the `PANELS` block at the top of the Groovy script.

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
   - Confirm the channel `idx:` order in `PANELS` matches your acquisition.
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
$env:IFQ_RECURSIVE = 'true'
$env:IFQ_INCLUDE_REGEX = '.*CC10_488.*20x 2k_Cycle.*G001_0001\.oir$'
$env:IFQ_MAX_IMAGES = '1'            # 0 means all matching files
$env:IFQ_TISSUE_MODE = 'whole_field' # count disconnected DAPI objects across the field
```

The 260719-CW filename convention is recognized automatically: mouse/date,
condition, section, and panel M/E/R are inferred. A `samplesheet.csv` remains
authoritative when supplied. Recursive runs add a stable suffix for duplicate
basenames so `Cycle` and `Cycle_01` data cannot overwrite each other.
Output folders and files use the concise pattern
`<mouse>_<condition>_<panel>_<section>`; the complete acquisition filename is
retained in `run_manifest.json` and each `__params.json` file.

Panel-R QC images display DAPI/T1A/tdTOM/mRAGE as blue/green/red/white.
Orange marks the analysis ROI, cyan marks counted nuclei, violet marks DAPI
particle candidates rejected by the minimum-area or image-edge rule, and
positive measurement rings use the corresponding marker color.

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
| `TISSUE_THRESH_METHOD` / `TISSUE_MIN_AREA_UM2` | auto tissue detection when no manual ROI |

**Positivity rule:** an object is positive for a marker if its mean raw
intensity ≥ (Otsu threshold computed **within that tissue region**) ×
`POS_SENSITIVITY[marker]`. Thresholds adapt per image/region; the *resolved*
numeric cutoff is exported so you can audit it.

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
    └── <stem>__KRT5_pod_mask.tif       # binary pod mask (255 = KRT5⁺ pod)
```

**Key `run_summary.csv` columns** (marker columns vary by panel):

- `mouse_id, section_id, genotype, condition, panel, region` — identifiers.
- `region_area_um2, n_nuclei` — denominators.
- `<marker>_pos_count`, `<marker>_density_per_mm2`, `<marker>_pos_threshold`.
- `KRT5_pod_area_um2`, `KRT5_pod_area_frac`, `KRT5_n_pods`,
  `KRT5_mean_pod_area_um2`, `KRT5_pod_threshold` — **primary endpoint**.
- `class_<rule>_count` — double +/− populations (e.g. `class_KRT5+_AGER-_count`).

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
- **Auto thresholds are placeholders.** They adapt per image but you must
  confirm positivity against QC overlays and freeze sensitivities before
  reporting. The automatic cutoff comes from tissue pixels but is applied to
  per-object ring/nuclear means; those distributions differ, so a sensitivity
  of `1.00` is not inherently a biologically neutral cutoff.
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
