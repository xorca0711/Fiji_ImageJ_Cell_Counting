# Universal Lung Marker Configuration

This document turns marker selection into a reusable research workflow for
different microscopes, panels, species, and lung-disease questions. The
pipeline measures declared staining patterns; it does not infer a diagnosis or
cell identity from one marker.

The marker catalog is in
[`config/lung_marker_registry.json`](../config/lung_marker_registry.json). A
working multi-context configuration is in
[`config/custom_panels.example.json`](../config/custom_panels.example.json).

## 1. Marker-selection hierarchy

Use this order before creating a panel:

1. **Research question:** lineage abundance, transient state, localization,
   regional burden, or cell-cell geography.
2. **Species and preparation:** mouse versus human; fresh/frozen/FFPE; IF,
   immunohistochemistry, RNA-ISH, reporter, or FACS. A transcript marker is not
   automatically a validated IF antigen.
3. **Anatomical context:** airway, alveolar, tumor, fibrotic, stromal, vascular,
   or immune ROI. Draw context ROIs without viewing the target channel.
4. **Lineage anchor:** choose a relatively stable cell-compartment marker.
5. **State marker:** add activation, transition, differentiation, or disease
   state only after the lineage and geography are defined.
6. **Exclusion or alternative-lineage marker:** distinguish the nearest
   biologically plausible alternative.
7. **Analytical geometry:** nuclear, nuclear ratio, perinuclear cytoplasm,
   membrane, apical cilia, or regional area.
8. **Controls and frozen decision parameters:** validate antibody specificity,
   threshold, spatial coverage, connectedness, projection, and minimum area.

This hierarchy prevents common errors such as treating KRT8 alone as a DATP,
PDPN alone as AT1, NKX2-1 alone as lung adenocarcinoma, or ACTA2-positive area
alone as a pathogenic myofibroblast count.

## 2. Separate five kinds of information

| Layer | Question | Where it belongs |
|---|---|---|
| Marker identity | What antigen or reporter is measured? | Marker registry |
| Image layout | Which acquisition channel contains it? | Panel JSON `idx` |
| Analytical role | What spatial support is biologically appropriate? | Panel JSON `role`, or registry default |
| Biological context | In which anatomy/state can this interpretation be made? | ROI names, `expectedCompartments`, study design |
| Decision parameters | What control-derived values define the call? | Environment variables or panel overrides |

Do not encode a disease diagnosis as a marker role. For example, `KRT17` has a
cytoplasmic-keratin role whether it is measured in airway basal cells, an IPF
region, or a tumor. The context and co-markers change the interpretation.

## 3. Supported analytical roles

| Role | Primary support | Suitable examples | Main limitation |
|---|---|---|---|
| `nuclear` | DAPI/Hoechst nucleus | DAPI | Segmentation only |
| `nuc_marker` | Connected DAPI-nuclear signal plus enrichment | TP63, SOX9, Ki-67/MKI67, FOXJ1, NKX2-1, FOXP3 | Projection and nuclear segmentation |
| `nuc_ratio` | Nuclear support plus cytoplasmic reference | YAP | Single plane or validated 3D workflow |
| `cyto` | Nucleus-associated cytoplasmic ring | KRT5, KRT8, SFTPC, SCGB1A1, NAPSA | Approximate cell boundary |
| `membrane` | Local thin/circumferential support | AGER, PDPN, CLDN4, EPCAM, ITGA2, PDGFRB, CD31 | Shared boundaries and long stromal processes can be indeterminate |
| `apical_cilia` | Apical ciliary patch/support | Acetylated tubulin, TUBB4A | Regional patch area is primary at modest resolution |
| `regional_area` | Independent connected positive-area mask | COL1A1, CTHRC1, ACTA2, PDGFRB at 20x, secreted mucins | Not a per-cell identity call |

Unknown markers can be used without editing Groovy source when the panel JSON
declares one of these roles. Registry membership is a convenience, not a
whitelist.

### Operational profiles for KRT8, ITGA2, PDGFRB, SOX9, Red2-Kras, pan-KRAS, and Ki-67

| Marker | Default role | Primary morphology | Interpretation boundary |
|---|---|---|---|
| KRT8 | `cyto` | Connected extranuclear/perinuclear filament support | Positivity is broad epithelial expression; a KRT8-high transitional state additionally needs a frozen abundance rule, alveolar/fibrotic geography, and co-markers |
| ITGA2/CD49b | `membrane` | Connected membrane-support coverage; assay-dependent cytoplasmic staining is retained only as an audit measurement | Not lineage-specific; require the intended epithelial, tumor, stromal, or immune ROI and co-markers |
| PDGFRB/CD140b | `membrane` when cell ownership is validated; otherwise `regional_area` | Perivascular/stromal connected signal and spatial relationship to vessels | PDGFRB alone does not establish pericyte or myofibroblast identity; regional area is preferred at 20x |
| SOX9 | `nuc_marker` | Connected DAPI-nuclear support plus nuclear:ring enrichment | Cytoplasmic-only staining is not positive; developmental, injury, fibrotic, or tumor meaning requires geography and co-markers |
| Red2-KrasG12D RFP | `cyto` plus `areaMarker: true` in the study panel | Connected RFP reporter support; filtered RFP-positive clone area is primary at 20x | RFP-positive marks the oncogene-coupled clone after model/induction verification; RFP-negative alone is not wild type |
| KRAS | `cyto` | Connected nucleus-associated cytoplasmic/inner-membrane protein support | Pan-KRAS staining cannot establish a mutation or allele; mutant-specific claims require allele-specific validation against genotyping |
| Ki-67/MKI67 | `nuc_marker` | Connected DAPI-nuclear support plus enrichment, using a lower coverage gate for granular/nucleolar patterns | Report a labeling index within a predeclared cell population/ROI; there is no universal lung high/low cutoff |

The canonical symbol is `ITGA2`; the registry accepts the user-supplied
`IGTA2` spelling only as an alias. The explicit pilot morphology defaults are
documented in [`MARKER_MORPHOLOGY_GUIDE.md`](MARKER_MORPHOLOGY_GUIDE.md).

## 4. Research-context profiles

### Acute respiratory injury and regeneration

Use an alveolar lineage anchor together with transition and terminal-state
markers. KRT8-high/CLDN4-positive cells can support a DATP/ADI/PATS-like state,
but KRT8 is broadly epithelial and is not sufficient alone. Useful questions
include:

- `SFTPC -> KRT8/CLDN4 -> AGER/PDPN` progression;
- airway-derived `KRT5/TP63` expansion into damaged alveolar regions;
- recovery of `SCGB1A1` secretory and `FOXJ1/AcTub` ciliated programs;
- regional persistence versus resolution across explicitly recorded timepoints.

The injury-state basis comes from AT2-lineage and repair studies describing
damage-associated or transitional states, including
[Choi et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7487779/),
[Strunz et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7366678/), and
[Kobayashi et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7461628/).

### IPF and fibrotic remodeling

Separate epithelial state from stromal burden:

- aberrant basaloid research: `KRT17+ / KRT5-`, supported by TP63, ITGB6,
  MMP7, senescence markers, and fibrotic alveolar geography;
- persistent transitional epithelium: `KRT8-high / CLDN4+` with altered
  SFTPC/AGER/PDPN context;
- fibrotic stroma: CTHRC1, ACTA2, COL1A1, PDGFRA/PDGFRB, measured as regional
  area or carefully segmented stromal objects;
- AT1/AT2 loss and remodeling: quantify AGER/PDPN/AQP5 and SFTPC/ABCA3 in
  separate lineage-aware endpoints.

Human pulmonary-fibrosis atlases identify KRT5-negative/KRT17-positive
pathologic epithelial populations and distinct mesenchymal contributions; see
[Habermann et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7439444/) and
[Adams et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7439502/). These are
population/state findings, not stand-alone diagnostic IF rules.

### Lung adenocarcinoma research

The attached lineage reference provides normal lung-lineage context but is not
a tumor-diagnostic manual. For tumor-oriented research:

- define tumor ROIs using histomorphology or an orthogonal annotation;
- use a pulmonary glandular-lineage combination such as nuclear NKX2-1/TTF-1
  plus granular cytoplasmic Napsin A, optionally with EPCAM or KRT7;
- retain adjacent normal AT2 tissue as a distinct region because normal AT2
  cells can express both NKX2-1 and Napsin A;
- include the relevant differential lineage panel when the scientific question
  requires it, rather than converting a positive fluorescence pair into a
  diagnosis.

Primary immunohistochemistry studies support complementary use of TTF-1 and
Napsin A while also documenting imperfect sensitivity and specificity; see
[Ye et al.](https://pubmed.ncbi.nlm.nih.gov/21464700/) and
[Gurda et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4417108/). This Fiji
pipeline is for quantitative research and cannot replace pathologist review,
histomorphology, or a validated clinical assay.

## 5. Custom panel JSON

Point the pipeline to a study-owned JSON file:

```powershell
$env:IFQ_MARKER_REGISTRY = "$PWD\config\lung_marker_registry.json"
$env:IFQ_PANEL_CONFIG = 'D:\study\panels.json'
$env:IFQ_PANEL = 'ARI_DATP'
```

Minimal panel:

```json
{
  "schema_version": "1.0.0",
  "panels": {
    "MY_PANEL": {
      "label": "Study-specific description",
      "channels": [
        {"idx": 1, "marker": "DAPI", "role": "nuclear"},
        {"idx": 2, "marker": "KRT8", "expectedCompartments": ["alveolar"]},
        {"idx": 4, "marker": "CLDN4", "expectedCompartments": ["alveolar"]}
      ],
      "classify": [
        {"KRT8": true, "CLDN4": true}
      ]
    }
  }
}
```

Channel indices may skip unused acquisition channels. The pipeline checks the
highest referenced index against the imported image.

If a marker is in the registry, `role` can be omitted. For an unregistered
marker, declare it explicitly. Useful optional channel fields are:

| Field | Meaning |
|---|---|
| `measurement` | Human-readable measurement model stored in provenance |
| `expectedCompartments` | Any accepted ROI context tags |
| `requiresSinglePlane` | Marks multi-plane projections indeterminate |
| `cellCall` | Set `false` to disable per-nucleus calls |
| `areaMarker` | Produce an independent positive-area mask |
| `areaMode` | `generic`, `membrane`, `ciliary`, `reporter`, or `pod` |
| `areaMinAreaUm2` | Required for `regional_area`; must be study-validated |
| `areaBlurSigmaPx` | Pre-threshold blur for the area mask |
| `minPositiveFraction` | Channel-specific support-fraction override |
| `minLargestComponentShare` | Channel-specific connected-pattern override |
| `requireOwnership` | Reject shared nucleus-associated support when true |
| `minNuclearEnrichment` | Nuclear marker enrichment override |
| `minNucCytoRatio` | Nuclear localization ratio override |
| `supportExpandUm` | Ciliary/support-zone distance override |

Custom panels cannot replace built-in panel keys. This prevents a study file
from silently changing a previously validated acquisition map.
Unknown default or samplesheet panel keys stop the run instead of silently
falling back to another channel map.

## 6. ROI context vocabulary

An ROI name may contain more than one tag, for example
`alveolar_fibrotic_01` or `tumor_stromal_02`. Recognized tags are:

- `airway` or `bronchial`;
- `alveolar`;
- `tumor`, `tumour`, or `luad`;
- `fibrotic`, `honeycomb`, or `uip`;
- `stromal` or `mesenchymal`;
- `vascular`, `vessel`, or `capillary`;
- `immune`, `inflammatory`, or `lymphoid`;
- `ambiguous`.

The complete set is exported as `region_tags`. `compartment` remains as a
single backward-compatible primary label. A marker with multiple
`expectedCompartments` is evaluable when at least one expected tag is present.

## 7. Control and validation checklist

Before a confirmatory batch:

1. Verify species, antigen, clone, fluorophore, fixation, retrieval, and known
   positive/negative tissue.
2. Confirm the Bio-Formats channel map and optical-sectioning strategy.
3. Draw blinded anatomical/context ROIs.
4. Decide whether the endpoint is per nucleus, localization, regional area, or
   spatial relationship.
5. Derive fixed thresholds from controls and set
   `IFQ_<NORMALIZED_MARKER>_THRESHOLD`.
6. Validate morphology and minimum-area gates on blinded images spanning the
   expected biological and technical range.
7. Freeze panel JSON, registry version, thresholds, projection, segmentation,
   ROI vocabulary, and exclusions.
8. Retain three-state calls: positive, evaluable negative, and indeterminate.
9. Aggregate accepted regions to mouse or patient level before inference.

## 8. Reference scope and caveats

The registry was structured from the user-supplied seven-page
`Lung_Cell_Lineages_Markers_Reference.pdf`, which synthesizes proximal and
distal epithelial, stromal, vascular, immune, and injury-transition markers.
Its central caveat remains authoritative: marker expression is graded,
species-dependent, spatially variable, and altered by disease. The
[human lung cell atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC7704697/) is a
useful cross-study reference, but transcript abundance still does not guarantee
an IF-compatible antibody pattern.
