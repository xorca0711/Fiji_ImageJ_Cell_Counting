/*
 * ============================================================================
 *  IF_Quant_Pipeline.groovy
 *  Morphology-first lung immunofluorescence quantification. The built-in
 *  panels preserve the IFN-gamma KO / PR8 influenza-injury KRT5-pod workflow;
 *  optional JSON panel maps support other lung research contexts.
 * ----------------------------------------------------------------------------
 *  Built-in antibody set: KRT5, Pro-SPC, AGER, PDPN, CD4, CD8, Sox2
 *  (+ p63, YAP, Aqp5, CC10, tdTomato, AcTub, T1A, mRAGE).
 *
 *  Built-in panels. PANEL keys are single tokens so they survive filename
 *  parsing; select per image via samplesheet. New panels are loaded through
 *  IFQ_PANEL_CONFIG and may map any available acquisition channels:
 *     A : DAPI | KRT5 | AGER      -> pod + AT1 boundary   (KRT5+/AGER-)   [Scheme1 x3]
 *     B : DAPI | KRT5 | Pro-SPC   -> regeneration readout (AT2)          [Scheme1 x2]
 *     C : DAPI | KRT5 | CD8       -> cytotoxic T infiltrate              [Scheme1 x1]
 *     D : DAPI | KRT5 | CD4       -> helper T infiltrate
 *     P : DAPI | KRT5 | PDPN      -> AT1 alt (T1-alpha)    (KRT5+/PDPN-)
 *     S : DAPI | KRT5 | Sox2      -> airway/epithelial (optional)
 *     S2: DAPI | KRT5 | p63 | YAP -> FUTURE Scheme 2 (mechanistic; see notes)
 *
 *  Feature coverage:
 *   [x] Bio-Formats import (metadata + calibration preserved)
 *   [x] Channel split, Z handling (projection configurable), calibration kept
 *   [x] Tissue / lesion ROIs (manual RoiSet.zip if present, else auto-tissue)
 *   [x] Consistent nucleus segmentation (StarDist, classic watershed fallback)
 *   [x] Perinuclear "cell" ring for cytoplasmic/membrane marker readout
 *   [x] KRT5+ pod AREA (independent threshold) + pod count/size distribution
 *   [x] KRT5+ cell counts, AGER/Pro-SPC/CD8 populations
 *   [x] Double +/- classification (e.g. KRT5+/AGER-)
 *   [x] Exports: per-cell CSV, per-image summary, masks, QC overlays
 *   [x] Full provenance: versions, every threshold/filter/parameter -> JSON
 * ----------------------------------------------------------------------------
 *  REQUIREMENTS
 *   - Fiji (Bio-Formats is bundled).
 *   - Optional but recommended update sites: CSBDeep + StarDist
 *     (Help > Update > Manage update sites). If absent, set
 *     SEGMENTER = "classic" below and the pipeline still runs.
 *   - Run from Fiji: File > New > Script..., language = Groovy, Run.
 *
 *  HOW TO USE
 *   1. Set INPUT_DIR / OUTPUT_DIR and PANEL below (or use a samplesheet).
 *   2. Set the channel ORDER for your acquisition in PANELS (1-based).
 *   3. (Optional) Put a RoiSet.zip / <image>.roi next to each image to use
 *      manually drawn lesion ROIs. Name ROIs "CTRL"/"KO" to split a slide
 *      that carries control and KO tissue side by side.
 *   4. Run once, open a QC overlay, then TUNE thresholds (see CAVEATS).
 * ----------------------------------------------------------------------------
 *  CAVEATS (read before trusting any number)
 *   * n IS COUNTED BY MICE, NOT SECTIONS. mouse_id/section_id travel through
 *     every export so you can aggregate to biological n downstream. Three
 *     technical sections from one animal are still n = 1.
 *   * DEFAULT THRESHOLDS AND MORPHOLOGY GATES ARE PILOT PLACEHOLDERS. Auto-Otsu
 *     adapts per image and therefore produces exploratory calls only. Derive
 *     fixed cutoffs from controls, validate the spatial gates, freeze all
 *     parameters once, then batch the study cohort.
 *   * AGER is a thin AT1 membrane signal; per-cell mean is weak. KRT5+/AGER-
 *     is most robust as an AREA relationship (pod area vs AT1 area). The
 *     per-object AGER call is provided but interpret it conservatively.
 *   * PROJECTION: default is MAX intensity (fine for pod AREA, standard).
 *     For Scheme 2 YAP the nuclear:cytoplasmic ratio is CORRUPTED by MIP —
 *     use a SINGLE representative plane (projection="single") for YAP.
 *   * Global IFN-g LIGAND KO changes injury severity differently from the
 *     epithelial RECEPTOR KO. This script does not correct for that; keep a
 *     viral-clearance control (NP stain / qPCR) outside the image analysis.
 * ============================================================================
 */

import ij.IJ
import ij.ImagePlus
import ij.Prefs
import ij.gui.Roi
import ij.gui.Overlay
import ij.gui.PolygonRoi
import ij.gui.ShapeRoi
import ij.measure.Calibration
import ij.measure.Measurements
import ij.measure.ResultsTable
import ij.process.ImageProcessor
import ij.process.ImageStatistics
import ij.process.AutoThresholder
import ij.process.ByteProcessor
import ij.process.ColorProcessor
import ij.process.ShortProcessor
import ij.plugin.ZProjector
import ij.plugin.ChannelSplitter
import ij.plugin.RoiEnlarger
import ij.plugin.filter.GaussianBlur
import ij.plugin.filter.ParticleAnalyzer
import ij.plugin.filter.ThresholdToSelection
import ij.io.RoiDecoder
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import loci.formats.FormatTools
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.awt.Color
import java.awt.Font
import java.awt.Rectangle
import java.util.zip.ZipInputStream

// ============================================================================
//  1. USER CONFIG
// ============================================================================

// Environment variables let unattended/test runs override these defaults
// without rewriting the analysis script (see README). Forward slashes are
// accepted by Java on Windows and avoid escaping drive paths.
def envOr = { String name, String fallback ->
  def value = System.getenv(name)
  return (value == null || value.trim().isEmpty()) ? fallback : value.trim()
}
def INPUT_DIR   = envOr("IFQ_INPUT_DIR", new File("ref_images").getAbsolutePath())
def OUTPUT_DIR  = envOr("IFQ_OUTPUT_DIR", new File("analysis_output").getAbsolutePath())
def PANEL       = envOr("IFQ_PANEL", "T")
// The registry is descriptive evidence and supplies safe role defaults. A
// separate panel JSON maps the actual acquisition channels for a study. This
// keeps biological identity, analytical geometry, and image layout independent.
def MARKER_REGISTRY_PATH = envOr("IFQ_MARKER_REGISTRY", new File("config/lung_marker_registry.json").getAbsolutePath())
def PANEL_CONFIG_PATH = envOr("IFQ_PANEL_CONFIG", "")
def FILE_GLOB   = ~/(?i).*\.(czi|lif|nd2|oir|oib|oif|ics|tif|tiff)$/
def RECURSIVE   = envOr("IFQ_RECURSIVE", "false").toBoolean()
def INCLUDE_REGEX = envOr("IFQ_INCLUDE_REGEX", ".*")
def MAX_IMAGES  = Integer.parseInt(envOr("IFQ_MAX_IMAGES", "0")) // 0 = all
def TISSUE_MODE = envOr("IFQ_TISSUE_MODE", "auto").toLowerCase() // auto | whole_field
def COMPARTMENT_MODE = envOr("IFQ_COMPARTMENT_MODE", "optional").toLowerCase() // optional | required
// Explicit override for a visually reviewed, anatomically homogeneous field.
// Never use this to force a mixed airway/alveolar image into one compartment.
def WHOLE_FIELD_COMPARTMENT = envOr("IFQ_WHOLE_FIELD_COMPARTMENT", "unassigned").toLowerCase()
def ALLOWED_COMPARTMENTS = ["unassigned", "airway", "alveolar", "tumor", "fibrotic",
                            "stromal", "vascular", "immune", "ambiguous"] as Set
if (!ALLOWED_COMPARTMENTS.contains(WHOLE_FIELD_COMPARTMENT)) {
  throw new IllegalArgumentException("IFQ_WHOLE_FIELD_COMPARTMENT must be one of " + ALLOWED_COMPARTMENTS)
}
def FIXED_POS_THRESHOLDS = [:]
// Any marker can use a control-derived fixed cutoff via IFQ_<MARKER>_THRESHOLD.
// Adaptive Otsu remains available for pilots but is explicitly reported as
// exploratory. Examples: IFQ_CC10_THRESHOLD, IFQ_TDTOM_THRESHOLD.
def thresholdMarkers = ["KRT5", "AGER", "PDPN", "ProSPC", "CD8", "CD4",
                        "Sox2", "Aqp5", "p63", "YAP", "CC10", "tdTOM",
                        "AcTub", "T1A", "mRAGE"]
thresholdMarkers.each { marker ->
  String token = marker.toUpperCase().replaceAll(/[^A-Z0-9]+/, "")
  def rawValue = System.getenv("IFQ_" + token + "_THRESHOLD")
  if (rawValue != null && !rawValue.trim().isEmpty()) {
    FIXED_POS_THRESHOLDS[marker] = Double.parseDouble(rawValue.trim())
  }
}
def MIN_RING_POS_FRACTION = [
  "T1A": Double.parseDouble(envOr("IFQ_T1A_MIN_RING_FRACTION", "0.30")),
  "mRAGE": Double.parseDouble(envOr("IFQ_MRAGE_MIN_RING_FRACTION", "0.30"))
]
// Acetylated alpha-tubulin is concentrated in apical motile cilia rather than
// the perinuclear cytoplasm. At 20x, individual axonemes are not reliably
// resolvable, so associate thresholded ciliary signal with a nucleus using a
// wider proximity support zone and report ciliary patches as the primary
// regional readout. These pilot defaults must still be frozen against controls.
def ACTUB_SUPPORT_EXPAND_UM = Double.parseDouble(envOr("IFQ_ACTUB_SUPPORT_EXPAND_UM", "6.0"))
def ACTUB_MIN_SUPPORT_FRACTION = Double.parseDouble(envOr("IFQ_ACTUB_MIN_SUPPORT_FRACTION", "0.10"))
def ACTUB_MIN_PATCH_AREA_UM2 = Double.parseDouble(envOr("IFQ_ACTUB_MIN_PATCH_AREA_UM2", "2.0"))

// Morphology is the authoritative marker call. Intensity thresholds define
// candidate pixels; a final positive additionally requires role-appropriate
// coverage, connected spatial support, compartment, and object ownership.
// Defaults are conservative PILOT values and must be calibrated from blinded
// positive/negative controls before a final cohort run.
def MORPHOLOGY_PRIMARY = envOr("IFQ_MORPHOLOGY_PRIMARY", "true").toBoolean()
// These are geometry-class pilot defaults, not biological truth. Marker-level
// validated values below take precedence, followed by explicit channel-level
// overrides from IFQ_PANEL_CONFIG. Unknown markers therefore never need to be
// added to source code merely to participate in the same measurement model.
def ROLE_MORPHOLOGY_DEFAULTS = [
  "cyto"        : [minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true],
  "membrane"    : [minFraction:0.25d, minLargestShare:0.40d, requireOwnership:true],
  "nuc_marker"  : [minFraction:0.40d, minLargestShare:0.60d, requireOwnership:false, minNuclearEnrichment:1.25d],
  "nuc_ratio"   : [minFraction:0.30d, minLargestShare:0.60d, requireOwnership:false, minNucCytoRatio:1.50d],
  "apical_cilia": [minFraction:ACTUB_MIN_SUPPORT_FRACTION, minLargestShare:0.30d, requireOwnership:true]
]
def MORPHOLOGY_RULES = [
  "KRT5" : [minFraction:0.20d, minLargestShare:0.50d, requireOwnership:true],
  "AGER" : [minFraction:0.25d, minLargestShare:0.40d, requireOwnership:true],
  "PDPN" : [minFraction:0.25d, minLargestShare:0.40d, requireOwnership:true],
  "ProSPC":[minFraction:0.15d, minLargestShare:0.40d, requireOwnership:true],
  "CD8"  : [minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true],
  "CD4"  : [minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true],
  "Sox2" : [minFraction:0.40d, minLargestShare:0.60d, requireOwnership:false, minNuclearEnrichment:1.25d],
  "Aqp5" : [minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true],
  "p63"  : [minFraction:0.40d, minLargestShare:0.60d, requireOwnership:false, minNuclearEnrichment:1.25d],
  "YAP"  : [minFraction:0.30d, minLargestShare:0.60d, requireOwnership:false, minNucCytoRatio:1.50d],
  "CC10" : [minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true],
  "tdTOM":[minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true],
  "AcTub":[minFraction:ACTUB_MIN_SUPPORT_FRACTION, minLargestShare:0.30d, requireOwnership:true],
  "T1A"  : [minFraction:MIN_RING_POS_FRACTION["T1A"], minLargestShare:0.40d, requireOwnership:true],
  "mRAGE":[minFraction:MIN_RING_POS_FRACTION["mRAGE"], minLargestShare:0.40d, requireOwnership:true]
]
MORPHOLOGY_RULES.each { marker, rule ->
  String token = marker.toUpperCase().replaceAll(/[^A-Z0-9]+/, "")
  rule.minFraction = Double.parseDouble(envOr("IFQ_" + token + "_MIN_POSITIVE_FRACTION", rule.minFraction.toString()))
  rule.minLargestShare = Double.parseDouble(envOr("IFQ_" + token + "_MIN_LARGEST_COMPONENT_SHARE", rule.minLargestShare.toString()))
  if (rule.containsKey("minNuclearEnrichment")) {
    rule.minNuclearEnrichment = Double.parseDouble(envOr("IFQ_" + token + "_MIN_NUCLEAR_ENRICHMENT", rule.minNuclearEnrichment.toString()))
  }
  if (rule.containsKey("minNucCytoRatio")) {
    rule.minNucCytoRatio = Double.parseDouble(envOr("IFQ_" + token + "_MIN_NUC_CYTO_RATIO", rule.minNucCytoRatio.toString()))
  }
}
def DAPI_METHOD = envOr("IFQ_DAPI_METHOD", "local_phansalkar").toLowerCase() // local_phansalkar | global_otsu
def DAPI_BACKGROUND_RADIUS_UM = Double.parseDouble(envOr("IFQ_DAPI_BACKGROUND_RADIUS_UM", "15.0"))
def DAPI_LOCAL_RADIUS_UM = Double.parseDouble(envOr("IFQ_DAPI_LOCAL_RADIUS_UM", "4.0"))
def DAPI_BLUR_SIGMA_PX = Double.parseDouble(envOr("IFQ_DAPI_BLUR_SIGMA_PX", "1.0"))
def DAPI_CONTRAST_SATURATION = Double.parseDouble(envOr("IFQ_DAPI_CONTRAST_SATURATION", "0.35"))

// If present, a samplesheet.csv in INPUT_DIR overrides per-file metadata.
// Columns: filename,mouse_id,section_id,genotype,condition,panel
def USE_SAMPLESHEET = true

// --- Segmentation ---
// PILOT: forced to "classic". StarDist/CSBDeep are not installed in this Fiji,
// and TensorFlow has no windows-arm64 native build, so "stardist" cannot run on
// this machine. Set back to "stardist" on an x86_64 Fiji with the update sites on.
def SEGMENTER   = "classic"     // "stardist" (preferred) | "classic"
def STARDIST_PROB = 0.50
def STARDIST_NMS  = 0.40
def STARDIST_TILES = 1          // raise (e.g. 4/9) for large images / low RAM

// --- Z handling ---
def PROJECTION  = "max"         // "max" | "sum" | "avg" | "single"
def SINGLE_PLANE = -1           // used only when PROJECTION="single"; -1 = middle

// --- Geometry (calibrated, micrometres) ---
def RING_EXPAND_UM      = 2.0   // perinuclear ring width -> "cell" proxy
def MIN_NUCLEUS_AREA_UM2 = Double.parseDouble(envOr("IFQ_MIN_NUCLEUS_AREA_UM2", "8.0"))
def POD_MIN_AREA_UM2    = 50.0  // a "pod" particle must exceed this
def POD_BLUR_SIGMA_PX   = 2.0
def POD_THRESH_METHOD   = "Otsu" // Otsu|Triangle|Li|Huang|MaxEntropy...

// --- Candidate-pixel intensity cutoffs ---
// Otsu(inTissue) * sensitivity supplies a pilot threshold for spatial-support
// measurements. The morphology gates above, not the object mean, authorize the
// final call. Fixed control-derived thresholds should replace Otsu for final use.
def POS_SENSITIVITY = [ "KRT5":1.00, "AGER":1.00, "PDPN":1.00, "ProSPC":1.00,
                        "CD8":1.00, "CD4":1.00, "Sox2":1.00, "Aqp5":1.00,
                        "p63":1.00, "YAP":1.00, "CC10":1.00, "tdTOM":1.00,
                        "AcTub":1.00, "T1A":1.00, "mRAGE":1.00 ]

// --- Tissue auto-detection (used only if no manual ROI is supplied) ---
def TISSUE_BLUR_SIGMA_PX = 4.0
def TISSUE_THRESH_METHOD  = "Triangle" // permissive, keeps sparse tissue
def TISSUE_MIN_AREA_UM2   = 2000.0     // drop debris specks

// ============================================================================
//  2. PANEL DEFINITIONS  (channel idx is 1-based ACQUISITION order)
//     role: "nuclear"   -> segmentation channel (DAPI)
//           "cyto"      -> measured in perinuclear ring (KRT5, Pro-SPC, CC10)
//           "membrane"  -> measured in ring (AGER/PDPN/CD4/CD8)
//           "nuc_marker"-> measured in the nucleus (p63)
//           "nuc_ratio" -> nucleus vs ring separately (YAP)  [needs single plane]
//     areaMarker: also run independent threshold AREA quantification (KRT5 pod)
// ============================================================================

def PANELS = [
  "A": [ label:"A_KRT5_AGER",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear"],
               [idx:2, marker:"KRT5",  role:"cyto",     areaMarker:true],
               [idx:3, marker:"AGER",  role:"membrane", expectedCompartment:"alveolar",
                areaMarker:true, areaMode:"membrane", areaMinAreaUm2:2.0d, areaBlurSigmaPx:0.7d] ],
    classify:[ ["KRT5":true,"AGER":false], ["KRT5":true,"AGER":true] ] ],

  "B": [ label:"B_KRT5_ProSPC",
    channels:[ [idx:1, marker:"DAPI",   role:"nuclear"],
               [idx:2, marker:"KRT5",   role:"cyto",    areaMarker:true],
               [idx:3, marker:"ProSPC", role:"cyto", expectedCompartment:"alveolar"] ],
    classify:[ ["KRT5":true,"ProSPC":false], ["KRT5":false,"ProSPC":true] ] ],

  "C": [ label:"C_KRT5_CD8",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:3, marker:"CD8",  role:"membrane"] ],
    classify:[ ["CD8":true], ["KRT5":true,"CD8":true] ] ],

  "D": [ label:"D_KRT5_CD4",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:3, marker:"CD4",  role:"membrane"] ],
    classify:[ ["CD4":true], ["KRT5":true,"CD4":true] ] ],

  // AT1 alternative via podoplanin (T1-alpha). Enables the KRT5+/PDPN- readout.
  "P": [ label:"P_KRT5_PDPN",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:3, marker:"PDPN", role:"membrane", expectedCompartment:"alveolar",
                areaMarker:true, areaMode:"membrane", areaMinAreaUm2:2.0d, areaBlurSigmaPx:0.7d] ],
    classify:[ ["KRT5":true,"PDPN":false], ["KRT5":true,"PDPN":true] ] ],

  // Optional airway/epithelial marker.
  "S": [ label:"S_KRT5_Sox2",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",       areaMarker:true],
               [idx:3, marker:"Sox2", role:"nuc_marker"] ],
    classify:[ ["Sox2":true], ["KRT5":true,"Sox2":true], ["KRT5":true,"Sox2":false] ] ],

  // 260719-CW Olympus OIR panels. Bio-Formats metadata confirms channel order
  // by emission bands: 470 (DAPI), 540 (488), 620 (tdTOM), 750 nm (647).
  // The 4x stitched mapping file contains only the first three channels.
  "M": [ label:"M_4x_CC10_tdTOM_mapping",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear", qcColor:"blue"],
               [idx:2, marker:"CC10",  role:"cyto",    measurement:"perinuclear_secretory_cytoplasm", qcColor:"green", fileLabel:"CC10-488"],
               [idx:3, marker:"tdTOM", role:"cyto",    measurement:"perinuclear_lineage_reporter", qcColor:"red",
                areaMarker:true, areaMode:"reporter", areaMinAreaUm2:8.0d] ],
    classify:[ ["tdTOM":true], ["CC10":true,"tdTOM":true] ] ],

  "E": [ label:"E_CC10_tdTOM_AcTub",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear", qcColor:"blue"],
               [idx:2, marker:"CC10",  role:"cyto",    measurement:"perinuclear_secretory_cytoplasm", qcColor:"green", fileLabel:"CC10-488"],
               [idx:3, marker:"tdTOM", role:"cyto",    measurement:"perinuclear_lineage_reporter", qcColor:"red",
                areaMarker:true, areaMode:"reporter", areaMinAreaUm2:8.0d],
               [idx:4, marker:"AcTub", role:"apical_cilia", measurement:"apical_cilia_proximity", expectedCompartment:"airway", qcColor:"white", fileLabel:"AcTub-647",
                areaMarker:true, areaMode:"ciliary", areaMinAreaUm2:ACTUB_MIN_PATCH_AREA_UM2, areaBlurSigmaPx:0.7d] ],
    classify:[ ["tdTOM":true], ["CC10":true,"tdTOM":true],
               ["AcTub":true,"tdTOM":true] ] ],

  "R": [ label:"R_T1A_tdTOM_mRAGE",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear",  qcColor:"blue"],
               [idx:2, marker:"T1A",   role:"membrane", expectedCompartment:"alveolar", qcColor:"green", fileLabel:"T1alpha-488",
                areaMarker:true, areaMode:"membrane", areaMinAreaUm2:2.0d, areaBlurSigmaPx:0.7d],
               [idx:3, marker:"tdTOM", role:"cyto", measurement:"perinuclear_lineage_reporter", qcColor:"red",
                areaMarker:true, areaMode:"reporter", areaMinAreaUm2:8.0d],
               [idx:4, marker:"mRAGE", role:"membrane", expectedCompartment:"alveolar", qcColor:"white", fileLabel:"mRAGE-647",
                areaMarker:true, areaMode:"membrane", areaMinAreaUm2:2.0d, areaBlurSigmaPx:0.7d] ],
    classify:[ ["tdTOM":true], ["T1A":true,"tdTOM":true],
               ["mRAGE":true,"tdTOM":true] ] ],

  // ---------------------------------------------------------------------
  // PILOT / PLUMBING TEST ONLY -- NOT a study panel. Delete when done.
  // Maps the OME sample file ND2/karl/sample_image.nd2 (Nikon CSU, 20x,
  // 5 channels: 1=far-red 2=red 3=green 4=blue/DAPI 5=brightfield) onto the
  // panel-A *shape* (nuclear + cyto/areaMarker + membrane).
  // The DAPI-equivalent counterstain is channel 4, NOT channel 1 -- this is
  // the only panel here whose nuclear idx is not 1.
  // The marker names below are structural placeholders: the green/red
  // channels are smFISH RNA probes from a skin sample, so KRT5/AGER
  // positivity numbers this produces are MEANINGLESS. See ref_images/README.md.
  "T": [ label:"T_PILOT_ome_nd2",
    channels:[ [idx:4, marker:"DAPI", role:"nuclear"],
               [idx:3, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:2, marker:"AGER", role:"membrane"] ],
    classify:[ ["KRT5":true,"AGER":false], ["KRT5":true,"AGER":true] ] ],

  // FUTURE: 4 channels exceeds the 3-marker slide limit; use single plane for YAP
  // (set PROJECTION="single") so the nuclear:cytoplasmic ratio is not MIP-corrupted.
  "S2": [ label:"S2_KRT5_p63_YAP",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",       areaMarker:true],
               [idx:3, marker:"p63",  role:"nuc_marker"],
               [idx:4, marker:"YAP",  role:"nuc_ratio"] ],
    classify:[ ["KRT5":true,"p63":true], ["KRT5":true,"YAP":true] ] ]
]

// ---------------------------------------------------------------------------
// Optional universal marker registry + study-specific panel configuration
// ---------------------------------------------------------------------------
// The registry never diagnoses a disease and never assigns a cell identity by
// itself. It records aliases, localization, modalities, and research contexts,
// and can supply a default analytical role when a custom panel omits one.
def normalizeMarkerToken = { value ->
  value == null ? "" : value.toString().toUpperCase().replaceAll(/[^A-Z0-9]+/, "")
}
def markerRegistryFile = new File(MARKER_REGISTRY_PATH)
def MARKER_REGISTRY = [schema_version:"unavailable", markers:[:], research_profiles:[:]]
if (markerRegistryFile.isFile()) {
  try {
    def parsed = new JsonSlurper().parse(markerRegistryFile)
    if (!(parsed instanceof Map) || !(parsed.markers instanceof Map)) {
      throw new IllegalArgumentException("registry root must contain a 'markers' object")
    }
    MARKER_REGISTRY = parsed
  } catch (Throwable t) {
    throw new IllegalArgumentException("Cannot parse IFQ_MARKER_REGISTRY '" + markerRegistryFile + "': " + t.message, t)
  }
} else {
  IJ.log("Marker registry not found; explicit channel roles remain fully supported: " + markerRegistryFile)
}

def markerProfileIndex = [:]
(MARKER_REGISTRY.markers ?: [:]).each { canonical, profile ->
  ([canonical] + (profile.aliases ?: [])).each { alias ->
    String token = normalizeMarkerToken(alias)
    if (!token.isEmpty()) markerProfileIndex[token] = [canonical:canonical, profile:profile]
  }
}
def markerProfileFor = { marker -> markerProfileIndex[normalizeMarkerToken(marker)] }

// A custom panel file is opt-in. It can add panels but cannot silently replace
// the validated built-ins. See config/custom_panels.example.json.
def CUSTOM_PANEL_KEYS = []
if (PANEL_CONFIG_PATH != null && !PANEL_CONFIG_PATH.trim().isEmpty()) {
  def panelConfigFile = new File(PANEL_CONFIG_PATH)
  if (!panelConfigFile.isFile()) {
    throw new IllegalArgumentException("IFQ_PANEL_CONFIG is not a file: " + panelConfigFile)
  }
  def customDoc
  try {
    customDoc = new JsonSlurper().parse(panelConfigFile)
  } catch (Throwable t) {
    throw new IllegalArgumentException("Cannot parse IFQ_PANEL_CONFIG '" + panelConfigFile + "': " + t.message, t)
  }
  if (!(customDoc instanceof Map) || !(customDoc.panels instanceof Map) || customDoc.panels.isEmpty()) {
    throw new IllegalArgumentException("IFQ_PANEL_CONFIG must contain a non-empty 'panels' object")
  }
  customDoc.panels.each { rawKey, rawPanel ->
    String panelKey = rawKey.toString()
    if (!(panelKey ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
      throw new IllegalArgumentException("Invalid custom panel key: " + panelKey)
    }
    if (PANELS.containsKey(panelKey)) {
      throw new IllegalArgumentException("Custom panel key would replace a built-in panel: " + panelKey)
    }
    if (!(rawPanel instanceof Map)) {
      throw new IllegalArgumentException("Custom panel '" + panelKey + "' must be an object")
    }
    if (!(rawPanel.channels instanceof Collection) || rawPanel.channels.isEmpty()) {
      throw new IllegalArgumentException("Custom panel '" + panelKey + "' needs a non-empty channels array")
    }
    def channels = rawPanel.channels.collect { rawChannel ->
      if (!(rawChannel instanceof Map)) {
        throw new IllegalArgumentException("Every channel in panel '" + panelKey + "' must be an object")
      }
      def c = [:]
      c.putAll(rawChannel)
      def matchedProfile = markerProfileFor(c.marker)
      if ((c.role == null || c.role.toString().trim().isEmpty()) && matchedProfile != null) {
        c.role = matchedProfile.profile.default_role
      }
      if (c.measurement == null && matchedProfile?.profile?.default_measurement != null) {
        c.measurement = matchedProfile.profile.default_measurement
      }
      if (matchedProfile != null) c.registryKey = matchedProfile.canonical
      if (c.role == "regional_area") {
        c.cellCall = false
        c.areaMarker = true
        if (c.areaMode == null) c.areaMode = "generic"
      }
      return c
    }
    PANELS[panelKey] = [label:(rawPanel.label ?: panelKey).toString(),
                        channels:channels, classify:(rawPanel.classify ?: [])]
    CUSTOM_PANEL_KEYS << panelKey
  }
}

def ALLOWED_CHANNEL_ROLES = ["nuclear", "cyto", "membrane", "nuc_marker",
                             "nuc_ratio", "apical_cilia", "regional_area"] as Set
PANELS.each { panelKey, panelDef ->
  if (!(panelDef.channels instanceof Collection) || panelDef.channels.isEmpty()) {
    throw new IllegalArgumentException("Panel '" + panelKey + "' has no channels")
  }
  panelDef.channels.each { c ->
    if (c.idx == null || (c.idx as int) < 1) {
      throw new IllegalArgumentException("Panel '" + panelKey + "' has an invalid channel idx")
    }
    if (c.marker == null || c.marker.toString().trim().isEmpty()) {
      throw new IllegalArgumentException("Panel '" + panelKey + "' has a channel without a marker")
    }
    if (c.role == null || !ALLOWED_CHANNEL_ROLES.contains(c.role.toString())) {
      throw new IllegalArgumentException("Panel '" + panelKey + "', marker '" + c.marker +
                                         "' needs one of roles " + ALLOWED_CHANNEL_ROLES)
    }
    if (c.role == "regional_area") {
      c.cellCall = false
      c.areaMarker = true
      if (c.areaMode == null) c.areaMode = "generic"
      if (c.areaMinAreaUm2 == null) {
        throw new IllegalArgumentException("Area-only marker '" + c.marker + "' in panel '" + panelKey +
                                           "' needs areaMinAreaUm2; no universal biological size cutoff exists")
      }
    }
    ["minPositiveFraction", "minLargestComponentShare"].each { field ->
      if (c[field] != null && (((double)c[field]) < 0.0d || ((double)c[field]) > 1.0d)) {
        throw new IllegalArgumentException(field + " must be between 0 and 1 for marker '" + c.marker + "'")
      }
    }
  }
  def indexes = panelDef.channels.collect { it.idx as int }
  if (indexes.unique().size() != indexes.size()) {
    throw new IllegalArgumentException("Panel '" + panelKey + "' repeats a channel idx")
  }
  def markerNames = panelDef.channels.collect { it.marker.toString() }
  if (markerNames.unique().size() != markerNames.size()) {
    throw new IllegalArgumentException("Panel '" + panelKey + "' repeats a marker name")
  }
  if (panelDef.channels.count { it.role == "nuclear" } != 1) {
    throw new IllegalArgumentException("Panel '" + panelKey + "' must contain exactly one nuclear segmentation channel")
  }
  def callableMarkers = panelDef.channels.findAll { it.role != "nuclear" && it.cellCall != false }
                                        .collect { it.marker.toString() } as Set
  (panelDef.classify ?: []).each { classRule ->
    if (!(classRule instanceof Map) || classRule.isEmpty()) {
      throw new IllegalArgumentException("Panel '" + panelKey + "' contains an empty classification rule")
    }
    classRule.each { marker, wanted ->
      if (!callableMarkers.contains(marker.toString())) {
        throw new IllegalArgumentException("Panel '" + panelKey + "' classifies unavailable/area-only marker '" + marker + "'")
      }
      if (!(wanted instanceof Boolean)) {
        throw new IllegalArgumentException("Classification values must be true/false in panel '" + panelKey + "'")
      }
    }
  }
}
if (!PANELS.containsKey(PANEL)) {
  throw new IllegalArgumentException("IFQ_PANEL '" + PANEL + "' is unknown. Available panels: " + PANELS.keySet().sort())
}

// Extend thresholds and morphology rules to every marker that appears in any
// panel. Explicit marker rules remain intact; new markers inherit only their
// geometry-class defaults until a study validates channel-specific settings.
def allAnalysisChannels = PANELS.values().collectMany { it.channels }
                               .findAll { it.role != "nuclear" }
allAnalysisChannels.each { c ->
  String marker = c.marker.toString()
  String token = normalizeMarkerToken(marker)
  def rawThreshold = System.getenv("IFQ_" + token + "_THRESHOLD")
  if (rawThreshold != null && !rawThreshold.trim().isEmpty()) {
    FIXED_POS_THRESHOLDS[marker] = Double.parseDouble(rawThreshold.trim())
  }
  if (!POS_SENSITIVITY.containsKey(marker)) POS_SENSITIVITY[marker] = 1.0d
  if (c.cellCall != false && !MORPHOLOGY_RULES.containsKey(marker)) {
    def roleRule = ROLE_MORPHOLOGY_DEFAULTS[c.role]
    if (roleRule == null) {
      throw new IllegalArgumentException("No cell-call morphology defaults for role '" + c.role + "'")
    }
    MORPHOLOGY_RULES[marker] = new LinkedHashMap(roleRule)
  }
}
// Apply environment overrides after custom markers have been registered.
MORPHOLOGY_RULES.each { marker, rule ->
  String token = normalizeMarkerToken(marker)
  rule.minFraction = Double.parseDouble(envOr("IFQ_" + token + "_MIN_POSITIVE_FRACTION", rule.minFraction.toString()))
  rule.minLargestShare = Double.parseDouble(envOr("IFQ_" + token + "_MIN_LARGEST_COMPONENT_SHARE", rule.minLargestShare.toString()))
  if (rule.containsKey("minNuclearEnrichment")) {
    rule.minNuclearEnrichment = Double.parseDouble(envOr("IFQ_" + token + "_MIN_NUCLEAR_ENRICHMENT", rule.minNuclearEnrichment.toString()))
  }
  if (rule.containsKey("minNucCytoRatio")) {
    rule.minNucCytoRatio = Double.parseDouble(envOr("IFQ_" + token + "_MIN_NUC_CYTO_RATIO", rule.minNucCytoRatio.toString()))
  }
}

// ============================================================================
//  3. PROVENANCE
// ============================================================================

def captureVersions() {
  def v = [:]
  v.imagej_version = IJ.getFullVersion()
  try { v.bioformats_version = FormatTools.VERSION } catch (e) { v.bioformats_version = "unknown" }
  v.java_version = System.getProperty("java.version")
  v.os = System.getProperty("os.name") + " " + System.getProperty("os.arch")
  // StarDist / CSBDeep versions are not reliably queryable; record from the
  // Updater if you need exact pins. Model + params below fully define the run.
  v.stardist_note = "record CSBDeep+StarDist versions from Help>Update if needed"
  v.timestamp = new Date().format("yyyy-MM-dd'T'HH:mm:ss")
  return v
}

// ============================================================================
//  4. HELPERS
// ============================================================================

def ensureDir(String p) { def d = new File(p); if (!d.exists()) d.mkdirs(); return d }

// Headless-safe replacement for Analyze Particles + ROI Manager. When an ROI
// is present, clear pixels outside it first: ParticleAnalyzer only honors the
// ROI bounds reliably and can otherwise count particles from excluded parts of
// a non-rectangular/composite ROI.
def particlesToRois(ImagePlus imp, double minAreaCal, boolean excludeEdges) {
  ImagePlus work = imp
  Roi restriction = imp.getRoi()
  if (restriction != null) {
    // ImagePlus.duplicate() may crop to the active ROI. Duplicate the processor
    // directly so the restriction remains in the original image coordinates.
    work = new ImagePlus(imp.getTitle() + "_roi_restricted",
                         imp.getProcessor().duplicate())
    work.setCalibration(imp.getCalibration())
    ImageProcessor wp = work.getProcessor()
    Rectangle rb = restriction.getBounds()
    ImageProcessor rm = restriction.getMask()
    int rx = (int)rb.x, ry = (int)rb.y
    int rw = (int)rb.width, rh = (int)rb.height
    // Clear outside explicitly. ImageProcessor.fillOutside(ShapeRoi) changes
    // processor ROI/mask state and is unreliable with these composite ROIs.
    for (int y = 0; y < wp.getHeight(); y++) {
      for (int x = 0; x < wp.getWidth(); x++) {
        boolean inBounds = x >= rx && x < rx + rw && y >= ry && y < ry + rh
        if (!inBounds || (rm != null && rm.get(x - rx, y - ry) == 0)) wp.set(x, y, 0)
      }
    }
  }

  int opts = ParticleAnalyzer.SHOW_ROI_MASKS
  if (excludeEdges) opts |= ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES
  def rt = new ResultsTable()
  // ParticleAnalyzer's Java constructor expects pixel counts (the interactive
  // dialog performs this conversion before constructing it). Convert the
  // public calibrated µm² setting explicitly.
  def workCal = work.getCalibration()
  double pixelArea = workCal.pixelWidth * workCal.pixelHeight
  double minAreaPixels = pixelArea > 0 ? minAreaCal / pixelArea : minAreaCal
  def pa = new ParticleAnalyzer(opts, Measurements.AREA, rt,
                                minAreaPixels, Double.MAX_VALUE)
  pa.setHideOutputImage(true)
  ImageProcessor src = work.getProcessor()
  // Every caller supplies a 0/255 binary mask. Ignore any threshold state left
  // behind by Convert to Mask/Watershed and always select non-zero foreground.
  src.setThreshold(128, 255, ImageProcessor.NO_LUT_UPDATE)
  if (!pa.analyze(work)) {
    if (!work.is(imp)) work.close()
    return []
  }
  ImagePlus labels = pa.getOutputImage()
  if (!work.is(imp)) work.close()
  if (labels == null) return []

  // SHOW_ROI_MASKS assigns consecutive integer labels. The former conversion
  // thresholded the complete 2048x2048 label image once per object, making ROI
  // extraction O(objects * image pixels). Find all label bounds in one pass,
  // then trace only each small cropped mask. This preserves the exact particle
  // shapes and coordinates while reducing a multi-minute hotspot to seconds.
  ImageProcessor lp = labels.getProcessor()
  int nLabels = rt.size()
  int[] minX = new int[nLabels + 1]
  int[] minY = new int[nLabels + 1]
  int[] maxX = new int[nLabels + 1]
  int[] maxY = new int[nLabels + 1]
  java.util.Arrays.fill(minX, lp.getWidth())
  java.util.Arrays.fill(minY, lp.getHeight())
  java.util.Arrays.fill(maxX, -1)
  java.util.Arrays.fill(maxY, -1)
  for (int y = 0; y < lp.getHeight(); y++) {
    for (int x = 0; x < lp.getWidth(); x++) {
      int label = lp.get(x, y)
      if (label < 1 || label > nLabels) continue
      if (x < minX[label]) minX[label] = x
      if (x > maxX[label]) maxX[label] = x
      if (y < minY[label]) minY[label] = y
      if (y > maxY[label]) maxY[label] = y
    }
  }
  def out = []
  for (int i = 1; i <= nLabels; i++) {
    if (maxX[i] < minX[i] || maxY[i] < minY[i]) continue
    int w = maxX[i] - minX[i] + 1
    int h = maxY[i] - minY[i] + 1
    def bp = new ByteProcessor(w, h)
    for (int yy = 0; yy < h; yy++) {
      for (int xx = 0; xx < w; xx++) {
        if (lp.get(minX[i] + xx, minY[i] + yy) == i) bp.set(xx, yy, 255)
      }
    }
    bp.setThreshold(128, 255, ImageProcessor.NO_LUT_UPDATE)
    def particle = new ImagePlus("particle_label_" + i, bp)
    def r = ThresholdToSelection.run(particle)
    if (r != null) {
      def rb = r.getBounds()
      r.setLocation(minX[i] + rb.x, minY[i] + rb.y)
      out << r
    }
    particle.close()
  }
  labels.close()
  return out
}

// Headless-safe replacement for RoiManager.runCommand("Open", ...).
def readRoiFile(File f) {
  def rois = []
  if (f.getName().toLowerCase().endsWith(".roi")) {
    def r = new RoiDecoder(f.getAbsolutePath()).getRoi()
    if (r != null) {
      if (r.getName() == null) r.setName(f.getName() - ".roi")
      rois << r
    }
    return rois
  }
  def zis = new ZipInputStream(new BufferedInputStream(new FileInputStream(f)))
  try {
    def entry
    while ((entry = zis.getNextEntry()) != null) {
      if (!entry.getName().toLowerCase().endsWith(".roi")) continue
      def baos = new ByteArrayOutputStream()
      byte[] buf = new byte[8192]
      int len
      while ((len = zis.read(buf)) > 0) baos.write(buf, 0, len)
      def r = new RoiDecoder(baos.toByteArray(), entry.getName()).getRoi()
      if (r != null) {
        r.setName(entry.getName() - ".roi")
        rois << r
      }
    }
  } finally {
    zis.close()
  }
  return rois
}

// Import via Bio-Formats keeping metadata + calibration; grayscale, no split yet.
def bfOpen(String path) {
  def opts = new ImporterOptions()
  opts.setId(path)
  opts.setSplitChannels(false)
  opts.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
  opts.setVirtual(false)
  opts.setAutoscale(false)
  def imps = BF.openImagePlus(opts)
  return imps[0]
}

// Project one single-channel stack to 2D per PROJECTION setting.
def projectChannel(ImagePlus ch, String projection, int singlePlane) {
  if (ch.getNSlices() <= 1) return ch
  if (projection == "single") {
    int z = (singlePlane >= 1) ? singlePlane : (int)Math.ceil(ch.getNSlices()/2.0)
    ch.setSlice(z)
    def ip = ch.getProcessor().duplicate()
    def out = new ImagePlus(ch.getTitle(), ip)
    out.setCalibration(ch.getCalibration())
    return out
  }
  def method = (projection == "sum") ? "sum" : (projection == "avg") ? "avg" : "max"
  def out = ZProjector.run(ch, method)
  out.setCalibration(ch.getCalibration())
  return out
}

// Iterate ROI pixels safely (respects non-rectangular mask).
def eachRoiPixel(Roi roi, Closure body) {
  Rectangle b = roi.getBounds()
  ImageProcessor mask = roi.getMask()
  int bx = (int)b.x, by = (int)b.y
  int bw = (int)b.width, bh = (int)b.height
  for (int y = 0; y < bh; y++) {
    for (int x = 0; x < bw; x++) {
      if (mask == null || mask.get(x, y) != 0) body((int)(bx + x), (int)(by + y))
    }
  }
}

// Bit-depth-robust Otsu (or named method) on raw intensities within an ROI.
def autoThresholdInRoi(ImagePlus imp, Roi roi, String method) {
  ImageProcessor ip = imp.getProcessor()
  ip.setRoi(roi)
  ImageStatistics st = ImageStatistics.getStatistics(ip, Measurements.MIN_MAX, imp.getCalibration())
  double mn = st.min, mx = st.max
  if (mx <= mn) return mx
  int[] hist = new int[256]
  eachRoiPixel(roi) { x, y ->
    double v = ip.getPixelValue(x, y)
    int bin = (int)Math.round((v - mn) / (mx - mn) * 255.0)
    if (bin < 0) bin = 0; if (bin > 255) bin = 255
    hist[bin]++
  }
  def m = AutoThresholder.Method.valueOf(method)
  int tbin = new AutoThresholder().getThreshold(m, hist)
  return mn + (tbin / 255.0) * (mx - mn)
}

// Mean raw intensity + area of one ROI on one channel.
def measureRoi(ImagePlus imp, Roi roi) {
  ImageProcessor ip = imp.getProcessor()
  ip.setRoi(roi)
  ImageStatistics st = ImageStatistics.getStatistics(
      ip, Measurements.MEAN | Measurements.AREA | Measurements.CENTROID,
      imp.getCalibration())
  return [mean: st.mean, area: st.area, cx: st.xCentroid, cy: st.yCentroid]
}

// Count positive (>0) pixels of a mask inside an ROI -> calibrated area.
def positiveAreaInRoi(ImagePlus maskImp, Roi roi) {
  ImageProcessor mp = maskImp.getProcessor()
  Calibration cal = maskImp.getCalibration()
  Rectangle b = roi.getBounds()
  ImageProcessor roiMask = roi.getMask()
  int bx = (int)b.x, by = (int)b.y
  int x0 = Math.max(0, bx), y0 = Math.max(0, by)
  int x1 = Math.min(mp.getWidth(), bx + (int)b.width)
  int y1 = Math.min(mp.getHeight(), by + (int)b.height)
  long pos = 0L
  for (int y = y0; y < y1; y++) {
    for (int x = x0; x < x1; x++) {
      if ((roiMask == null || roiMask.get(x - bx, y - by) != 0) &&
          mp.get(x, y) != 0) pos++
    }
  }
  return pos * cal.pixelWidth * cal.pixelHeight
}

// Fraction of pixels in an object/ring whose raw intensity reaches threshold.
// This separates a spatially supported membrane pattern from a bright speck
// that happens to raise the object's mean.
def fractionAboveThreshold(ImagePlus imp, Roi roi, double threshold) {
  ImageProcessor ip = imp.getProcessor()
  long total = 0L, positive = 0L
  eachRoiPixel(roi) { x, y ->
    if (x >= 0 && y >= 0 && x < ip.getWidth() && y < ip.getHeight()) {
      total++
      if (ip.getPixelValue(x, y) >= threshold) positive++
    }
  }
  return total > 0 ? positive / (double)total : 0.0d
}

// Morphology support inside an object ROI. In addition to total positive
// coverage, report how much of the positive signal belongs to the largest
// 8-connected component. A high largest-component share distinguishes a
// coherent nuclear/cytoplasmic/membrane pattern from scattered bright specks.
def spatialSupportStats(ImagePlus imp, Roi roi, double threshold) {
  if (roi == null) return [total:0L, positive:0L, fraction:0.0d,
                           components:0, largest:0L, largestShare:0.0d]
  ImageProcessor ip = imp.getProcessor()
  Rectangle b = roi.getBounds()
  ImageProcessor rm = roi.getMask()
  int bw = (int)b.width, bh = (int)b.height
  if (bw <= 0 || bh <= 0) return [total:0L, positive:0L, fraction:0.0d,
                                  components:0, largest:0L, largestShare:0.0d]
  boolean[] positiveMask = new boolean[bw * bh]
  long total = 0L, positive = 0L
  for (int yy = 0; yy < bh; yy++) {
    for (int xx = 0; xx < bw; xx++) {
      if (rm != null && rm.get(xx, yy) == 0) continue
      int x = (int)b.x + xx, y = (int)b.y + yy
      if (x < 0 || y < 0 || x >= ip.getWidth() || y >= ip.getHeight()) continue
      total++
      if (ip.getPixelValue(x, y) >= threshold) {
        positiveMask[yy * bw + xx] = true
        positive++
      }
    }
  }

  boolean[] visited = new boolean[bw * bh]
  int[] queue = new int[bw * bh]
  int components = 0
  long largest = 0L
  int[] dx = [-1, 0, 1, -1, 1, -1, 0, 1] as int[]
  int[] dy = [-1, -1, -1, 0, 0, 1, 1, 1] as int[]
  for (int start = 0; start < positiveMask.length; start++) {
    if (!positiveMask[start] || visited[start]) continue
    components++
    int head = 0, tail = 0
    queue[tail++] = start
    visited[start] = true
    long size = 0L
    while (head < tail) {
      int idx = queue[head++]
      size++
      int x0 = idx % bw, y0 = (int)(idx / bw)
      for (int k = 0; k < 8; k++) {
        int nx = x0 + dx[k], ny = y0 + dy[k]
        if (nx < 0 || ny < 0 || nx >= bw || ny >= bh) continue
        int ni = ny * bw + nx
        if (positiveMask[ni] && !visited[ni]) {
          visited[ni] = true
          queue[tail++] = ni
        }
      }
    }
    if (size > largest) largest = size
  }
  return [total:total, positive:positive,
          fraction:(total > 0 ? positive / (double)total : 0.0d),
          components:components, largest:largest,
          largestShare:(positive > 0 ? largest / (double)positive : 0.0d)]
}

// Strict ownership screen for nucleus-associated measurements. If a support
// territory encloses another nucleus centroid, pixels cannot be assigned to one
// cell unambiguously without a full membrane/cell segmentation; leave the final
// call indeterminate instead of double-counting shared signal.
def supportHasOtherNucleus(Roi support, int currentIndex, nuclei) {
  if (support == null) return true
  for (int j = 0; j < nuclei.size(); j++) {
    if (j == currentIndex) continue
    Rectangle nb = nuclei[j].getBounds()
    int cx = (int)Math.round(nb.getCenterX())
    int cy = (int)Math.round(nb.getCenterY())
    if (support.contains(cx, cy)) return true
  }
  return false
}

def buildMaskAtThreshold(ImagePlus ch, double threshold) {
  ImagePlus dup = ch.duplicate()
  double upper = ch.getBitDepth() == 8 ? 255.0d : (ch.getBitDepth() == 16 ? 65535.0d : Double.MAX_VALUE)
  dup.getProcessor().setThreshold(threshold, upper, ImageProcessor.NO_LUT_UPDATE)
  IJ.run(dup, "Convert to Mask", "")
  dup.setCalibration(ch.getCalibration())
  dup.setProperty("thresholdValue", threshold)
  return dup
}

// Build a binary mask ImagePlus (255 = signal) from a channel by threshold.
// The numeric lower threshold is stashed as a property for provenance export.
// Requires Prefs.blackBackground=true (set in main) so foreground = 255.
def buildThresholdMask(ImagePlus ch, double blurSigma, String method, Double fixedThreshold = null) {
  ImagePlus dup = ch.duplicate()
  if (blurSigma > 0) new GaussianBlur().blurGaussian(dup.getProcessor(), blurSigma)
  double upper = dup.getBitDepth() == 8 ? 255.0d : (dup.getBitDepth() == 16 ? 65535.0d : Double.MAX_VALUE)
  double thr
  if (fixedThreshold != null) {
    thr = fixedThreshold
    dup.getProcessor().setThreshold(thr, upper, ImageProcessor.NO_LUT_UPDATE)
  } else {
    IJ.setAutoThreshold(dup, method + " dark")
    thr = dup.getProcessor().getMinThreshold()   // raw intensity, -1 if unset
  }
  IJ.run(dup, "Convert to Mask", "")
  dup.setCalibration(ch.getCalibration())
  dup.setProperty("thresholdValue", thr)
  return dup
}

// Remove connected foreground components below the declared physical area.
// Applying this before both area measurement and mask export keeps the numeric
// endpoint, component table, and saved QC mask internally consistent.
def filterBinaryMaskByArea(ImagePlus mask, double minAreaUm2) {
  def work = mask.duplicate()
  work.setCalibration(mask.getCalibration())
  def accepted = particlesToRois(work, minAreaUm2, false)
  work.close()
  def bp = new ByteProcessor(mask.getWidth(), mask.getHeight())
  bp.setValue(255)
  accepted.each { bp.fill(it) }
  def out = new ImagePlus(mask.getTitle() + "_area_filtered", bp)
  out.setCalibration(mask.getCalibration())
  out.setProperty("thresholdValue", mask.getProperty("thresholdValue"))
  out.setProperty("minimumComponentAreaUm2", minAreaUm2)
  return out
}

// Perinuclear cytoplasm only = enlarged cell ROI minus the nucleus ROI.
// Used for a true nuclear:cytoplasmic ratio (e.g. YAP). null if degenerate.
def ringOnly(Roi nuc, Roi cell) {
  try {
    ShapeRoi s = new ShapeRoi(cell).not(new ShapeRoi(nuc))
    Rectangle b = s.getBounds()
    if (b == null || b.width <= 0 || b.height <= 0) return null
    return s
  } catch (Throwable e) { return null }
}

// Tissue ROI: manual if RoiSet/.roi present, else auto from DAPI.
def resolveTissueRois(String imgPath, ImagePlus dapi, cfg) {
  def base = new File(imgPath)
  def stem = base.name.replaceFirst(/\.[^.]+$/, "")
  def parent = base.getParent()
  def candidates = [ new File(parent, stem + "_RoiSet.zip"),
                     new File(parent, stem + ".zip"),
                     new File(parent, stem + ".roi"),
                     new File(parent, "RoiSet.zip") ]
  def hit = candidates.find { it.exists() }
  if (hit != null) {
    def rois = readRoiFile(hit)
    def named = []
    rois.eachWithIndex { r, i -> named << [name: (r.getName() ?: ("region" + (i+1))), roi: r] }
    return [source: hit.name, regions: named]
  }
  if (cfg.tissueMode == "whole_field") {
    return [source: "whole_field", regions: [[name: "whole_field",
             roi: new Roi(0, 0, dapi.getWidth(), dapi.getHeight())]]]
  }
  // Auto tissue from DAPI
  def mask = buildThresholdMask(dapi, cfg.tissueBlur, cfg.tissueMethod)
  IJ.run(mask, "Options...", "iterations=2 count=1 do=Close")
  def rois = particlesToRois(mask, cfg.tissueMinArea, false)
  mask.close()
  if (rois.isEmpty()) {  // fall back to whole field
    return [source: "whole_field", regions: [[name: "tissue", roi: new Roi(0, 0, dapi.getWidth(), dapi.getHeight())]]]
  }
  // merge all particles into one "tissue" region
  ShapeRoi merged = null
  rois.each { r -> def s = new ShapeRoi(r); merged = (merged == null) ? s : merged.or(s) }
  return [source: "auto_dapi", regions: [[name: "tissue", roi: merged]]]
}

// Nucleus ROIs within a given tissue region.
def segmentNuclei(ImagePlus dapi, Roi region, cfg) {
  ImagePlus crop = dapi.duplicate()
  crop.setRoi(region)
  if (cfg.segmenter == "stardist") {
    // StarDist's ROI Manager output is interactive-only. Keep it isolated from
    // the classic path so classic segmentation remains fully headless-safe.
    def rm = ij.plugin.frame.RoiManager.getInstance() ?: new ij.plugin.frame.RoiManager()
    rm.reset()
    crop.setTitle("DAPI_seg")
    crop.show()   // StarDist is happiest with a shown image (interactive mode)
    IJ.run(crop, "Command From Macro",
      "command=[de.csbdresden.stardist.StarDist2D], " +
      "args=['input':'DAPI_seg', 'modelChoice':'Versatile (fluorescent nuclei)', " +
      "'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', " +
      "'probThresh':'" + cfg.prob + "', 'nmsThresh':'" + cfg.nms + "', " +
      "'outputType':'ROI Manager', 'nTiles':'" + cfg.tiles + "', " +
      "'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', " +
      "'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]")
    def rois = rm.getRoisAsArray().collect { it }
    crop.changes = false; crop.close()
    // keep only nuclei whose centroid lies inside the region
    def included = rois.findAll { region.contains((int)it.getBounds().getCenterX(), (int)it.getBounds().getCenterY()) }
    return [included: included, rejected: []]
  } else {
    // Classic/local-threshold watershed fallback. The local mode is intended
    // for uneven DAPI illumination and is fully recorded in params.json.
    def m = crop.duplicate()
    if (cfg.dapiMethod == "local_phansalkar") {
      double px = Math.max(dapi.getCalibration().pixelWidth, 1.0e-9d)
      int backgroundRadiusPx = Math.max(3, (int)Math.round(cfg.dapiBackgroundRadiusUm / px))
      int localRadiusPx = Math.max(3, (int)Math.round(cfg.dapiLocalRadiusUm / px))
      IJ.run(m, "Subtract Background...", "rolling=" + backgroundRadiusPx + " sliding")
      if (cfg.dapiBlurSigmaPx > 0) new GaussianBlur().blurGaussian(m.getProcessor(), cfg.dapiBlurSigmaPx)
      IJ.run(m, "Enhance Contrast...", "saturated=" + cfg.dapiContrastSaturation + " normalize")
      if (m.getBitDepth() != 8) IJ.run(m, "8-bit", "")
      IJ.run(m, "Auto Local Threshold",
             "method=Phansalkar radius=" + localRadiusPx + " parameter_1=0 parameter_2=0 white")
    } else {
      new GaussianBlur().blurGaussian(m.getProcessor(), 2.0)
      IJ.setAutoThreshold(m, "Otsu dark")
      IJ.run(m, "Convert to Mask", "")
    }
    IJ.run(m, "Fill Holes", "")
    IJ.run(m, "Watershed", "")
    m.setCalibration(dapi.getCalibration())
    m.setRoi(region)
    // Keep rejected candidates for QC instead of silently dropping them.
    // A candidate is the same Otsu/watershed particle used for counting; the
    // only post-segmentation rejection rules are calibrated area and image edge.
    def candidates = particlesToRois(m, 0.0d, false)
    // Preserve ParticleAnalyzer's calibrated size/edge decision as the
    // authoritative count. Diagnostics are derived afterward and cannot alter it.
    def rois = particlesToRois(m, cfg.minNucArea, true)
    def roiKey = { r ->
      def b = r.getBounds()
      return b.x + ":" + b.y + ":" + b.width + ":" + b.height
    }
    def includedKeys = rois.collect { roiKey(it) } as Set
    def rejected = candidates.findAll { !includedKeys.contains(roiKey(it)) }.collect { r ->
      def b = r.getBounds()
      double area = measureRoi(dapi, r).area
      boolean edge = b.x <= 0 || b.y <= 0 || b.x + b.width >= dapi.getWidth() || b.y + b.height >= dapi.getHeight()
      return [roi: r,
              reason: (edge ? "image_edge" : (area < cfg.minNucArea ? "area_below_minimum" : "particle_filter")),
              area_um2: area]
    }
    def candidateMask = m.duplicate(); candidateMask.setCalibration(dapi.getCalibration())
    m.close()
    crop.close()
    return [included: rois, rejected: rejected, candidateMask: candidateMask]
  }
}

// ============================================================================
//  5. PER-IMAGE PROCESSING
// ============================================================================

def processImage(String imgPath, String outputKey, panelKey, panelDef, meta, cfg, outDir) {
  IJ.log("---- " + new File(imgPath).name + "  [panel " + panelKey + "] ----")
  def sourceStem = new File(imgPath).name.replaceFirst(/\.[^.]+$/, "")
  def imgOut = ensureDir(outDir + "/" + outputKey)
  def fileSafe = { value ->
    def s = (value == null || value.toString().trim().isEmpty()) ? "NA" : value.toString().trim()
    return s.replaceAll(/[^A-Za-z0-9._-]+/, "-").replaceAll(/^-+|-+$/, "")
  }
  // Use the marker map as the per-image filename prefix. The containing
  // directory already identifies the specimen, so avoiding that duplication
  // keeps Windows paths short enough to open reliably.
  def channelSignature = panelDef.channels.sort { it.idx }.collect { c ->
    "C" + c.idx + "-" + fileSafe(c.fileLabel ?: c.marker)
  }.join("_")
  def fileKey = channelSignature

  def raw = bfOpen(imgPath)
  Calibration cal = raw.getCalibration()
  // Channel maps may intentionally skip an unused acquisition channel. The
  // highest referenced index, not the number of mapped channels, defines the
  // minimum acquisition size.
  def nChExpected = panelDef.channels.collect { it.idx as int }.max()
  def channels = ChannelSplitter.split(raw)
  if (channels.length < nChExpected) {
    IJ.log("  WARNING: found " + channels.length + " channels, panel expects " + nChExpected + " -> skipped")
    raw.close()
    return null
  }

  // Map marker -> projected 2D channel image
  def markerImg = [:]
  def nuclearMarker = null
  panelDef.channels.each { c ->
    def ch = channels[c.idx - 1]
    ch.setCalibration(cal)
    def proj = projectChannel(ch, cfg.projection, cfg.singlePlane)
    proj.setCalibration(cal)
    markerImg[c.marker] = proj
    if (c.role == "nuclear") nuclearMarker = c.marker
  }
  def dapi = markerImg[nuclearMarker]
  def nonNuclearChannels = panelDef.channels.findAll { it.role != "nuclear" }
  def cellChannels = nonNuclearChannels.findAll { it.cellCall != false }
  def expectedCompartmentsFor = { c ->
    def expected = []
    if (c.expectedCompartments instanceof Collection) {
      expected.addAll(c.expectedCompartments.collect { it.toString().toLowerCase() })
    } else if (c.expectedCompartment != null) {
      expected << c.expectedCompartment.toString().toLowerCase()
    }
    return expected.unique()
  }

  // Precompute per-marker positivity thresholds (Otsu within whole field first;
  // refined per tissue region below).
  def tissue = resolveTissueRois(imgPath, dapi, cfg)

  // Pre-build marker-specific area masks. The minimum component area is applied
  // to the saved mask itself so reported area and component counts agree.
  def areaMasks = [:]
  def areaThresholdSources = [:]
  panelDef.channels.findAll { it.areaMarker }.each { c ->
    double areaBlur = c.containsKey("areaBlurSigmaPx") ? (double)c.areaBlurSigmaPx : cfg.podBlur
    String areaMethod = c.containsKey("areaThresholdMethod") ? c.areaThresholdMethod : cfg.podMethod
    double minArea = c.containsKey("areaMinAreaUm2") ? (double)c.areaMinAreaUm2 : cfg.podMinArea
    boolean fixedAreaThreshold = cfg.fixedThresholds.containsKey(c.marker)
    Double resolvedAreaThreshold = fixedAreaThreshold ? (double)cfg.fixedThresholds[c.marker] : null
    def rawAreaMask = buildThresholdMask(markerImg[c.marker], areaBlur, areaMethod, resolvedAreaThreshold)
    areaMasks[c.marker] = filterBinaryMaskByArea(rawAreaMask, minArea)
    areaThresholdSources[c.marker] = fixedAreaThreshold ? "fixed_predeclared" : "adaptive_otsu_exploratory"
    rawAreaMask.close()
  }

  def cellRows = []      // per-object records (all regions)
  def summaryRows = []   // per-region summary
  def qcRegionOverlays = [:]

  tissue.regions.eachWithIndex { reg, ri ->
    def region = reg.roi
    def regName = reg.name
    def regionLower = regName.toLowerCase()
    def regionTags = []
    if (regionLower.contains("alveol")) regionTags << "alveolar"
    if (regionLower.contains("airway") || regionLower.contains("bronch")) regionTags << "airway"
    if (regionLower.contains("tumor") || regionLower.contains("tumour") || regionLower.contains("luad")) regionTags << "tumor"
    if (regionLower.contains("fibrot") || regionLower.contains("honeycomb") || regionLower.contains("uip")) regionTags << "fibrotic"
    if (regionLower.contains("strom") || regionLower.contains("mesench")) regionTags << "stromal"
    if (regionLower.contains("vascul") || regionLower.contains("vessel") || regionLower.contains("capillar")) regionTags << "vascular"
    if (regionLower.contains("immune") || regionLower.contains("inflamm") || regionLower.contains("lymph")) regionTags << "immune"
    if (regionLower.contains("ambig")) regionTags = ["ambiguous"]
    regionTags = regionTags.unique()
    if (regionTags.isEmpty() && cfg.wholeFieldCompartment != "unassigned") {
      regionTags << cfg.wholeFieldCompartment
    }
    def compartment = regionTags.contains("ambiguous") ? "ambiguous" :
                      (regionTags.contains("alveolar") ? "alveolar" :
                       (regionTags.contains("airway") ? "airway" :
                        (regionTags.isEmpty() ? "unassigned" : regionTags[0])))
    if (cfg.compartmentMode == "required" && regionTags.isEmpty()) {
      throw new IllegalArgumentException("Morphology classification required: use a recognizable anatomical/context ROI name; found '" + regName + "'")
    }
    def regStat = measureRoi(dapi, region)   // area of region
    double regionAreaUm2 = regStat.area

    // per-marker channel thresholds inside THIS region (adaptive)
    def chThresh = [:]
    def chThreshSource = [:]
    nonNuclearChannels.each { c ->
      boolean fixed = cfg.fixedThresholds.containsKey(c.marker)
      double t = fixed ? (double)cfg.fixedThresholds[c.marker] : autoThresholdInRoi(markerImg[c.marker], region, "Otsu")
      double sens = fixed ? 1.0d : (cfg.sensitivity[c.marker] ?: 1.0)
      chThresh[c.marker] = t * sens
      chThreshSource[c.marker] = fixed ? "fixed_predeclared" : "adaptive_otsu_exploratory"
    }

    // Region boundaries for each fluorescence channel; these replace the
    // former per-cell marker-positive circles in the QC overlay.
    def qcMasks = [:]
    nonNuclearChannels.each { c ->
      qcMasks[c.marker] = buildMaskAtThreshold(markerImg[c.marker], chThresh[c.marker])
    }

    // ---- Marker-positive area/components (pods, reporter fields, ciliary patches) ----
    def areaStats = [:]
    panelDef.channels.findAll { it.areaMarker }.each { c ->
      def mask = areaMasks[c.marker]
      String areaMode = c.areaMode ?: "pod"
      double minComponentArea = c.containsKey("areaMinAreaUm2") ? (double)c.areaMinAreaUm2 : cfg.podMinArea
      double positiveArea = positiveAreaInRoi(mask, region)
      // Particle analysis is intentionally mode-specific: large KRT5 pods,
      // reporter-positive cell/clusters, or subcellular ciliary patches.
      def maskReg = mask.duplicate(); maskReg.setCalibration(cal); maskReg.setRoi(region)
      def componentRois = particlesToRois(maskReg, minComponentArea, false)
      def componentAreas = componentRois.collect { measureRoi(mask, it).area }
      def areaThr = mask.getProperty("thresholdValue")
      areaStats[c.marker] = [ mode: areaMode, area_um2: positiveArea,
                              frac_of_region: (regionAreaUm2 > 0 ? positiveArea/regionAreaUm2 : 0),
                              n_components: componentRois.size(),
                              mean_component_area_um2: (componentAreas.isEmpty()? 0 : componentAreas.sum()/componentAreas.size()),
                              min_component_area_um2: minComponentArea,
                              threshold: (areaThr != null ? areaThr : -1),
                              threshold_source: areaThresholdSources[c.marker] ]
      maskReg.close()
    }

    // ---- Nuclei -> cells ----
    def segmentation = segmentNuclei(dapi, region, cfg)
    def nuclei = segmentation.included
    def rejectedNuclei = segmentation.rejected ?: []
    def posCount = [:].withDefault { 0 }
    def finalPosCount = [:].withDefault { 0 }
    def finalNegCount = [:].withDefault { 0 }
    def indeterminateCount = [:].withDefault { 0 }
    def classCount = [:].withDefault { 0 }
    def classEvaluableCount = [:].withDefault { 0 }
    def allNucRois = []
    def finalPositiveRois = [:].withDefault { [] }
    def indeterminateRois = [:].withDefault { [] }

    nuclei.eachWithIndex { nuc, ni ->
      allNucRois << nuc
      def cellRoi = RoiEnlarger.enlarge(nuc, (int)Math.round(cfg.ringExpandUm / cal.pixelWidth))
      def row = [ image: sourceStem, panel: panelKey, region: regName,
                  compartment: compartment, region_tags: regionTags.join("|"), cell_id: (ni + 1) ]
      row.mouse_id = meta.mouse_id; row.section_id = meta.section_id
      row.genotype = meta.genotype; row.condition = meta.condition
      def cs = measureRoi(dapi, nuc)
      row.centroid_x_um = cs.cx; row.centroid_y_um = cs.cy; row.nucleus_area_um2 = cs.area

      def calls = [:]       // morphology-authoritative three-state calls: 1, 0, or ""
      def rawCalls = [:]    // legacy mean-intensity calls retained for audit only
      cellChannels.each { c ->
        def m = c.marker
        def img = markerImg[m]
        def rule = new LinkedHashMap(cfg.roleMorphologyDefaults[c.role] ?:
                                     [minFraction:0.20d, minLargestShare:0.40d, requireOwnership:true])
        if (cfg.morphologyRules[m] != null) rule.putAll(cfg.morphologyRules[m])
        // Per-channel overrides are the final authority because the same
        // antigen can require different geometry in different assays.
        if (c.minPositiveFraction != null) rule.minFraction = c.minPositiveFraction as double
        if (c.minLargestComponentShare != null) rule.minLargestShare = c.minLargestComponentShare as double
        if (c.requireOwnership != null) rule.requireOwnership = c.requireOwnership as boolean
        if (c.minNuclearEnrichment != null) rule.minNuclearEnrichment = c.minNuclearEnrichment as double
        if (c.minNucCytoRatio != null) rule.minNucCytoRatio = c.minNucCytoRatio as double
        def spatialRoi = nuc
        def ownershipSupport = null
        double val, nucVal = 0.0d, cytoVal = 0.0d, enrichmentRatio = 0.0d
        boolean projectionValid = true
        if (c.role == "nuc_marker") {
          val = measureRoi(img, nuc).mean
          nucVal = val
          def ring = ringOnly(nuc, cellRoi)
          cytoVal = ring != null ? measureRoi(img, ring).mean : 0.0d
          enrichmentRatio = cytoVal > 0 ? nucVal / cytoVal : (nucVal > 0 ? Double.POSITIVE_INFINITY : 0.0d)
          row[m + "_nuc_mean"] = nucVal
          row[m + "_reference_ring_mean"] = cytoVal
          row[m + "_nuclear_enrichment_ratio"] = enrichmentRatio
        } else if (c.role == "nuc_ratio") {
          nucVal = measureRoi(img, nuc).mean
          // true cytoplasmic ring = (enlarged cell) MINUS (nucleus)
          def ring = ringOnly(nuc, cellRoi)
          cytoVal = (ring != null) ? measureRoi(img, ring).mean : measureRoi(img, cellRoi).mean
          spatialRoi = (ring != null) ? ring : cellRoi
          // Morphology support is measured in the nucleus; the cytoplasmic
          // ring supplies the localization ratio, not the positive pixels.
          spatialRoi = nuc
          val = nucVal
          enrichmentRatio = cytoVal > 0 ? nucVal / cytoVal : (nucVal > 0 ? Double.POSITIVE_INFINITY : 0.0d)
          row[m + "_nuc_mean"] = nucVal
          row[m + "_cyto_mean"] = cytoVal
          row[m + "_nuc_cyto_ratio"] = enrichmentRatio
          projectionValid = raw.getNSlices() <= 1 || cfg.projection == "single"
        } else if (c.role == "apical_cilia") {
          // Ciliary axonemes sit on the luminal/apical surface and commonly lie
          // beyond a 2-um cytoplasmic ring in tissue sections. Use a wider
          // support zone, and call proximity by the spatial fraction above the
          // image/region threshold rather than by a diluted support-zone mean.
          double supportExpandUm = c.supportExpandUm != null ? (double)c.supportExpandUm : cfg.actubSupportExpandUm
          def support = RoiEnlarger.enlarge(nuc, (int)Math.round(supportExpandUm / cal.pixelWidth))
          spatialRoi = support ?: cellRoi
          ownershipSupport = spatialRoi
          val = measureRoi(img, spatialRoi).mean
          row[m + "_support_expand_um"] = supportExpandUm
          row[m + "_measurement_model"] = c.measurement ?: "apical_cilia_proximity"
        } else { // cyto / membrane: measure the perinuclear ring, not nucleus + ring
          def ring = ringOnly(nuc, cellRoi)
          spatialRoi = (ring != null) ? ring : cellRoi
          ownershipSupport = cellRoi
          val = measureRoi(img, spatialRoi).mean
          row[m + "_measurement_model"] = c.measurement ?: "perinuclear_ring"
        }

        if (c.requiresSinglePlane == true) {
          projectionValid = raw.getNSlices() <= 1 || cfg.projection == "single"
        }
        if (!row.containsKey(m + "_measurement_model")) {
          row[m + "_measurement_model"] = c.measurement ?: c.role
        }

        def supportStats = spatialSupportStats(img, spatialRoi, chThresh[m])
        double minFraction = (double)rule.minFraction
        double minLargestShare = (double)rule.minLargestShare
        boolean fractionPass = supportStats.fraction >= minFraction
        boolean connectedPass = supportStats.largestShare >= minLargestShare
        boolean ownershipClear = !(rule.requireOwnership ?: false) ||
                                 !supportHasOtherNucleus(ownershipSupport, ni, nuclei)
        boolean enrichmentPass = true
        if (c.role == "nuc_marker") {
          enrichmentPass = enrichmentRatio >= (double)(rule.minNuclearEnrichment ?: 1.0d)
        } else if (c.role == "nuc_ratio") {
          enrichmentPass = enrichmentRatio >= (double)(rule.minNucCytoRatio ?: 1.0d)
        }

        def expectedCompartments = expectedCompartmentsFor(c)
        boolean compartmentAssigned = !regionTags.isEmpty() && !regionTags.contains("ambiguous")
        boolean compartmentPass = expectedCompartments.isEmpty() ||
                                  (compartmentAssigned && expectedCompartments.any { regionTags.contains(it) })
        boolean evaluable = supportStats.total > 0 && projectionValid
        def indeterminateReasons = []
        if (!expectedCompartments.isEmpty() && !compartmentAssigned) {
          evaluable = false
          indeterminateReasons << "compartment_unassigned"
        }
        if (!ownershipClear) {
          evaluable = false
          indeterminateReasons << "shared_perinuclear_support"
        }
        if (!projectionValid) {
          indeterminateReasons << (c.role == "nuc_ratio" ?
                                    "projection_invalid_for_nuclear_ratio" :
                                    "projection_invalid_for_marker")
        }
        if (supportStats.total <= 0) {
          indeterminateReasons << "empty_spatial_support"
        }

        row[m + "_mean"] = val
        boolean intensityPos = val >= chThresh[m]
        rawCalls[m] = intensityPos ? 1 : 0
        row[m + "_pos"] = intensityPos ? 1 : 0  // legacy/raw audit field
        row[m + "_threshold_source"] = chThreshSource[m]
        row[m + "_support_fraction_above_threshold"] = supportStats.fraction
        row[m + "_minimum_support_fraction"] = minFraction
        row[m + "_positive_component_count"] = supportStats.components
        row[m + "_largest_positive_component_share"] = supportStats.largestShare
        row[m + "_minimum_largest_component_share"] = minLargestShare
        row[m + "_fraction_pass"] = fractionPass ? 1 : 0
        row[m + "_connected_pattern_pass"] = connectedPass ? 1 : 0
        row[m + "_ownership_clear"] = ownershipClear ? 1 : 0
        row[m + "_projection_valid"] = projectionValid ? 1 : 0
        row[m + "_expected_compartment"] = expectedCompartments.isEmpty() ? "none" : expectedCompartments.join("|")
        row[m + "_compartment_pass"] = compartmentPass ? 1 : 0
        row[m + "_enrichment_pass"] = enrichmentPass ? 1 : 0
        if (c.role == "membrane" || c.role == "cyto") {
          row[m + "_ring_fraction_above_threshold"] = supportStats.fraction
          row[m + "_minimum_ring_fraction"] = minFraction
        }
        row[m + "_pattern_pos"] = (fractionPass && connectedPass) ? 1 : 0
        row[m + "_compartment_consistent"] = expectedCompartments.isEmpty() ? 1 :
                                              (compartmentAssigned ? (compartmentPass ? 1 : 0) : "")

        if (intensityPos) {
          posCount[m] = posCount[m] + 1
        }

        boolean morphologyPass = fractionPass && connectedPass && enrichmentPass && compartmentPass
        def finalCall = ""
        String callStatus
        def failureReasons = []
        if (!fractionPass) failureReasons << "insufficient_spatial_coverage"
        if (!connectedPass) failureReasons << "fragmented_spatial_pattern"
        if (!enrichmentPass) failureReasons << (c.role == "nuc_ratio" ? "nuc_cyto_ratio_below_minimum" : "nuclear_enrichment_below_minimum")
        if (!compartmentPass && compartmentAssigned) failureReasons << "wrong_compartment"

        if (!cfg.morphologyPrimary) {
          finalCall = intensityPos ? 1 : 0
          callStatus = intensityPos ? "legacy_intensity_positive" : "legacy_intensity_negative"
        } else if (!evaluable) {
          finalCall = ""
          callStatus = "indeterminate"
          indeterminateCount[m] = indeterminateCount[m] + 1
        } else {
          finalCall = morphologyPass ? 1 : 0
          boolean fixedThreshold = chThreshSource[m] == "fixed_predeclared"
          callStatus = morphologyPass ?
                       (fixedThreshold ? "positive" : "exploratory_positive") :
                       (fixedThreshold ? "negative" : "exploratory_negative")
        }
        if (finalCall == 1) {
          finalPosCount[m] = finalPosCount[m] + 1
          finalPositiveRois[m] << nuc
        } else if (finalCall == 0) {
          finalNegCount[m] = finalNegCount[m] + 1
        } else {
          indeterminateRois[m] << nuc
        }
        calls[m] = finalCall
        row[m + "_morphology_pass"] = (evaluable && morphologyPass) ? 1 : 0
        row[m + "_final_call"] = finalCall
        row[m + "_true_pos"] = finalCall       // compatibility alias
        row[m + "_call_status"] = callStatus
        row[m + "_call_reason"] = (indeterminateReasons + failureReasons).unique().join(";")
      }
      // classifications
      panelDef.classify.each { rule ->
        def key = rule.collect { mk, want -> mk + (want ? "+" : "-") }.join("_")
        boolean classEvaluable = rule.every { mk, want -> calls[mk] == 0 || calls[mk] == 1 }
        if (!classEvaluable) {
          row["class_" + key] = ""
          row["class_" + key + "_status"] = "indeterminate"
        } else {
          boolean ok = rule.every { mk, want -> calls[mk] == (want ? 1 : 0) }
          row["class_" + key] = ok ? 1 : 0
          row["class_" + key + "_status"] = ok ? "positive" : "negative"
          classEvaluableCount[key] = classEvaluableCount[key] + 1
          if (ok) classCount[key] = classCount[key] + 1
        }
      }
      cellRows << row
    }

    // ---- QC overlay for this region ----
    def qc = buildQcOverlay(markerImg, panelDef, region, allNucRois, qcMasks)
    def qcPath = imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__QC.png"
    IJ.saveAs(qc, "PNG", qcPath); qc.close()
    def dapiQc = buildDapiQc(dapi, region, allNucRois, rejectedNuclei, cfg)
    IJ.saveAs(dapiQc, "PNG", imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__DAPI_QC.png")
    dapiQc.close()
    cellChannels.each { c ->
      def callQc = buildCallDecisionQc(dapi, region, allNucRois,
                                       finalPositiveRois[c.marker], indeterminateRois[c.marker],
                                       c.marker, chThreshSource[c.marker], cfg)
      IJ.saveAs(callQc, "PNG", imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__" + c.marker + "_CALL_QC.png")
      callQc.close()
    }
    if (segmentation.candidateMask != null) {
      IJ.saveAs(segmentation.candidateMask, "Tiff",
                imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__DAPI_candidate_mask.tif")
      segmentation.candidateMask.close()
    }

    // ---- region summary row ----
    def srow = [ image: sourceStem, panel: panelKey, region: regName,
                 mouse_id: meta.mouse_id, section_id: meta.section_id,
                 genotype: meta.genotype, condition: meta.condition, compartment: compartment,
                 region_tags: regionTags.join("|"),
                 region_area_um2: regionAreaUm2,
                 dapi_segmentation_method: cfg.dapiMethod,
                 n_nuclei: nuclei.size(),
                 n_rejected_nucleus_candidates: rejectedNuclei.size(),
                 n_rejected_below_min_area: rejectedNuclei.count { it.reason == "area_below_minimum" },
                 n_rejected_at_image_edge: rejectedNuclei.count { it.reason == "image_edge" },
                 n_rejected_by_particle_filter: rejectedNuclei.count { it.reason == "particle_filter" } ]
    cellChannels.each { c ->
      // The conventional summary fields follow the authoritative final call.
      // Raw object-mean decisions remain available under explicit audit names.
      srow[c.marker + "_pos_count"] = finalPosCount[c.marker]
      srow[c.marker + "_density_per_mm2"] = (regionAreaUm2 > 0 ? finalPosCount[c.marker] / (regionAreaUm2/1e6) : 0)
      srow[c.marker + "_raw_mean_pos_count"] = posCount[c.marker]
      srow[c.marker + "_raw_mean_density_per_mm2"] = (regionAreaUm2 > 0 ? posCount[c.marker] / (regionAreaUm2/1e6) : 0)
      srow[c.marker + "_pos_threshold"] = chThresh[c.marker]   // resolved raw-intensity cutoff
      srow[c.marker + "_threshold_source"] = chThreshSource[c.marker]
      srow[c.marker + "_measurement_model"] = c.measurement ?: c.role
      srow[c.marker + "_call_authority"] = cfg.morphologyPrimary ? "morphology_primary" : "legacy_mean_intensity"
      srow[c.marker + "_morphology_pos_count"] = finalPosCount[c.marker]
      srow[c.marker + "_morphology_negative_count"] = finalNegCount[c.marker]
      srow[c.marker + "_indeterminate_count"] = indeterminateCount[c.marker]
      srow[c.marker + "_morphology_evaluable_count"] = finalPosCount[c.marker] + finalNegCount[c.marker]
      srow[c.marker + "_morphology_density_per_mm2"] =
        (regionAreaUm2 > 0 ? finalPosCount[c.marker] / (regionAreaUm2/1e6) : 0)
      srow[c.marker + "_true_pos_count"] = finalPosCount[c.marker] // compatibility alias
      def expectedForSummary = expectedCompartmentsFor(c)
      srow[c.marker + "_expected_compartment"] = expectedForSummary.isEmpty() ? "none" : expectedForSummary.join("|")
    }
    areaStats.each { m, as ->
      srow[m + "_positive_area_um2"] = as.area_um2
      srow[m + "_positive_area_frac"] = as.frac_of_region
      srow[m + "_n_components"] = as.n_components
      srow[m + "_mean_component_area_um2"] = as.mean_component_area_um2
      srow[m + "_min_component_area_um2"] = as.min_component_area_um2
      srow[m + "_area_threshold"] = as.threshold
      srow[m + "_area_mode"] = as.mode
      srow[m + "_area_threshold_source"] = as.threshold_source
      def areaChannel = panelDef.channels.find { it.marker == m }
      def areaExpectedCompartments = areaChannel == null ? [] : expectedCompartmentsFor(areaChannel)
      if (!areaExpectedCompartments.isEmpty() &&
          (compartment == "unassigned" || compartment == "ambiguous")) {
        srow[m + "_area_call_status"] = as.threshold_source == "fixed_predeclared" ?
                                         "indeterminate_compartment_unassigned" :
                                         "exploratory_compartment_unassigned"
      } else if (!areaExpectedCompartments.isEmpty() &&
                 !areaExpectedCompartments.any { regionTags.contains(it) }) {
        srow[m + "_area_call_status"] = "wrong_compartment_not_interpretable"
      } else {
        srow[m + "_area_call_status"] = as.threshold_source == "fixed_predeclared" ?
                                         "fixed_threshold_area" :
                                         "exploratory_adaptive_threshold"
      }
      // Keep the historical KRT5 pod fields for downstream compatibility.
      if (as.mode == "pod") {
        srow[m + "_pod_area_um2"] = as.area_um2
        srow[m + "_pod_area_frac"] = as.frac_of_region
        srow[m + "_n_pods"] = as.n_components
        srow[m + "_mean_pod_area_um2"] = as.mean_component_area_um2
        srow[m + "_pod_threshold"] = as.threshold
      }
    }
    (classCount.keySet() + classEvaluableCount.keySet()).unique().each { k ->
      srow["class_" + k + "_count"] = classCount[k]
      srow["class_" + k + "_evaluable_count"] = classEvaluableCount[k]
      srow["class_" + k + "_indeterminate_count"] = nuclei.size() - classEvaluableCount[k]
    }
    summaryRows << srow
    qcMasks.each { k, v -> v.close() }

    // save nuclei mask for the region
    saveLabelMask(dapi, allNucRois, imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__nuclei_mask.tif")
    saveLabelMask(dapi, rejectedNuclei.collect { it.roi }, imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__rejected_nuclei_mask.tif")
    cellChannels.each { c ->
      saveLabelMask(dapi, finalPositiveRois[c.marker],
                    imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__" + c.marker + "_morphology_positive_nuclei_mask.tif")
      saveLabelMask(dapi, indeterminateRois[c.marker],
                    imgOut.getAbsolutePath() + "/" + fileKey + "__" + regName + "__" + c.marker + "_indeterminate_nuclei_mask.tif")
    }
  }

  // Save morphology-specific binary masks with names that describe the unit.
  panelDef.channels.findAll { it.areaMarker }.each { c ->
    def mask = areaMasks[c.marker]
    String areaMode = c.areaMode ?: "pod"
    String suffix = areaMode == "pod" ? "pod_mask" :
                    (areaMode == "ciliary" ? "ciliary_mask" :
                     (areaMode == "membrane" ? "membrane_positive_mask" :
                      (areaMode == "reporter" ? "reporter_positive_mask" : "positive_area_mask")))
    IJ.saveAs(mask, "Tiff", imgOut.getAbsolutePath() + "/" + fileKey + "__" + c.marker + "_" + suffix + ".tif")
  }

  // per-image params/provenance
  def params = [
    image: new File(imgPath).name, output_key: outputKey, channel_signature: channelSignature,
    panel: panelKey, panel_label: panelDef.label,
    calibration: [ pixel_width_um: cal.pixelWidth, pixel_height_um: cal.pixelHeight,
                   pixel_depth_um: cal.pixelDepth, unit: cal.getUnit(),
                   n_slices: raw.getNSlices(), n_channels: channels.length ],
    projection: cfg.projection, single_plane: cfg.singlePlane,
    segmenter: cfg.segmenter, stardist_prob: cfg.prob, stardist_nms: cfg.nms, stardist_tiles: cfg.tiles,
    dapi_preprocessing: [ method: cfg.dapiMethod,
                          background_radius_um: cfg.dapiBackgroundRadiusUm,
                          local_radius_um: cfg.dapiLocalRadiusUm,
                          blur_sigma_px: cfg.dapiBlurSigmaPx,
                          contrast_saturation_percent: cfg.dapiContrastSaturation ],
    ring_expand_um: cfg.ringExpandUm, min_nucleus_area_um2: cfg.minNucArea,
    acetylated_tubulin_model: [ measurement: "apical_cilia_proximity_and_regional_patches",
                                support_expand_um: cfg.actubSupportExpandUm,
                                minimum_support_positive_fraction: cfg.actubMinSupportFraction,
                                minimum_ciliary_patch_area_um2: cfg.actubMinPatchAreaUm2 ],
    pod_min_area_um2: cfg.podMinArea, pod_blur_sigma_px: cfg.podBlur, pod_thresh_method: cfg.podMethod,
    pos_sensitivity: cfg.sensitivity, black_background: cfg.blackBackground,
    fixed_pos_thresholds: cfg.fixedThresholds,
    decision_hierarchy: [ authority: cfg.morphologyPrimary ? "morphology_primary" : "legacy_mean_intensity",
                          call_states: ["positive", "negative", "indeterminate"],
                          intensity_role: "candidate-pixel threshold and audit field; not final-call authority",
                          fixed_threshold_requirement: "confirmatory calls require predeclared control-derived thresholds" ],
    morphology_rules: cfg.morphologyRules,
    role_morphology_defaults: cfg.roleMorphologyDefaults,
    marker_registry: [path:cfg.markerRegistryPath, schema_version:cfg.markerRegistrySchema],
    custom_panel_config: cfg.panelConfigPath,
    custom_panel_keys: cfg.customPanelKeys,
    compartment_mode: cfg.compartmentMode,
    whole_field_compartment: cfg.wholeFieldCompartment,
    minimum_ring_positive_fraction: cfg.minRingPosFraction,
    tissue_mode: cfg.tissueMode, tissue_roi_source: tissue.source,
    tissue_thresh_method: cfg.tissueMethod,
    rejected_nucleus_rules: [minimum_area_um2: cfg.minNucArea, exclude_image_edge: true],
    channel_map: panelDef.channels
  ]
  new File(imgOut, fileKey + "__params.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(params))

  // write per-image cell CSV
  writeCsv(cellRows, imgOut.getAbsolutePath() + "/" + fileKey + "__cells.csv")

  // cleanup
  markerImg.each { k, v -> v.close() }
  areaMasks.each { k, v -> v.close() }
  raw.close()
  return [summary: summaryRows, cells: cellRows.size(), tissue_source: tissue.source,
          channel_signature: channelSignature]
}

// ============================================================================
//  6. OUTPUT HELPERS
// ============================================================================

def markerOutlineColor(String marker) {
  switch (marker) {
    case "T1A": case "CC10": return new Color(0, 255, 0)
    case "tdTOM": return new Color(255, 40, 40)
    case "mRAGE": case "AcTub": return Color.WHITE
    case "KRT5": return Color.MAGENTA
    default: return new Color(255, 180, 0)
  }
}

// Additively merge display-normalized channels into a headless-safe RGB image.
// For panel R: DAPI=blue, T1A=green, tdTOM=red, mRAGE=white.
def buildQcComposite(markerImg, panelDef) {
  def first = markerImg.values().iterator().next()
  int w = first.getWidth(), h = first.getHeight()
  def layers = panelDef.channels.collect { c ->
    def ip = markerImg[c.marker].getProcessor().duplicate()
    ip.resetMinAndMax()
    return [color: (c.qcColor ?: "white"), ip: ip.convertToByte(true)]
  }
  def out = new ColorProcessor(w, h)
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      int rr = 0, gg = 0, bb = 0
      layers.each { layer ->
        int v = layer.ip.get(x, y)
        switch (layer.color) {
          case "red": rr = Math.max(rr, v); break
          case "green": gg = Math.max(gg, v); break
          case "blue": bb = Math.max(bb, v); break
          default:
            rr = Math.max(rr, v); gg = Math.max(gg, v); bb = Math.max(bb, v)
        }
      }
      out.set(x, y, (rr << 16) | (gg << 8) | bb)
    }
  }
  return new ImagePlus("four_channel_QC", out)
}

def buildQcOverlay(markerImg, panelDef, Roi region, nucRois, channelMasks) {
  def rgb = buildQcComposite(markerImg, panelDef)
  ImageProcessor cp = rgb.getProcessor()

  // Continuous fluorescence-region boundaries replace the former per-cell
  // marker circles. Color follows acquisition: green/red/white.
  panelDef.channels.findAll { it.role != "nuclear" }.each { c ->
    def pm = channelMasks[c.marker].duplicate()
    pm.getProcessor().setThreshold(128, 255, ImageProcessor.NO_LUT_UPDATE)
    def signalRoi = ThresholdToSelection.run(pm)
    cp.setLineWidth(c.marker == "mRAGE" ? 2 : 1)
    cp.setColor(markerOutlineColor(c.marker))
    if (signalRoi != null) {
      def clipped = new ShapeRoi(signalRoi).and(new ShapeRoi(region))
      clipped.drawPixels(cp)
    }
    pm.close()
  }

  // Orange is analysis-region membership; it is never an exclusion symbol.
  cp.setLineWidth(2); cp.setColor(new Color(255, 150, 0)); region.drawPixels(cp)

  // Cyan is the only per-object circle: nuclei included in DAPI-based counts.
  cp.setLineWidth(2); cp.setColor(new Color(0, 220, 255))
  nucRois.each { it.drawPixels(cp) }

  // Burn a compact legend into the exported PNG so colors remain auditable.
  cp.setMask(null); cp.setColor(Color.BLACK); cp.setRoi(0, 0, Math.min(cp.getWidth(), 1190), 64); cp.fill()
  cp.resetRoi(); cp.setMask(null); cp.setFont(new Font("SansSerif", Font.BOLD, 16)); cp.setColor(Color.WHITE)
  def rawLegend = panelDef.channels.collect { c -> c.marker + " " + (c.qcColor ?: "gray") }.join(" | ")
  cp.drawString("RAW: " + rawLegend, 10, 22)
  cp.drawString("OUTLINES: counted DAPI nuclei cyan | ROI orange | thresholded marker regions", 10, 48)
  return rgb
}

// DAPI-only audit view. Contrast balancing affects this PNG only; segmentation
// uses the separately recorded preprocessing path and candidate mask.
def buildDapiQc(ImagePlus dapi, Roi region, nucRois, rejectedNuclei, cfg) {
  def display = dapi.duplicate()
  IJ.run(display, "Enhance Contrast...", "saturated=" + cfg.dapiContrastSaturation + " normalize")
  if (display.getBitDepth() != 8) IJ.run(display, "8-bit", "")
  def rgb = new ImagePlus("DAPI_segmentation_QC", display.getProcessor().convertToRGB())
  display.close()
  ImageProcessor cp = rgb.getProcessor()
  cp.setLineWidth(2); cp.setColor(new Color(255, 150, 0)); region.drawPixels(cp)
  cp.setLineWidth(2); cp.setColor(new Color(210, 0, 255))
  rejectedNuclei.each { it.roi.drawPixels(cp) }
  cp.setLineWidth(2); cp.setColor(new Color(0, 220, 255))
  nucRois.each { it.drawPixels(cp) }
  cp.setMask(null); cp.setColor(Color.BLACK); cp.setRoi(0, 0, Math.min(cp.getWidth(), 1050), 58); cp.fill()
  cp.resetRoi(); cp.setMask(null); cp.setFont(new Font("SansSerif", Font.BOLD, 16)); cp.setColor(Color.WHITE)
  cp.drawString("DAPI ONLY (display-balanced): counted cyan | rejected candidate violet | ROI orange", 10, 23)
  cp.drawString("Segmentation method: " + cfg.dapiMethod, 10, 47)
  return rgb
}

// Marker-specific decision audit: negative/evaluable nuclei are cyan, final
// positives are green, and indeterminate nuclei are magenta. This directly
// visualizes the hierarchy that downstream classifications consume.
def buildCallDecisionQc(ImagePlus dapi, Roi region, allNucRois,
                        positiveRois, indeterminateRois,
                        String marker, String thresholdSource, cfg) {
  def display = dapi.duplicate()
  IJ.run(display, "Enhance Contrast...", "saturated=" + cfg.dapiContrastSaturation + " normalize")
  if (display.getBitDepth() != 8) IJ.run(display, "8-bit", "")
  def rgb = new ImagePlus(marker + "_call_QC", display.getProcessor().convertToRGB())
  display.close()
  ImageProcessor cp = rgb.getProcessor()
  cp.setLineWidth(2); cp.setColor(new Color(255, 150, 0)); region.drawPixels(cp)
  cp.setLineWidth(2); cp.setColor(new Color(0, 210, 255)); allNucRois.each { it.drawPixels(cp) }
  cp.setLineWidth(3); cp.setColor(new Color(60, 255, 70)); positiveRois.each { it.drawPixels(cp) }
  cp.setLineWidth(3); cp.setColor(new Color(230, 0, 255)); indeterminateRois.each { it.drawPixels(cp) }
  cp.setMask(null); cp.setColor(Color.BLACK); cp.setRoi(0, 0, Math.min(cp.getWidth(), 1200), 61); cp.fill()
  cp.resetRoi(); cp.setMask(null); cp.setFont(new Font("SansSerif", Font.BOLD, 16)); cp.setColor(Color.WHITE)
  cp.drawString(marker + " FINAL CALLS: positive green | negative cyan | indeterminate magenta | ROI orange", 10, 23)
  cp.drawString("Threshold source: " + thresholdSource + " | authority: morphology", 10, 48)
  return rgb
}

def saveLabelMask(ImagePlus ref, nucRois, String path) {
  // 16-bit labels so >255 nuclei per region do not collide.
  def ip = new ShortProcessor(ref.getWidth(), ref.getHeight())
  nucRois.eachWithIndex { r, i ->
    ip.setValue(i + 1); ip.fill(r)
  }
  def lab = new ImagePlus("labels", ip); lab.setCalibration(ref.getCalibration())
  IJ.saveAs(lab, "Tiff", path); lab.close()
}

def writeCsv(rows, String path) {
  if (rows == null || rows.isEmpty()) { new File(path).text = ""; return }
  def cols = [] as LinkedHashSet
  rows.each { r -> cols.addAll(r.keySet()) }
  cols = cols as List
  def sb = new StringBuilder()
  sb.append(cols.join(",")).append("\n")
  rows.each { r ->
    sb.append(cols.collect { c ->
      def v = r.containsKey(c) ? r[c] : ""
      def s = (v == null) ? "" : v.toString()
      (s.contains(",") || s.contains("\"")) ? "\"" + s.replace("\"", "\"\"") + "\"" : s
    }.join(",")).append("\n")
  }
  new File(path).text = sb.toString()
}

// ============================================================================
//  7. SAMPLESHEET / METADATA
// ============================================================================

def parseMeta(String fname, sheet, defaultPanel, String sourceContext = "") {
  // sheet: map filename -> [mouse_id, section_id, genotype, condition, panel]
  if (sheet != null && sheet.containsKey(fname)) return sheet[fname]
  // fallback: try token convention mouseID_condition_panel_section.ext
  def stem = fname.replaceFirst(/\.[^.]+$/, "")
  // 260719-CW convention. Use the parent path as context because stitched
  // overview names (for example Stitch_A01_G001.oir) omit mouse/stain details.
  def identityText = stem + " " + sourceContext
  def cw = (identityText =~ /(?i).*?\b(\d{6})\s+([MF]\d+)\s+(pr8_(?:bleo|PBS))\b/)
  def section = (stem =~ /(?i).*(A\d+_G\d+(?:_\d+)?)$/)
  if (cw.find() && section.matches()) {
    def lower = identityText.toLowerCase()
    def inferredPanel = lower.contains("4x mapping") && lower.contains("cc10_488") ? "M" :
                        lower.contains("cc10_488") ? "E" :
                        (lower.contains("t1a_488") ? "R" : defaultPanel)
    return [ mouse_id: cw.group(1) + "_" + cw.group(2),
             section_id: section.group(1), genotype: "krt5-creERT2;tdTOM",
             condition: cw.group(3).toLowerCase(), panel: inferredPanel ]
  }
  def toks = stem.split("_")
  return [ mouse_id: (toks.length > 0 ? toks[0] : "NA"),
           condition: (toks.length > 1 ? toks[1] : "NA"),
           panel: (toks.length > 2 ? toks[2] : defaultPanel),
           section_id: (toks.length > 3 ? toks[3] : "NA"),
           genotype: "NA" ]
}

def loadSamplesheet(String dir) {
  def f = new File(dir, "samplesheet.csv")
  if (!f.exists()) return null
  def lines = f.readLines()
  if (lines.isEmpty()) return null
  def header = lines[0].split(",").collect { it.trim() }
  def idx = { name -> header.indexOf(name) }
  def map = [:]
  lines.drop(1).each { ln ->
    if (ln.trim().isEmpty() || ln.trim().startsWith("#")) return   // skip blanks + comments
    def p = ln.split(",")
    def get = { n -> def i = idx(n); (i >= 0 && i < p.length) ? p[i].trim() : "NA" }
    def fn = get("filename")
    if (fn == null || fn.isEmpty() || fn == "NA") return
    map[fn] = [ mouse_id: get("mouse_id"), section_id: get("section_id"),
                genotype: get("genotype"), condition: get("condition"),
                panel: get("panel") ]
  }
  return map
}

// ============================================================================
//  8. MAIN
// ============================================================================

// Force a consistent binary convention: thresholded objects become 255 on a 0
// background. Without this, "Convert to Mask" can invert per the user's Fiji
// prefs and silently corrupt every area/count. Recorded in provenance below.
Prefs.blackBackground = true

ensureDir(OUTPUT_DIR)
def cfg = [ segmenter: SEGMENTER, prob: STARDIST_PROB, nms: STARDIST_NMS, tiles: STARDIST_TILES,
           dapiMethod: DAPI_METHOD, dapiBackgroundRadiusUm: DAPI_BACKGROUND_RADIUS_UM,
           dapiLocalRadiusUm: DAPI_LOCAL_RADIUS_UM, dapiBlurSigmaPx: DAPI_BLUR_SIGMA_PX,
           dapiContrastSaturation: DAPI_CONTRAST_SATURATION,
           blackBackground: true,
           projection: PROJECTION, singlePlane: SINGLE_PLANE,
           ringExpandUm: RING_EXPAND_UM, minNucArea: MIN_NUCLEUS_AREA_UM2,
           podMinArea: POD_MIN_AREA_UM2, podBlur: POD_BLUR_SIGMA_PX, podMethod: POD_THRESH_METHOD,
           sensitivity: POS_SENSITIVITY, fixedThresholds: FIXED_POS_THRESHOLDS,
           morphologyPrimary: MORPHOLOGY_PRIMARY, morphologyRules: MORPHOLOGY_RULES,
           roleMorphologyDefaults: ROLE_MORPHOLOGY_DEFAULTS,
           markerRegistryPath: markerRegistryFile.isFile() ? markerRegistryFile.getAbsolutePath() : "unavailable",
           markerRegistrySchema: MARKER_REGISTRY.schema_version ?: "unavailable",
           panelConfigPath: PANEL_CONFIG_PATH ?: "built_in_only", customPanelKeys: CUSTOM_PANEL_KEYS,
           minRingPosFraction: MIN_RING_POS_FRACTION, compartmentMode: COMPARTMENT_MODE,
           wholeFieldCompartment: WHOLE_FIELD_COMPARTMENT,
           actubSupportExpandUm: ACTUB_SUPPORT_EXPAND_UM,
           actubMinSupportFraction: ACTUB_MIN_SUPPORT_FRACTION,
           actubMinPatchAreaUm2: ACTUB_MIN_PATCH_AREA_UM2,
           tissueMode: TISSUE_MODE, tissueBlur: TISSUE_BLUR_SIGMA_PX,
           tissueMethod: TISSUE_THRESH_METHOD, tissueMinArea: TISSUE_MIN_AREA_UM2 ]

def versions = captureVersions()
IJ.log("ImageJ " + versions.imagej_version + " | Bio-Formats " + versions.bioformats_version)

def inDir = new File(INPUT_DIR)
if (!inDir.isDirectory()) {
  IJ.error("INPUT_DIR is not a folder:\n" + INPUT_DIR +
           "\nEdit INPUT_DIR / OUTPUT_DIR at the top of the script.")
  return
}

def sheet = USE_SAMPLESHEET ? loadSamplesheet(INPUT_DIR) : null

def listed = []
if (RECURSIVE) inDir.eachFileRecurse { f -> if (f.isFile()) listed << f }
else listed = (inDir.listFiles() ?: [] as File[]).toList()
def includePattern = ~/(?i)${INCLUDE_REGEX}/
def files = listed.findAll { it.isFile() && (it.name ==~ FILE_GLOB) && (it.getAbsolutePath() ==~ includePattern) }
                  .sort { it.getAbsolutePath() }
if (MAX_IMAGES > 0) files = files.take(MAX_IMAGES)
IJ.log("Found " + files.size() + " image(s).")
if (files.isEmpty()) IJ.log("  (nothing matched FILE_GLOB in INPUT_DIR)")

def masterSummary = []
def manifest = [ run_timestamp: versions.timestamp, versions: versions, config: cfg,
                 input_dir: INPUT_DIR, output_dir: OUTPUT_DIR, recursive: RECURSIVE,
                 include_regex: INCLUDE_REGEX, max_images: MAX_IMAGES, images: [] ]

// Duplicate basenames occur in folders such as Cycle and Cycle_01. Give only
// duplicates a stable parent-path suffix so a recursive run cannot overwrite
// an earlier image's output directory.
def basenameCounts = files.countBy { it.name }

files.each { f ->
  def m = parseMeta(f.name, sheet, PANEL, f.parentFile.absolutePath)
  String requestedPanel = m.panel == null ? "" : m.panel.toString().trim()
  def panelKey = (requestedPanel.isEmpty() || requestedPanel == "NA") ? PANEL : requestedPanel
  if (!PANELS.containsKey(panelKey)) {
    throw new IllegalArgumentException("Image '" + f.name + "' requests unknown panel '" + panelKey +
                                       "'. Available panels: " + PANELS.keySet().sort())
  }
  def panelDef = PANELS[panelKey]
  // Keep Windows paths manageable. The source acquisition name is preserved
  // in params/manifest; output paths use concise, analysis-relevant metadata.
  def safeToken = { value ->
    def s = (value == null || value.toString().trim().isEmpty()) ? "NA" : value.toString().trim()
    return s.replaceAll(/[^A-Za-z0-9._-]+/, "-").replaceAll(/^-+|-+$/, "")
  }
  def baseStem = [m.mouse_id, m.condition, panelKey, m.section_id].collect { safeToken(it) }.join("_")
  def outputKey = (basenameCounts[f.name] > 1) ?
                  baseStem + "__" + Integer.toHexString(f.parentFile.absolutePath.hashCode()) : baseStem
  try {
    def res = processImage(f.getAbsolutePath(), outputKey, panelKey, panelDef, m, cfg, OUTPUT_DIR)
    if (res != null) {
      masterSummary.addAll(res.summary)
      manifest.images << [ file: f.name, relative_path: inDir.toPath().relativize(f.toPath()).toString(),
                           output_key: outputKey, panel: panelKey, channel_signature: res.channel_signature,
                           tissue_source: res.tissue_source, n_cells: res.cells ]
    }
  } catch (Throwable t) {
    IJ.log("  ERROR on " + f.name + ": " + t)
    def sw = new java.io.StringWriter(); t.printStackTrace(new java.io.PrintWriter(sw)); IJ.log(sw.toString())
    manifest.images << [ file: f.name, panel: panelKey, error: t.getMessage() ]
  }
}

// master summary + manifest
writeCsv(masterSummary, OUTPUT_DIR + "/run_summary.csv")
new File(OUTPUT_DIR, "run_manifest.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(manifest))

IJ.log("DONE. Wrote run_summary.csv and run_manifest.json to " + OUTPUT_DIR)
IJ.log("Reminder: aggregate run_summary.csv to MOUSE level before stats (n = mice, not sections).")

// ImageJ starts non-daemon UI/event threads even with --headless. Exit after
// synchronous exports so command-line and cluster jobs do not hang at DONE.
if (java.awt.GraphicsEnvironment.isHeadless()) System.exit(0)
