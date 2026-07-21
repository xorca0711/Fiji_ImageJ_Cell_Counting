/*
 * ============================================================================
 *  IF_Quant_Pipeline.groovy
 *  Confocal immunofluorescence quantification for the IFN-gamma KO / PR8
 *  influenza injury project (KRT5 pod remodeling readout).
 * ----------------------------------------------------------------------------
 *  Antibody set: KRT5, Pro-SPC, AGER, PDPN, CD4, CD8, Sox2  (+ p63, YAP, Aqp5)
 *
 *  Panels (max 3 markers/slide = DAPI + 2 primaries). PANEL keys are single
 *  tokens so they survive filename parsing; pick per slide via samplesheet:
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
 *   * DEFAULT THRESHOLDS ARE PLACEHOLDERS. Auto-Otsu adapts per image, but you
 *     MUST confirm positivity calls against the QC overlays and set per-marker
 *     sensitivity before reporting. Freeze parameters once, then batch.
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
def OUTPUT_DIR  = envOr("IFQ_OUTPUT_DIR", new File("pilot_output").getAbsolutePath())
def PANEL       = envOr("IFQ_PANEL", "T")
def FILE_GLOB   = ~/(?i).*\.(czi|lif|nd2|oir|oib|oif|ics|tif|tiff)$/
def RECURSIVE   = envOr("IFQ_RECURSIVE", "false").toBoolean()
def INCLUDE_REGEX = envOr("IFQ_INCLUDE_REGEX", ".*")
def MAX_IMAGES  = Integer.parseInt(envOr("IFQ_MAX_IMAGES", "0")) // 0 = all
def TISSUE_MODE = envOr("IFQ_TISSUE_MODE", "auto").toLowerCase() // auto | whole_field

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
def MIN_NUCLEUS_AREA_UM2 = 15.0 // classic segmenter min object size
def POD_MIN_AREA_UM2    = 50.0  // a "pod" particle must exceed this
def POD_BLUR_SIGMA_PX   = 2.0
def POD_THRESH_METHOD   = "Otsu" // Otsu|Triangle|Li|Huang|MaxEntropy...

// --- Positivity ---
// Object is positive if its mean (raw) >= channelOtsu(inTissue) * sensitivity.
// Tune per marker: >1 = stricter, <1 = more permissive.
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
//           "cyto"      -> measured in perinuclear ring (KRT5, Pro-SPC, CD8, CD4)
//           "membrane"  -> measured in ring; interpret as area too (AGER/PDPN)
//           "nuc_marker"-> measured in the nucleus (p63)
//           "nuc_ratio" -> nucleus vs ring separately (YAP)  [needs single plane]
//     areaMarker: also run independent threshold AREA quantification (KRT5 pod)
// ============================================================================

def PANELS = [
  "A": [ label:"A_KRT5_AGER",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear"],
               [idx:2, marker:"KRT5",  role:"cyto",     areaMarker:true],
               [idx:3, marker:"AGER",  role:"membrane"] ],
    classify:[ ["KRT5":true,"AGER":false], ["KRT5":true,"AGER":true] ] ],

  "B": [ label:"B_KRT5_ProSPC",
    channels:[ [idx:1, marker:"DAPI",   role:"nuclear"],
               [idx:2, marker:"KRT5",   role:"cyto",    areaMarker:true],
               [idx:3, marker:"ProSPC", role:"cyto"] ],
    classify:[ ["KRT5":true,"ProSPC":false], ["KRT5":false,"ProSPC":true] ] ],

  "C": [ label:"C_KRT5_CD8",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:3, marker:"CD8",  role:"cyto"] ],
    classify:[ ["CD8":true], ["KRT5":true,"CD8":true] ] ],

  "D": [ label:"D_KRT5_CD4",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:3, marker:"CD4",  role:"cyto"] ],
    classify:[ ["CD4":true], ["KRT5":true,"CD4":true] ] ],

  // AT1 alternative via podoplanin (T1-alpha). Enables the KRT5+/PDPN- readout.
  "P": [ label:"P_KRT5_PDPN",
    channels:[ [idx:1, marker:"DAPI", role:"nuclear"],
               [idx:2, marker:"KRT5", role:"cyto",     areaMarker:true],
               [idx:3, marker:"PDPN", role:"membrane"] ],
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
               [idx:2, marker:"CC10",  role:"cyto",    qcColor:"green"],
               [idx:3, marker:"tdTOM", role:"cyto",    qcColor:"red", areaMarker:true] ],
    classify:[ ["tdTOM":true], ["CC10":true,"tdTOM":true] ] ],

  "E": [ label:"E_CC10_tdTOM_AcTub",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear", qcColor:"blue"],
               [idx:2, marker:"CC10",  role:"cyto",    qcColor:"green"],
               [idx:3, marker:"tdTOM", role:"cyto",    qcColor:"red", areaMarker:true],
               [idx:4, marker:"AcTub", role:"cyto",    qcColor:"white"] ],
    classify:[ ["tdTOM":true], ["CC10":true,"tdTOM":true],
               ["AcTub":true,"tdTOM":true] ] ],

  "R": [ label:"R_T1A_tdTOM_mRAGE",
    channels:[ [idx:1, marker:"DAPI",  role:"nuclear",  qcColor:"blue"],
               [idx:2, marker:"T1A",   role:"membrane", qcColor:"green"],
               [idx:3, marker:"tdTOM", role:"cyto",     qcColor:"red", areaMarker:true],
               [idx:4, marker:"mRAGE", role:"membrane", qcColor:"white"] ],
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
  def pa = new ParticleAnalyzer(opts, Measurements.AREA, rt,
                                minAreaCal, Double.MAX_VALUE)
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

  ImageProcessor lp = labels.getProcessor()
  def out = []
  for (int i = 1; i <= rt.size(); i++) {
    lp.setThreshold((double)i, (double)i, ImageProcessor.NO_LUT_UPDATE)
    def r = ThresholdToSelection.run(new ImagePlus("particle_labels", lp))
    if (r != null) out << r
  }
  lp.resetThreshold()
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

// Build a binary mask ImagePlus (255 = signal) from a channel by threshold.
// The numeric lower threshold is stashed as a property for provenance export.
// Requires Prefs.blackBackground=true (set in main) so foreground = 255.
def buildThresholdMask(ImagePlus ch, double blurSigma, String method) {
  ImagePlus dup = ch.duplicate()
  if (blurSigma > 0) new GaussianBlur().blurGaussian(dup.getProcessor(), blurSigma)
  IJ.setAutoThreshold(dup, method + " dark")
  double thr = dup.getProcessor().getMinThreshold()   // raw intensity, -1 if unset
  IJ.run(dup, "Convert to Mask", "")
  dup.setCalibration(ch.getCalibration())
  dup.setProperty("thresholdValue", thr)
  return dup
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
    // classic watershed fallback
    def m = crop.duplicate()
    new GaussianBlur().blurGaussian(m.getProcessor(), 2.0)
    IJ.setAutoThreshold(m, "Otsu dark")
    IJ.run(m, "Convert to Mask", "")
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
    m.close()
    crop.close()
    return [included: rois, rejected: rejected]
  }
}

// ============================================================================
//  5. PER-IMAGE PROCESSING
// ============================================================================

def processImage(String imgPath, String outputKey, panelKey, panelDef, meta, cfg, outDir) {
  IJ.log("---- " + new File(imgPath).name + "  [panel " + panelKey + "] ----")
  def sourceStem = new File(imgPath).name.replaceFirst(/\.[^.]+$/, "")
  def imgOut = ensureDir(outDir + "/" + outputKey)

  def raw = bfOpen(imgPath)
  Calibration cal = raw.getCalibration()
  def nChExpected = panelDef.channels.size()
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

  // Precompute per-marker positivity thresholds (Otsu within whole field first;
  // refined per tissue region below).
  def tissue = resolveTissueRois(imgPath, dapi, cfg)

  // Pre-build KRT5 (and any areaMarker) threshold masks once for pod area.
  def areaMasks = [:]
  panelDef.channels.findAll { it.areaMarker }.each { c ->
    areaMasks[c.marker] = buildThresholdMask(markerImg[c.marker], cfg.podBlur, cfg.podMethod)
  }

  def cellRows = []      // per-object records (all regions)
  def summaryRows = []   // per-region summary
  def qcRegionOverlays = [:]

  tissue.regions.eachWithIndex { reg, ri ->
    def region = reg.roi
    def regName = reg.name
    def regStat = measureRoi(dapi, region)   // area of region
    double regionAreaUm2 = regStat.area

    // per-marker channel thresholds inside THIS region (adaptive)
    def chThresh = [:]
    panelDef.channels.findAll { it.role != "nuclear" }.each { c ->
      double t = autoThresholdInRoi(markerImg[c.marker], region, "Otsu")
      double sens = (cfg.sensitivity[c.marker] ?: 1.0)
      chThresh[c.marker] = t * sens
    }

    // ---- Pod area (areaMarkers, e.g. KRT5) ----
    def podStats = [:]
    panelDef.channels.findAll { it.areaMarker }.each { c ->
      def mask = areaMasks[c.marker]
      double podArea = positiveAreaInRoi(mask, region)
      // discrete pods: particle analysis of the mask restricted to region
      def maskReg = mask.duplicate(); maskReg.setCalibration(cal); maskReg.setRoi(region)
      def podRois = particlesToRois(maskReg, cfg.podMinArea, false)
      def podAreas = podRois.collect { measureRoi(mask, it).area }
      def podThr = mask.getProperty("thresholdValue")
      podStats[c.marker] = [ area_um2: podArea,
                             frac_of_region: (regionAreaUm2 > 0 ? podArea/regionAreaUm2 : 0),
                             n_pods: podRois.size(),
                             mean_pod_area_um2: (podAreas.isEmpty()? 0 : podAreas.sum()/podAreas.size()),
                             threshold: (podThr != null ? podThr : -1) ]
      maskReg.close()
    }

    // ---- Nuclei -> cells ----
    def segmentation = segmentNuclei(dapi, region, cfg)
    def nuclei = segmentation.included
    def rejectedNuclei = segmentation.rejected ?: []
    def posCount = [:].withDefault { 0 }
    def classCount = [:].withDefault { 0 }
    def posRois = [:].withDefault { [] }
    def allNucRois = []

    nuclei.eachWithIndex { nuc, ni ->
      allNucRois << nuc
      def cellRoi = RoiEnlarger.enlarge(nuc, (int)Math.round(cfg.ringExpandUm / cal.pixelWidth))
      def row = [ image: sourceStem, panel: panelKey, region: regName, cell_id: (ni + 1) ]
      row.mouse_id = meta.mouse_id; row.section_id = meta.section_id
      row.genotype = meta.genotype; row.condition = meta.condition
      def cs = measureRoi(dapi, nuc)
      row.centroid_x_um = cs.cx; row.centroid_y_um = cs.cy; row.nucleus_area_um2 = cs.area

      def calls = [:]
      panelDef.channels.findAll { it.role != "nuclear" }.each { c ->
        def m = c.marker
        def img = markerImg[m]
        double val, nucVal = 0, ringVal = 0
        if (c.role == "nuc_marker") {
          val = measureRoi(img, nuc).mean
        } else if (c.role == "nuc_ratio") {
          nucVal = measureRoi(img, nuc).mean
          // true cytoplasmic ring = (enlarged cell) MINUS (nucleus)
          def ring = ringOnly(nuc, cellRoi)
          double cytoVal = (ring != null) ? measureRoi(img, ring).mean : measureRoi(img, cellRoi).mean
          val = nucVal   // positivity uses nuclear signal for YAP
          row[m + "_nuc_mean"] = nucVal
          row[m + "_cyto_mean"] = cytoVal
          row[m + "_nuc_cyto_ratio"] = (cytoVal > 0 ? nucVal/cytoVal : 0)
        } else { // cyto / membrane: measure the perinuclear ring, not nucleus + ring
          def ring = ringOnly(nuc, cellRoi)
          val = (ring != null) ? measureRoi(img, ring).mean : measureRoi(img, cellRoi).mean
        }
        row[m + "_mean"] = val
        boolean pos = val >= chThresh[m]
        calls[m] = pos
        row[m + "_pos"] = pos ? 1 : 0
        if (pos) {
          posCount[m] = posCount[m] + 1
          posRois[m] << cellRoi
        }
      }
      // classifications
      panelDef.classify.each { rule ->
        boolean ok = rule.every { mk, want -> (calls[mk] == want) }
        def key = rule.collect { mk, want -> mk + (want ? "+" : "-") }.join("_")
        row["class_" + key] = ok ? 1 : 0
        if (ok) classCount[key] = classCount[key] + 1
      }
      cellRows << row
    }

    // ---- QC overlay for this region ----
    def qc = buildQcOverlay(markerImg, panelDef, region, allNucRois, rejectedNuclei, posRois,
                            (areaMasks.isEmpty()? null : areaMasks[areaMasks.keySet().iterator().next()]))
    def qcPath = imgOut.getAbsolutePath() + "/" + outputKey + "__" + regName + "__QC.png"
    IJ.saveAs(qc, "PNG", qcPath); qc.close()

    // ---- region summary row ----
    def srow = [ image: sourceStem, panel: panelKey, region: regName,
                 mouse_id: meta.mouse_id, section_id: meta.section_id,
                 genotype: meta.genotype, condition: meta.condition,
                 region_area_um2: regionAreaUm2,
                 n_nuclei: nuclei.size(),
                 n_rejected_nucleus_candidates: rejectedNuclei.size(),
                 n_rejected_below_min_area: rejectedNuclei.count { it.reason == "area_below_minimum" },
                 n_rejected_at_image_edge: rejectedNuclei.count { it.reason == "image_edge" },
                 n_rejected_by_particle_filter: rejectedNuclei.count { it.reason == "particle_filter" } ]
    panelDef.channels.findAll { it.role != "nuclear" }.each { c ->
      srow[c.marker + "_pos_count"] = posCount[c.marker]
      srow[c.marker + "_density_per_mm2"] = (regionAreaUm2 > 0 ? posCount[c.marker] / (regionAreaUm2/1e6) : 0)
      srow[c.marker + "_pos_threshold"] = chThresh[c.marker]   // resolved raw-intensity cutoff
    }
    podStats.each { m, ps ->
      srow[m + "_pod_area_um2"] = ps.area_um2
      srow[m + "_pod_area_frac"] = ps.frac_of_region
      srow[m + "_n_pods"] = ps.n_pods
      srow[m + "_mean_pod_area_um2"] = ps.mean_pod_area_um2
      srow[m + "_pod_threshold"] = ps.threshold
    }
    classCount.each { k, v -> srow["class_" + k + "_count"] = v }
    summaryRows << srow

    // save nuclei mask for the region
    saveLabelMask(dapi, allNucRois, imgOut.getAbsolutePath() + "/" + outputKey + "__" + regName + "__nuclei_mask.tif")
    saveLabelMask(dapi, rejectedNuclei.collect { it.roi }, imgOut.getAbsolutePath() + "/" + outputKey + "__" + regName + "__rejected_nuclei_mask.tif")
  }

  // save pod masks
  areaMasks.each { m, mask ->
    IJ.saveAs(mask, "Tiff", imgOut.getAbsolutePath() + "/" + outputKey + "__" + m + "_pod_mask.tif")
  }

  // per-image params/provenance
  def params = [
    image: new File(imgPath).name, output_key: outputKey,
    panel: panelKey, panel_label: panelDef.label,
    calibration: [ pixel_width_um: cal.pixelWidth, pixel_height_um: cal.pixelHeight,
                   pixel_depth_um: cal.pixelDepth, unit: cal.getUnit(),
                   n_slices: raw.getNSlices(), n_channels: channels.length ],
    projection: cfg.projection, single_plane: cfg.singlePlane,
    segmenter: cfg.segmenter, stardist_prob: cfg.prob, stardist_nms: cfg.nms, stardist_tiles: cfg.tiles,
    ring_expand_um: cfg.ringExpandUm, min_nucleus_area_um2: cfg.minNucArea,
    pod_min_area_um2: cfg.podMinArea, pod_blur_sigma_px: cfg.podBlur, pod_thresh_method: cfg.podMethod,
    pos_sensitivity: cfg.sensitivity, black_background: cfg.blackBackground,
    tissue_mode: cfg.tissueMode, tissue_roi_source: tissue.source,
    tissue_thresh_method: cfg.tissueMethod,
    rejected_nucleus_rules: [minimum_area_um2: cfg.minNucArea, exclude_image_edge: true],
    channel_map: panelDef.channels
  ]
  new File(imgOut, outputKey + "__params.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(params))

  // write per-image cell CSV
  writeCsv(cellRows, imgOut.getAbsolutePath() + "/" + outputKey + "__cells.csv")

  // cleanup
  markerImg.each { k, v -> v.close() }
  areaMasks.each { k, v -> v.close() }
  raw.close()
  return [summary: summaryRows, cells: cellRows.size(), tissue_source: tissue.source]
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

def buildQcOverlay(markerImg, panelDef, Roi region, nucRois, rejectedNuclei, posRois, podMask) {
  def rgb = buildQcComposite(markerImg, panelDef)
  ImageProcessor cp = rgb.getProcessor()

  // Yellow is the continuous area-marker mask (tdTOM in panels E/M/R).
  if (podMask != null) {
    def pm = podMask.duplicate(); pm.setRoi(region)
    cp.setLineWidth(1); cp.setColor(Color.YELLOW)
    particlesToRois(pm, 0.0d, false).each { it.drawPixels(cp) }
    pm.close()
  }

  // Orange is analysis-region membership; it is never an exclusion symbol.
  cp.setLineWidth(2); cp.setColor(new Color(255, 150, 0)); region.drawPixels(cp)

  // Violet candidates were found by the DAPI Otsu/watershed mask but rejected
  // by minimum calibrated area or image-edge rules.
  cp.setLineWidth(2); cp.setColor(new Color(220, 0, 255))
  rejectedNuclei.each { it.roi.drawPixels(cp) }

  // Cyan nuclei are the objects included in the cell table and summary count.
  cp.setLineWidth(1); cp.setColor(new Color(0, 220, 255))
  nucRois.each { it.drawPixels(cp) }

  // Positive perinuclear measurement rings use their acquisition colors.
  panelDef.channels.findAll { it.role != "nuclear" }.each { c ->
    cp.setLineWidth(c.marker == "mRAGE" ? 3 : 2)
    cp.setColor(markerOutlineColor(c.marker))
    (posRois[c.marker] ?: []).each { it.drawPixels(cp) }
  }

  // Burn a compact legend into the exported PNG so colors remain auditable.
  cp.setColor(Color.BLACK); cp.setRoi(0, 0, Math.min(cp.getWidth(), 1190), 64); cp.fill()
  cp.resetRoi(); cp.setFont(new Font("SansSerif", Font.BOLD, 16)); cp.setColor(Color.WHITE)
  cp.drawString("RAW: DAPI blue | T1A green | tdTOM red | mRAGE white", 10, 22)
  cp.drawString("OUTLINES: ROI orange | counted nucleus cyan | rejected DAPI violet | positive rings use marker color", 10, 48)
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
           blackBackground: true,
           projection: PROJECTION, singlePlane: SINGLE_PLANE,
           ringExpandUm: RING_EXPAND_UM, minNucArea: MIN_NUCLEUS_AREA_UM2,
           podMinArea: POD_MIN_AREA_UM2, podBlur: POD_BLUR_SIGMA_PX, podMethod: POD_THRESH_METHOD,
           sensitivity: POS_SENSITIVITY,
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
  def panelKey = (m.panel && PANELS.containsKey(m.panel)) ? m.panel : PANEL
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
                           output_key: outputKey, panel: panelKey,
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
