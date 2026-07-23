# Configuration assets

- [`lung_marker_registry.json`](lung_marker_registry.json) stores descriptive
  marker aliases, localization, analytical-role defaults, lineage/state notes,
  and research-context cautions. It is not a positivity cutoff table or
  diagnostic classifier.
- [`custom_panels.example.json`](custom_panels.example.json) demonstrates
  opt-in image channel maps for acute injury, IPF, stromal fibrosis,
  Red2-KrasG12D RFP/Ki-67/SOX9 analysis, and lung adenocarcinoma lineage
  research.

Copy the example to a study-controlled location before editing it. Load that
copy with `IFQ_PANEL_CONFIG`; do not alter the built-in panel definitions for a
new acquisition.

See
[`docs/UNIVERSAL_MARKER_CONFIGURATION.md`](../docs/UNIVERSAL_MARKER_CONFIGURATION.md)
for schema, ROI vocabulary, and validation rules.
