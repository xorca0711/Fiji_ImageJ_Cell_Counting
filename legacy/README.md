# Legacy Archive

This folder contains superseded development artifacts. Nothing here is an
authoritative analysis entry point.

The complete repository layout immediately before organization is preserved on
Git branch:

`codex/legacy-pre-reorganization`

## Contents

- `figures/`: earlier workflow and overlay diagrams. These predate the final
  morphology-authoritative three-state hierarchy.
- `scripts/README.md`: routes historical script access to the snapshot branch
  and commit history; no legacy script is an active entry point.
- `test_runs/`: local smoke, baseline, and intermediate Fiji outputs.
- `pilot_output/`: earliest sample-image pilot output.
- `ref_images_analysis_output/`: output formerly stored beside the reference
  ND2 image.

Generated run folders and test scripts remain ignored by Git. This README and
the archived tracked figures document what was moved.

## Run chronology

- `TestRun1` through `TestRun9`: import, channel-map, DAPI, naming, and initial
  T1alpha/mRAGE development.
- `TestRun10`: CC10/AcTub baseline before morphology authority.
- `TestRun11`: initial morphology implementation.
- `TestRun12`: two-image morphology-primary validation before final area-mask
  consistency improvements.
- `FinalPilot_*_morphology_primary`: first one-image confirmation runs.
- The final `*_v2` runs are not legacy; they remain in `../test_runs/current/`.

For current work, start at `../WORKFLOW.md` and
`../IF_Quant_Pipeline.groovy`.
