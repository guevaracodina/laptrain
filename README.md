# Longitudinal fNIRS of cognitive workload in laparoscopic training across a 24-hour shift

Analysis code for the Neurophotonics data descriptor **NPH-260108-1**,
*A Longitudinal fNIRS Dataset of Cognitive Workload in Laparoscopic Training
Across a 24-Hour Clinical Shift* (Guevara, Torres Cuevas, Avalos Martínez,
Kolosovas-Machuca, Martínez-Jiménez).

The dataset is on Zenodo: <https://doi.org/10.5281/zenodo.15186569>

Thirty General Surgery residents performed a standardized laparoscopic training
task (Origami Box Folding Exercise) and an eyes-closed resting-state recording
at 0, 12 and 24 h of a single continuous duty cycle, while prefrontal
hemodynamics were recorded with a 24-channel Artinis Brite MKII.

---

## What is in here

| folder | contents |
|---|---|
| `code/` | the analysis pipeline as originally submitted |
| `code/rev01/` | the revised pipeline; **this is what reproduces the published numbers** |
| `code/BIDS/` | conversion scripts used to build the released dataset |
| `data/` | empty; the Zenodo archive is unpacked here (see `data/README.md`) |

`code/rev01/` takes precedence over `code/` wherever a filename appears in both.
Put it **first** on the MATLAB path; `rev01_run_all` refuses to start otherwise.

---

## Requirements

MATLAB R2023a or later, with the **Statistics and Machine Learning Toolbox**
(`fitlme`, `fitrm`, `multcompare`, `ttest`, `ranksum`).

External toolboxes:

| toolbox | version used | needed for | source |
|---|---|---|---|
| Homer3 | 1.87.0 | preprocessing, SNIRF I/O | <https://github.com/BUNPC/Homer3> |
| Brain Connectivity Toolbox | 2019-03-03 | graph metrics | <https://sites.google.com/site/bctnet/> |
| QT-NIRS | — | signal-quality metrics | <https://github.com/lpollonini/qt-nirs> |
| NBS | 1.2 | network-based statistic (GUI) | <https://www.nitrc.org/projects/nbs/> |
| shadedErrorBar | — | figures only | <https://github.com/raacampbell/shadedErrorBar> |
| FieldTrip | 20250402 | `code/BIDS/` only | <https://www.fieldtriptoolbox.org/> |

Only Homer3, BCT and the Statistics Toolbox are required to reproduce the
reported statistics.

### Configure once

```matlab
cd code
copyfile('setup_paths_config_example.m', 'setup_paths_config.m')
edit setup_paths_config.m        % point each field at your own copy
setup_paths                      % prints which toolboxes it can find
```

`setup_paths_config.m` is gitignored, so local paths are never committed.

---

## Reproducing the published results

Unpack the Zenodo archive into `data/` and extract the six Homer3 derivative
archives (`data/BIDSsource/<cond>derivatives.7z`) so that
`data/<cond>/derivatives/homer/sub-XX/nirs/*.mat` exists. See `data/README.md`
for the full expected tree.

Then, from `code/`:

```matlab
reproduce_all
```

That runs the whole pipeline in order, checking each stage's inputs before it
starts. To run the stages by hand instead, the order matters and is:

| # | command | produces | notes |
|---|---|---|---|
| 0 | `setup_paths` | — | then `addpath(fullfile(pwd,'rev01'),'-begin')` |
| 1 | `qt_nirs_script` | `data/qtnirs/<cond>/QC_reportTable_<cond>.mat` | slow; needs QT-NIRS |
| 2 | `qc_summary_export` | `data/qc_summary.mat` | condenses stage 1 (~142 MB → ~60 kB) |
| 3 | `load_motion_artifacts_script` | `data/motion_artifacts_summary.mat` | slow |
| 4 | `rev01_run_all` | `data/rev01_results.mat`, `.txt` | the reported statistics |
| 5 | `rev01_pruning_report` | pruning and exclusion counts | printed |
| 6 | `NBS_script` | `data/conn_Mat_HbO_NBS_<cond>.mat` + design and exchange | then run NBS in its GUI |

**`rev01` must be first on the MATLAB path.** It shadows `conn_regress_global`
and `NBS_script` in `code/`, and the corrected versions are the ones that
reproduce the published numbers. `rev01_run_all` and `reproduce_all` both
refuse to start otherwise.

Stages 1 and 3 read every derivative file and are the slow ones; they only need
running once. Inside `rev01_run_all`, `doRecompute = true` re-derives the
connectivity matrices from the derivatives — necessary the first time, and
skippable afterwards.

The network-based statistic is not scripted end to end: NBS has no documented
batch entry point, so stage 6 writes its inputs and the test is run in the NBS
GUI with the F-test, 5000 permutations, and the contrast saved in
`data/NBS_subjects_<cond>.mat`.

### Runtime

Stages 1 and 3 take tens of minutes each. Stage 4 is dominated by the
permutation tests: six condition-by-chromophore combinations at 5000
permutations, roughly an hour on a current desktop.

### What you should get

`data/rev01_results.txt` reports, for laparoscopic training, a binary
clustering coefficient AUC of 0.098, 0.114 and 0.100 at 0, 12 and 24 h
(F = 5.74, p = 0.008 for HbO; F = 6.84, p = 0.005 for HbR; not significant for
HbT), and 20 of 30 participants complete at all three time points. Permutation
p-values vary in the third decimal between runs; the rng seed is fixed at 42 in
`rev01_run_all`, so an identical MATLAB version should reproduce them exactly.

## What changed in revision 01, and why it matters for reuse

Four corrections were made during peer review. Anyone reusing the submitted
code should be aware of all four.

**1. Short-separation channel indices.** `conn_regress_global.m` and
`NBS_script.m` used channels `[3, 14]`, which are the short channels of a
*motor* montage used in an earlier study. In this prefrontal montage they are
**3 (S03D01) and 24 (S10D08)**. The indices are now derived from
`get_channels_from_template` rather than hard-coded.

**2. The clustering coefficient was not computed on binarized graphs.** BCT's
`clustering_coef_bu` requires binary input; it was being passed the weighted
output of `threshold_proportional`, so `sum(S)/(k²−k)` summed correlation
weights instead of counting triangles, and the metric scaled linearly with
connectivity strength. Scaling every edge by 1.05 changed it by exactly 5.00%.
`rev01_network_measures.m` now binarizes first, and returns both the corrected
metric (`AUCclust`) and the original quantity (`AUCclustW`) so the difference
is auditable.

**3. Independent-sample tests on a within-subject design.** `permutationTest3`
and `permutationTest` pool all 90 observations and permute freely. The three
visits are repeated measures on the same 30 residents.
`permutationTestRM3` permutes condition labels *within* each participant, and
`permutationTestRMpair` sign-flips within-subject differences. In simulation
with realistic between-subject variance the free scheme had power 0.005 against
0.990 for the restricted one, with type I error correctly calibrated at 0.040.

**4. Incomplete subjects in the NBS design.** The submitted `NBS_script` built
a 30-subject design matrix before loading any data, then asserted every cell
held a 24×24 matrix — which fails as soon as a run is dropped by the 120 s
minimum. It now selects subjects complete at all three time points and builds
the array, design matrix and exchange blocks together (N = 20 for training,
N = 30 for resting-state).

### Known limitations of this code

- `swi.m` computes the small-world index on **weighted** graphs against a
  random reference that is neither degree- nor density-preserving. It is kept
  unchanged for continuity with the submitted analysis and should not be
  interpreted quantitatively.
- `conn_mat_subject_level.m` contains a `doPreWhitening` branch that calls
  `prewhiten`, which is not included here. The branch is disabled
  (`doPreWhitening = false`) and was never used.
- The scripts in `code/` (outside `rev01/`) use Windows path separators and
  were run on Windows. The `rev01/` scripts use `fullfile` throughout and are
  platform-independent.

---

## Montage

24 channels; **3 (S03D01)** and **24 (S10D08)** are short-separation (1.5 cm),
the other 22 are long (3 cm). `code/rev01/get_channel_info_prefrontal.m` is the
single source of truth for channel index, source, detector, MNI coordinate and
anatomical label, and `code/rev01/optodes_MNI_coordinates.csv` is generated
from it. Channel order is detector-major, matching the SNIRF measurement list.

---

## Third-party files included here

`permutationTest.m` and `permutationTest3.m` are by Laurens R. Krol; `swtest.m`
is the Shapiro-Wilk test from the MATLAB File Exchange. Their original
copyright headers are intact. They are redistributed for completeness of the
analysis record; their own licence terms apply to them, not the licence below.

## Licence

MIT, see `LICENSE`. The dataset itself is CC BY 4.0.

## Citation

See `CITATION.cff`.
