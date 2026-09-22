# Data

Not tracked in git. Download from Zenodo:
<https://doi.org/10.5281/zenodo.15186569>

## Expected layout

The scripts are run from `code/` and reference `../data`. After unpacking,
and after extracting the Homer3 derivative archives, the tree should be:

```
data/
  lap01/  lap02/  lap03/
  resting01/  resting02/  resting03/
      derivatives/homer/sub-XX/nirs/sub-XX_task-<cond>_nirs.mat
  BIDSsource/
      <cond>/sub-XX/nirs/sub-XX_task-<cond>_nirs.snirf
      participants.tsv
      optodes_MNI_coordinates.csv
      residents_lap_train_tasks_timing.csv
  participants/
      participants.tsv
      demographics.xlsx
      residents_lap_train_tasks_timing.csv
  qtnirs/<cond>/QC_reportTable_<cond>.mat
  motion_artifacts_summary.mat
```

The derivatives ship as `<cond>derivatives.7z` inside `BIDSsource/`; extract
each to `data/<cond>/derivatives/homer/`.

## Generated during analysis

`qc_summary.mat`, `rev01_results.mat`, `rev01_results.txt`,
`*Conn{HbO,HbR,HbT}.mat`, `*_all_significant_connections_*.mat`,
`conn_Mat_HbO_NBS_*.mat`, `design_NBS_*.mat`, `exchange_NBS_*.mat`,
`NBS_subjects_*.mat`. All are gitignored.

## participants.tsv

Columns: `participant_id`, `age`, `sex`, `handedness`, `residency_year`,
`PSCscore01-03`, `GSRscore01-03`.

`residency_year` is the year of the General Surgery programme (1-4).
Laparoscopic training begins in year 1, so it indexes cumulative exposure.

`GSRscore01-03` is the Global Rating Scale. It was recorded **once per
participant** as a measure of general operative skill and is therefore
identical across the three columns by design, not by transcription error.
