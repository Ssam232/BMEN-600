# Data (not in git)

ADNI data cannot be redistributed, so nothing in `data/raw/` or `data/derived/` is committed.
Each team member downloads the files from ADNI/IDA and places them here.

## `data/raw/` - downloads

| File | Source |
|---|---|
| `BMEN_AD_10_09_2026.csv` | IDA Advanced Image Search export: ADNI3, AD, year-1 visit, MPRAGE |
| `BMEN_MCI_10_09_2026.csv` | same search, MCI group |
| `BMEN_CN_10_09_2026.csv` | same search, CN group |
| `UPENNBIOMK_ROCHE_ELECSYS_06Oct2026.csv` | ADNI download: UPenn Roche Elecsys CSF biomarkers |

File names are set in `code/config.py`; update them there if you export new versions.

## `data/derived/` - created by the scripts

| File | Made by | Contents |
|---|---|---|
| `subjects.csv` | step 1 | every participant with KEEP/REMOVE, matched controls and the reason |
| `matches_AD.csv`, `matches_MCI.csv` | step 1 | each case with its controls |
| `volumes.csv` *or* a FreeSurfer SUBJECTS_DIR | step 2 (preprocessing) | subfield volumes and eTIV - format in `code/03_merge.py` |
| `qc.csv` | step 2 | `Subject, QC_pass` (1 = pass, 0 = fail) |
| `analysis_dataset.csv` | step 3 | one row per participant: demographics, volumes, CSF |
