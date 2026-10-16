# Data (not in git)

ADNI data cannot be redistributed, so nothing in this folder except this README is committed.
Each team member downloads the files from ADNI/IDA and places them here. All paths are set in `code/config.py`.

## `data/raw/` - downloads and team files

| File | Source |
|---|---|
| `ADNI_Participants.xlsx` | Final study sample (team): Image Data ID, Acq Date, Subject, ID (sub-001...), Group, Sex, Age |
| `UPENNBIOMK_ROCHE_ELECSYS_06Oct2026.csv` | ADNI download: UPenn Roche Elecsys CSF biomarkers |
| `BMEN_AD_10_09_2026.csv`, `BMEN_MCI_10_09_2026.csv`, `BMEN_CN_10_09_2026.csv` | IDA Advanced Image Search exports (ADNI3, year-1, MPRAGE) - step 1 only |
| `etiv.csv` *or* an ADNI `UCSFFSX7` table | Intracranial volume: `Subject,eTIV`, or ADNI FreeSurfer table with `PTID, IMAGEUID, ST10CV` |

## Image folders (step 2) - can live on another drive (set in `config.py`)

| Folder | Contents |
|---|---|
| `data/dicom/` (`DICOM_DIR`) | `<Group>/<sub-XXX>/.../*.dcm` as downloaded (or already-converted `sub-XXX_*T1w.nii.gz`) |
| `data/bids/` (`BIDS_DIR`) | HippUnfold input, created by step 2: `sub-XXX/anat/sub-XXX_T1w.nii.gz` |
| `data/hippunfold/` (`HIPPUNFOLD_DIRS`) | HippUnfold output |

## `data/derived/` - created by the scripts

| File | Made by | Contents |
|---|---|---|
| `subjects.csv`, `matches_AD.csv`, `matches_MCI.csv` | step 1 | matching results with reasons |
| `volumes.csv` | step 2 | `Subject, ID, <region>_L, <region>_R, eTIV` |
| `qc.csv` | step 2 (filled in by hand) | `Subject, ID, Group, QC_pass (1/0), QC_notes` |
| `analysis_dataset.csv` | step 3 | one row per participant: demographics, volumes, CSF |
