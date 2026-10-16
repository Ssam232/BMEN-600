# Analysis code

Pipeline for: *Do CSF Abeta42, t-tau and p-tau181 correlate with hippocampal subfield volumes, and does the
association differ between CN, MCI and AD?* (ADNI3, year-1 T1-weighted MPRAGE, HippUnfold subfields).

Study sample: `data/raw/ADNI_Participants.xlsx` - 24 AD, 41 MCI, 47 CN, one scan each.

| Step | File | What it does | Output |
|---|---|---|---|
| 1 | `01_select_and_match.py` | Matches each AD participant to 2 CN controls (same sex, age +/-3 y); documents how the controls were chosen | `data/derived/subjects.csv`, `outputs/step1_matching/` |
| 2 | `02_preprocessing.ipynb` (+ `preprocess_hippunfold.py`) | DICOM -> NIfTI -> BIDS, HippUnfold command, subfield volumes from HippUnfold label maps, eTIV, QC sheet | `data/derived/volumes.csv`, `data/derived/qc.csv` |
| 3 | `03_merge.py` | One row per participant (volumes, eTIV, QC, nearest CSF sample); Table 1 | `data/derived/analysis_dataset.csv`, `outputs/table1_demographics.md`, `outputs/merge_report.md` |
| 4 | `04_analysis.py` | Statistics, tables, heatmap | `outputs/results/*.csv`, `outputs/figures/*.png` |

All settings (file names, image folders, subfields, markers, covariates, bootstrap count, seed) are in `config.py`.
If the images are on another drive, change `DICOM_DIR`, `BIDS_DIR` and `HIPPUNFOLD_DIRS` there.

## Run

```bash
pip install -r requirements.txt
```

Then either open `run_pipeline.ipynb` (steps 1, 3, 4) and `code/02_preprocessing.ipynb` (step 2) and run all
cells, or from the `BMEN-600_code` folder:

```bash
python code/01_select_and_match.py
# step 2: run code/02_preprocessing.ipynb; run HippUnfold in a terminal with the command it prints;
#         then rerun its cells D-F to collect volumes, add eTIV and create the QC sheet
python code/03_merge.py
python code/04_analysis.py
```

## What each step does

**1 - Matching.** Collapses the IDA image lists to one row per participant and matches each AD participant to
2 CN controls of the same sex within +/-3 years, solved exactly as an integer linear program
(`scipy.optimize.milp`, HiGHS). No control is reused; AD 032_S_6602 is kept with its only eligible control. Every
pair is checked automatically. The team then fixed the final sample in `ADNI_Participants.xlsx`: the 24 AD
participants, their 47 matched controls, and 41 MCI participants compared with the same controls.

**2 - Preprocessing** For every participant in the spreadsheet:
converts DICOMs (or copies an already-converted NIfTI) to `BIDS_DIR/sub-XXX/anat/sub-XXX_T1w.nii.gz` without moving
or deleting source files; prints the HippUnfold command (`--modality T1w`); measures each subfield's volume per
hemisphere by counting voxels in HippUnfold's `space-cropT1w` subfield label map (bigbrain atlas: subiculum,
CA1, CA2, CA3, CA4, dentate gyrus; whole hippocampus = these + SRLM, cysts excluded); adds eTIV from a
`Subject,eTIV` CSV (`ST10CV`, matched on image ID); creates `qc.csv` for the
visual check of HippUnfold's QC images.

**3 - Merge.** One row per participant: demographics, left/right/bilateral volumes, eTIV, QC flag, and the CSF
sample closest in time to the MRI (gap recorded). Hemolyzed samples dropped, out-of-range values flagged, t-tau
and p-tau181 log-transformed. Writes Table 1 for the final sample.

**4 - Analysis.** Within-group partial correlations (age, sex, eTIV adjusted) with 5,000-resample bootstrap CIs
and Benjamini-Hochberg FDR; pooled `volume ~ marker * group + age + sex + eTIV` models testing whether
associations differ by stage; whole-hippocampus comparison; ANCOVA of group differences; exploratory
marker x sex interaction; sensitivity analyses (in-range CSF only, left/right hemispheres, Spearman).
Participants marked `QC_pass = 0` are excluded.

## Data privacy

Only aggregate outputs (counts, means, test statistics, figures) are committed. Images, the participant
spreadsheet and all subject-level files stay in `data/` and are git-ignored (ADNI Data Use Agreement). Clear
notebook outputs before committing if a cell displays participant-level data.

## AI assistance

Steps 1, 3 and 4, the step 2 helper module and both notebooks were written with Claude (Anthropic; model
Claude Opus 5.5, `claude-opus-5-5`), building on the team's original preprocessing notebook. Step 1
reproduces the participant lists the team used; Table 1 from step 3 matches the team's `demographics.md`.
Steps 2-4 were tested end-to-end on synthetic images and volumes only and must be checked on the first real
HippUnfold outputs (in particular the label numbers in `config.HIPPUNFOLD_LABELS`).
