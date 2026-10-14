# Analysis code

Pipeline for: *Do CSF Abeta42, t-tau and p-tau181 correlate with hippocampal subfield volumes, and does
the association differ between CN, MCI and AD?* (ADNI3, year-1 T1 MPRAGE).

| Step | Script | Status | Output |
|---|---|---|---|
| 1 | `01_select_and_match.py` | done | keep/remove lists (`data/derived/`), `outputs/table1_demographics.md`, `outputs/matching_summary.md` |
| 2 | preprocessing (FreeSurfer, subfield segmentation, QC) | to be added by the preprocessing lead | `data/derived/volumes.csv` (or SUBJECTS_DIR) + `qc.csv` |
| 3 | `03_merge.py` | ready; waiting on step 2 | `data/derived/analysis_dataset.csv`, `outputs/merge_report.md` |
| 4 | `04_analysis.py` | ready; waiting on step 2 | `outputs/results/*.csv`, `outputs/figures/partial_correlation_heatmap.png` |

All settings (file names, age caliper, subfields, markers, covariates, bootstrap count, seed) live in `config.py`.

## Run

```bash
pip install -r requirements.txt
python code/01_select_and_match.py
# step 2 (preprocessing) goes here
python code/03_merge.py                       # or: python code/03_merge.py --fs-dir /path/to/SUBJECTS_DIR
python code/04_analysis.py
```

## What each step does

**1 - Selection and matching.** Collapses the IDA image lists to one row per participant, then matches
each AD and each MCI participant to 2 CN controls of the same sex within +/-3 years. Matching is solved
exactly as an integer linear program (`scipy.optimize.milp`, HiGHS): keep as many cases fully matched as
possible, prefer cases with CSF data, then minimise total age difference. No control is reused within a
comparison; one control can serve one AD and one MCI participant. AD 032_S_6602 is kept with its only
eligible control (team decision). Every pair is checked automatically, and each participant gets a written
reason for being kept or removed.

**3 - Merge.** One row per participant: demographics, bilateral and per-hemisphere subfield volumes, eTIV,
QC flag, and the CSF sample closest in time to the MRI (gap recorded). Hemolyzed samples are dropped,
out-of-range CSF values are flagged, t-tau and p-tau181 are log-transformed.

**4 - Analysis.** Within-group partial correlations (age, sex, eTIV adjusted) with 5,000-resample bootstrap
CIs and Benjamini-Hochberg FDR; pooled `volume ~ marker * group + age + sex + eTIV` models testing whether
associations differ by stage; whole-hippocampus comparison; ANCOVA of group differences; exploratory
marker x sex interaction; sensitivity analyses (unmatched sample, in-range CSF only, left/right
hemispheres, Spearman).

## Data privacy

Only aggregate outputs (tables of counts, means, test statistics, figures) are committed. Subject-level
files stay in `data/` and are git-ignored, as required by the ADNI Data Use Agreement.

## AI assistance

Steps 1, 3 and 4 were written with Claude (Anthropic; model Claude Opus 5.5, `claude-opus-5-5`). Step 1
reproduces the participant lists the team used to edit the IDA collections; steps 3 and 4 were tested
end-to-end on synthetic volumes only and must be re-checked once real FreeSurfer outputs exist.
