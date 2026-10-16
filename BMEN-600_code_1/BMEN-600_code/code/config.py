"""Shared settings for the BMEN 600 ADNI pipeline.

Every script and notebook imports this file, so paths and analysis choices are set in one place.
Data live under data/ (never committed - see data/README.md). If your images are on another drive,
change DICOM_DIR / BIDS_DIR / HIPPUNFOLD_DIRS below (e.g. Path("D:/BMEN_600/ADNI")).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"          # ADNI/IDA downloads (git-ignored)
DERIVED = ROOT / "data" / "derived"  # subject-level outputs (git-ignored)
OUTPUTS = ROOT / "outputs"           # aggregate tables and figures (safe to commit)

# ---- Final study sample (one row per participant / scan), set by the team ----
# Columns used: Image Data ID, Acq Date, Subject, ID (sub-001...), Group, Sex, Age
PARTICIPANTS_FILE = RAW / "ADNI_Participants.xlsx"

# ---- IDA exports (Advanced Image Search -> CSV), one per diagnostic group (step 1 only) ----
IDA_EXPORTS = {
    "AD": RAW / "BMEN_AD_10_09_2026.csv",
    "MCI": RAW / "BMEN_MCI_10_09_2026.csv",
    "CN": RAW / "BMEN_CN_10_09_2026.csv",
}
# ---- UPenn Roche Elecsys CSF biomarkers (ADNI download) ----
CSF_FILE = RAW / "UPENNBIOMK_ROCHE_ELECSYS_06Oct2026.csv"

# ---- Step 2: image folders (change to your drive if needed) ----
DICOM_DIR = ROOT / "data" / "dicom"            # DICOM_DIR/<Group>/<sub-XXX or ADNI subject>/.../*.dcm
BIDS_DIR = ROOT / "data" / "bids"              # HippUnfold input: BIDS_DIR/sub-XXX/anat/sub-XXX_T1w.nii.gz
HIPPUNFOLD_DIRS = [ROOT / "data" / "hippunfold"]  # HippUnfold output folder(s); list several if run in batches
# Intracranial volume: CSV with columns Subject,eTIV  (or an ADNI UCSFFSX7 table, see 02_preprocessing)
ETIV_FILE = RAW / "etiv.csv"

# ---- Matching rules (step 1) ----
AGE_CALIPER = 3            # controls must be within +/- this many years
CONTROLS_PER_CASE = 2      # 2 controls per AD and per MCI participant
# Participants kept with a single control because only one eligible control exists (team decision, Oct 2026).
ALLOW_ONE_CONTROL = ["032_S_6602"]

# ---- HippUnfold subfields (bigbrain atlas) ----
# Default label numbers of the HippUnfold subfield map; overridden by the label table HippUnfold writes
# next to the map (*_dseg.tsv) when present.
HIPPUNFOLD_LABELS = {1: "Sub", 2: "CA1", 3: "CA2", 4: "CA3", 5: "CA4", 6: "DG", 7: "SRLM", 8: "Cyst"}
# Analysis name -> HippUnfold label name
SUBFIELDS = {"subiculum": "Sub", "CA1": "CA1", "CA2": "CA2", "CA3": "CA3", "CA4": "CA4", "DG": "DG"}
# Whole hippocampus = all hippocampal labels (cysts excluded)
WHOLE = {"whole_hippocampus": ["Sub", "CA1", "CA2", "CA3", "CA4", "DG", "SRLM"]}

# ---- Analysis choices (steps 3 and 4) ----
# CSF markers analysed: column name -> label. t-tau and p-tau181 are log-transformed in step 3.
MARKERS = {"ABETA42": "Abeta42", "log_TAU": "log t-tau", "log_PTAU": "log p-tau181"}
COVARIATES = ["Age", "sex_male", "eTIV"]

# Roche Elecsys measuring ranges (pg/mL); values outside are flagged in the CSF file COMMENT column.
ASSAY_RANGE = {"ABETA42": (200, 1700), "TAU": (80, 1300), "PTAU": (8, 120)}

N_BOOT = 5000
SEED = 600
FDR_Q = 0.05
MIN_N = 8   # skip a group-level estimate with fewer participants than this
