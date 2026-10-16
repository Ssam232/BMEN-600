"""Shared settings for the BMEN 600 ADNI pipeline.

Every script imports this file, so paths and analysis choices are set in one place.
Raw ADNI files go in data/raw/ (never committed - see data/README.md).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"          # ADNI/IDA downloads (git-ignored)
DERIVED = ROOT / "data" / "derived"  # subject-level outputs (git-ignored)
OUTPUTS = ROOT / "outputs"           # aggregate tables and figures (safe to commit)

# ---- IDA exports (Advanced Image Search -> CSV), one per diagnostic group ----
IDA_EXPORTS = {
    "AD": RAW / "BMEN_AD_10_09_2026.csv",
    "MCI": RAW / "BMEN_MCI_10_09_2026.csv",
    "CN": RAW / "BMEN_CN_10_09_2026.csv",
}
# ---- UPenn Roche Elecsys CSF biomarkers (ADNI download) ----
CSF_FILE = RAW / "UPENNBIOMK_ROCHE_ELECSYS_06Oct2026.csv"

# ---- Matching rules (script 01) ----
AGE_CALIPER = 3            # controls must be within +/- this many years
CONTROLS_PER_CASE = 2      # 2 controls per AD and per MCI participant
# Participants kept with a single control because only one eligible control exists
# (team decision, Oct 2026).
ALLOW_ONE_CONTROL = ["032_S_6602"]

# ---- Analysis choices (scripts 03 and 04) ----
# FreeSurfer hippocampal subfield labels used (head/body pieces are summed).
SUBFIELDS = {
    "CA1": "CA1",
    "CA3": "CA3",            # FreeSurfer's CA3 label includes CA2
    "CA4": "CA4",
    "DG": "GC-ML-DG",        # dentate gyrus (granule cell + molecular layer)
    "subiculum": "subiculum",
    "presubiculum": "presubiculum",
}
WHOLE = {"whole_hippocampus": "Whole_hippocampus"}

# CSF markers analysed: column name -> label. t-tau and p-tau181 are log-transformed in 03.
MARKERS = {"ABETA42": "Abeta42", "log_TAU": "log t-tau", "log_PTAU": "log p-tau181"}
COVARIATES = ["Age", "sex_male", "eTIV"]

# Roche Elecsys measuring ranges (pg/mL); values outside are flagged in the CSF file COMMENT column.
ASSAY_RANGE = {"ABETA42": (200, 1700), "TAU": (80, 1300), "PTAU": (8, 120)}

N_BOOT = 5000
SEED = 600
FDR_Q = 0.05
MIN_N = 8   # skip a group-level estimate with fewer participants than this
