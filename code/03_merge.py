"""Step 3 - merge hippocampal subfield volumes and CSF biomarkers into one analysis table.

Runs after step 1 (matching) and step 2 (preprocessing, owned by another team member).

Volumes can come from either source (whichever step 2 produces):
  A) a CSV (default data/derived/volumes.csv) with one row per participant:
        Subject, eTIV, CA1_L, CA1_R, CA3_L, CA3_R, ..., whole_hippocampus_L, whole_hippocampus_R
     (column stems = keys of config.SUBFIELDS and config.WHOLE)
  B) a FreeSurfer SUBJECTS_DIR (--fs-dir), one folder per ADNI subject ID, containing
        mri/lh.hippoSfVolumes-T1.v2*.txt, mri/rh.hippoSfVolumes-T1.v2*.txt  (segmentHA_T1.sh)
        stats/aseg.stats                                                    (eTIV)
Optional QC file (default data/derived/qc.csv): columns Subject, QC_pass (1 = pass, 0 = fail).

CSF rules:
  * one CSF sample per participant: the exam closest in time to the year-1 MRI scan
    (time gap recorded in CSF_gap_days)
  * samples marked "Sample Hemolyzed" are not used
  * values outside the Roche Elecsys measuring range are flagged (CSF_out_of_range);
    below-range t-tau/p-tau (missing in the file) are set to the lower limit;
    above-range Abeta42 keeps the extrapolated value supplied by UPenn
  * t-tau and p-tau181 are log-transformed (log_TAU, log_PTAU)

Output: data/derived/analysis_dataset.csv  - exactly one row per participant
        outputs/merge_report.md            - aggregate counts only

Run:  python code/03_merge.py                      (volumes CSV)
      python code/03_merge.py --fs-dir /path/to/SUBJECTS_DIR
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

import config as cfg

REGIONS = {**cfg.SUBFIELDS, **cfg.WHOLE}


# --------------------------------------------------------------------------- volumes
def parse_hippo_file(path):
    """FreeSurfer hippoSfVolumes file ('label volume' per line) -> {region_key: mm3}."""
    vals = {}
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            vals[parts[0]] = float(parts[1])
    out = {}
    for key, label in REGIONS.items():
        if label in vals:                       # e.g. Whole_hippocampus, or non-split labels
            out[key] = vals[label]
        else:                                   # head/body-split files: sum the pieces
            pieces = [v for k, v in vals.items() if re.sub(r"-(head|body)$", "", k) == label]
            out[key] = sum(pieces) if pieces else np.nan
    return out


def read_etiv(aseg_stats):
    m = re.search(r"EstimatedTotalIntraCranialVol, eTIV, [^,]*, ([\d.]+)", Path(aseg_stats).read_text())
    return float(m.group(1)) if m else np.nan


def volumes_from_freesurfer(fs_dir, subjects):
    rows = []
    for sid in subjects:
        sdir = Path(fs_dir) / sid
        if not sdir.exists():
            continue
        row = {"Subject": sid}
        for hemi, suf in (("lh", "L"), ("rh", "R")):
            files = [f for f in sorted((sdir / "mri").glob(f"{hemi}.hippoSfVolumes-T1.v2*.txt"))
                     if re.fullmatch(rf"{hemi}\.hippoSfVolumes-T1\.v2\d+\.txt", f.name)]
            if files:
                for k, v in parse_hippo_file(files[-1]).items():
                    row[f"{k}_{suf}"] = v
        aseg = sdir / "stats" / "aseg.stats"
        row["eTIV"] = read_etiv(aseg) if aseg.exists() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- CSF
def nearest_csf(csf, subjects_scan):
    csf = csf[~csf.COMMENT.fillna("").str.contains("hemolyz", case=False)].copy()
    csf["EXAMDATE"] = pd.to_datetime(csf.EXAMDATE)
    out = []
    for sid, scan in subjects_scan.items():
        x = csf[csf.PTID == sid]
        if x.empty:
            continue
        gap = (x.EXAMDATE - scan).dt.days
        r = x.loc[gap.abs().idxmin()]
        out.append({"Subject": sid, "CSF_examdate": r.EXAMDATE.date(), "CSF_visit": f"{r.PHASE} {r.VISCODE2}",
                    "CSF_gap_days": int((r.EXAMDATE - scan).days), "ABETA42": r.ABETA42, "TAU": r.TAU,
                    "PTAU": r.PTAU, "CSF_comment": r.COMMENT if isinstance(r.COMMENT, str) else ""})
    d = pd.DataFrame(out)
    com = d.CSF_comment.str.lower()
    lo, hi = cfg.ASSAY_RANGE["TAU"][0], cfg.ASSAY_RANGE["PTAU"][0]
    d["TAU"] = np.where(d.TAU.isna() & com.str.contains(r"(?:^|[^p])tau<", regex=True), lo, d.TAU)
    d["PTAU"] = np.where(d.PTAU.isna() & com.str.contains("ptau<"), hi, d.PTAU)
    d["CSF_out_of_range"] = com.str.contains("[<>]", regex=True)
    d["log_TAU"], d["log_PTAU"] = np.log(d.TAU), np.log(d.PTAU)
    return d


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--volumes", default=cfg.DERIVED / "volumes.csv", type=Path)
    ap.add_argument("--fs-dir", type=Path, help="FreeSurfer SUBJECTS_DIR instead of a volumes CSV")
    ap.add_argument("--qc", default=cfg.DERIVED / "qc.csv", type=Path)
    a = ap.parse_args()
    cfg.OUTPUTS.mkdir(parents=True, exist_ok=True)

    subj = pd.read_csv(cfg.DERIVED / "subjects.csv", parse_dates=["ScanDate"])
    subj["in_matched_sample"] = subj.Action == "KEEP"
    subj["sex_male"] = (subj.Sex == "M").astype(int)
    assert subj.Subject.is_unique, "a participant appears more than once in subjects.csv"

    vols = volumes_from_freesurfer(a.fs_dir, subj.Subject) if a.fs_dir else pd.read_csv(a.volumes)
    assert vols.Subject.is_unique, "more than one volume row for a participant - keep one scan per person"
    for k in REGIONS:
        if f"{k}_L" in vols and f"{k}_R" in vols:
            vols[f"{k}_bilat"] = vols[f"{k}_L"] + vols[f"{k}_R"]

    csf = nearest_csf(pd.read_csv(cfg.CSF_FILE), dict(zip(subj.Subject, subj.ScanDate)))
    df = subj.merge(vols, on="Subject", how="left").merge(csf, on="Subject", how="left")

    if a.qc.exists():
        qc = pd.read_csv(a.qc)
        df = df.merge(qc[["Subject", "QC_pass"]], on="Subject", how="left")
    else:
        print(f"WARNING: no QC file at {a.qc}; QC_pass left empty (treated as not yet checked)")
        df["QC_pass"] = np.nan
    df["has_volumes"] = df[[f"{k}_bilat" for k in REGIONS if f"{k}_bilat" in df]].notna().all(axis=1) & df.eTIV.notna()
    df["has_CSF_sample"] = df.ABETA42.notna() | df.TAU.notna()
    assert df.Subject.is_unique

    df.to_csv(cfg.DERIVED / "analysis_dataset.csv", index=False)
    m = df[df.in_matched_sample]
    rep = ["# Merge report (aggregate counts)", "",
           "| Group | Matched | With volumes | QC passed | With CSF | Volumes + CSF (+QC if done) |",
           "|:---|:---:|:---:|:---:|:---:|:---:|"]
    for g in ["AD", "MCI", "CN"]:
        x = m[m.Group == g]
        qc_ok = x.QC_pass.fillna(1).astype(bool)
        rep.append(f"| {g} | {len(x)} | {int(x.has_volumes.sum())} | {int((x.QC_pass == 1).sum())} | "
                   f"{int(x.has_CSF_sample.sum())} | {int((x.has_volumes & x.has_CSF_sample & qc_ok).sum())} |")
    g = m.CSF_gap_days.dropna()
    if len(g):
        rep += ["", f"CSF-to-MRI gap (days, negative = CSF before MRI): median {g.median():.0f}, "
                    f"range {g.min():.0f} to {g.max():.0f}; within 6 months: {(g.abs() <= 183).sum()} of {len(g)}.",
                f"CSF values outside the assay range: {int(m.CSF_out_of_range.fillna(False).sum())}."]
    rep.append("")
    (cfg.OUTPUTS / "merge_report.md").write_text("\n".join(rep))
    print("\n".join(rep))


if __name__ == "__main__":
    main()
