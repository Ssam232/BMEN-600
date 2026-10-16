"""Step 3 - merge subfield volumes, eTIV, QC and CSF biomarkers into one analysis table, and write Table 1.

Runs after step 2 (02_preprocessing.ipynb). Can also run before step 2 is finished: volumes are then left
empty and only Table 1 and the CSF summary are produced.

Inputs
  * Study sample: config.PARTICIPANTS_FILE (ADNI_Participants.xlsx, one row per participant/scan).
    If that file is absent, the KEEP list from step 1 (data/derived/subjects.csv) is used instead.
  * data/derived/volumes.csv from step 2: Subject, eTIV, <region>_L, <region>_R
    (regions = keys of config.SUBFIELDS and config.WHOLE)
  * data/derived/qc.csv from step 2 (optional): Subject, QC_pass (1 = pass, 0 = fail, blank = not checked)
  * config.CSF_FILE (UPenn Roche Elecsys)

CSF rules
  * one CSF sample per participant: the exam closest in time to the MRI scan (gap kept in CSF_gap_days)
  * samples marked "Sample Hemolyzed" are not used
  * values outside the Roche Elecsys measuring range are flagged (CSF_out_of_range); below-range t-tau/p-tau
    (missing in the file) are set to the lower limit; above-range Abeta42 keeps UPenn's extrapolated value
  * t-tau and p-tau181 are log-transformed (log_TAU, log_PTAU)

Outputs
  data/derived/analysis_dataset.csv  one row per participant (subject-level, git-ignored)
  outputs/table1_demographics.md/.csv, outputs/merge_report.md  (aggregate, safe to commit)

Run:  python code/03_merge.py
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import config as cfg

REGIONS = list(cfg.SUBFIELDS) + list(cfg.WHOLE)
GROUPS = ["AD", "MCI", "CN"]


# --------------------------------------------------------------------------- sample
def load_sample():
    if cfg.PARTICIPANTS_FILE.exists():
        from preprocess_hippunfold import load_participants
        s = load_participants()
        src = cfg.PARTICIPANTS_FILE.name
    else:
        s = pd.read_csv(cfg.DERIVED / "subjects.csv", parse_dates=["ScanDate"])
        s = s[s.Action == "KEEP"]
        src = "step 1 KEEP list (data/derived/subjects.csv)"
    assert s.Subject.is_unique, "a participant appears more than once in the study sample"
    s = s.assign(in_sample=True, sex_male=(s.Sex == "M").astype(int))
    print(f"Study sample from {src}: {s.Group.value_counts().reindex(GROUPS).to_dict()}")
    return s


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


# --------------------------------------------------------------------------- Table 1
def table1(df):
    """Demographics by group (AD, MCI, CN, Total) + group comparisons, as GitHub markdown."""
    g = {k: df[df.Group == k] for k in GROUPS}
    cols = [(k, g[k]) for k in GROUPS] + [("Total", df)]
    pct = lambda k, n: f"{k} ({100 * k / n:.1f})"
    rows = []
    def add(label, fn): rows.append([label] + [fn(x) for _, x in cols])
    add("Participants, n", lambda x: str(len(x)))
    add("Sex: Male, n (%)", lambda x: pct((x.Sex == "M").sum(), len(x)))
    add("Sex: Female, n (%)", lambda x: pct((x.Sex == "F").sum(), len(x)))
    add("Age (y), mean (SD)", lambda x: f"{x.Age.mean():.1f} ({x.Age.std():.1f})")
    add("Age (y), median [IQR]", lambda x: f"{x.Age.median():.0f} [{x.Age.quantile(.25):.0f}-{x.Age.quantile(.75):.0f}]")
    add("Age (y), range", lambda x: f"{x.Age.min()}-{x.Age.max()}")
    for lo, hi, lab in [(0, 69, "under 70"), (70, 79, "70-79"), (80, 200, "80+")]:
        add(f"Age group {lab}, n (%)", lambda x, lo=lo, hi=hi: pct(x.Age.between(lo, hi).sum(), len(x)))
    add("Male age (y), mean (SD)", lambda x: f"{x[x.Sex == 'M'].Age.mean():.1f} ({x[x.Sex == 'M'].Age.std():.1f})")
    add("Female age (y), mean (SD)", lambda x: f"{x[x.Sex == 'F'].Age.mean():.1f} ({x[x.Sex == 'F'].Age.std():.1f})")
    add("CSF (Roche Elecsys) available, n (%)", lambda x: pct(int(x.Has_CSF.sum()), len(x)))
    T = pd.DataFrame(rows, columns=["Characteristic"] + [f"{k} (n={len(x)})" for k, x in cols])

    fis = lambda a, b, c, v: stats.fisher_exact([[(a[c] == v).sum(), (a[c] != v).sum()],
                                                 [(b[c] == v).sum(), (b[c] != v).sum()]])[1]
    P = [[f"{k} vs CN", stats.ttest_ind(g[k].Age, g["CN"].Age, equal_var=False).pvalue,
          fis(g[k], g["CN"], "Sex", "M"), fis(g[k], g["CN"], "Has_CSF", True)] for k in ("AD", "MCI")]
    P.append(["AD vs MCI vs CN", stats.f_oneway(*(g[k].Age for k in GROUPS)).pvalue,
              stats.chi2_contingency(pd.crosstab(df.Group, df.Sex))[1],
              stats.chi2_contingency(pd.crosstab(df.Group, df.Has_CSF))[1]])
    P = pd.DataFrame(P, columns=["Comparison", "Age p", "Sex p", "CSF availability p"])

    fp = lambda v: "< 0.001" if v < 0.001 else f"{v:.3f}"
    md = ["## Table 1. Demographics of the study sample", "",
          f"ADNI3, year-1 T1-weighted MPRAGE scans; one scan per participant ({len(df)} participants, {len(df)} scans). "
          f"Source: {cfg.PARTICIPANTS_FILE.name}.", "",
          "| " + " | ".join(T.columns) + " |", "|" + "|".join([":---"] + [":---:"] * (len(T.columns) - 1)) + "|"]
    md += ["| " + " | ".join(r) + " |" for r in T.itertuples(index=False)]
    md += ["", "### Group comparisons (p-values)", "", "| Comparison | Age | Sex | CSF availability |", "|:---|:---:|:---:|:---:|"]
    md += [f"| {r[0]} | {fp(r[1])} | {fp(r[2])} | {fp(r[3])} |" for r in P.itertuples(index=False)]
    md += ["", "**Notes**", "",
           "- Age is age at the year-1 MRI scan, as recorded in IDA. Percentages are of each column's n.",
           f"- Each AD participant was matched to {cfg.CONTROLS_PER_CASE} CN participants of the same sex, aged within "
           f"+/-{cfg.AGE_CALIPER} years (AD {', '.join(cfg.ALLOW_ONE_CONTROL)} has 1 control because only one exists).",
           "- \"CSF available\" means the participant has at least one Roche Elecsys CSF result (UPENNBIOMK_ROCHE_ELECSYS, any visit).",
           "- Tests: age, Welch t-test (pairwise) and one-way ANOVA (three groups); sex and CSF availability, "
           "Fisher exact test (pairwise) and chi-square (three groups).", ""]
    return T, P, "\n".join(md)


# --------------------------------------------------------------------------- main
def main(argv=None):
    """argv: list of command-line style options (handy from a notebook)."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--volumes", default=cfg.DERIVED / "volumes.csv", type=Path)
    ap.add_argument("--qc", default=cfg.DERIVED / "qc.csv", type=Path)
    a = ap.parse_args(argv)
    cfg.OUTPUTS.mkdir(parents=True, exist_ok=True); cfg.DERIVED.mkdir(parents=True, exist_ok=True)

    subj = load_sample()
    csf_raw = pd.read_csv(cfg.CSF_FILE)
    subj["Has_CSF"] = subj.Subject.isin(set(csf_raw.PTID))

    T, P, md = table1(subj)
    T.to_csv(cfg.OUTPUTS / "table1_demographics.csv", index=False)
    P.to_csv(cfg.OUTPUTS / "table1_group_comparisons.csv", index=False)
    (cfg.OUTPUTS / "table1_demographics.md").write_text(md)

    if a.volumes.exists():
        vols = pd.read_csv(a.volumes).drop(columns=["ID"], errors="ignore")
        assert vols.Subject.is_unique, "more than one volume row for a participant - keep one scan per person"
    else:
        print(f"NOTE: no volumes yet ({a.volumes}); run step 2 first. Writing CSF/demographics only.")
        vols = pd.DataFrame({"Subject": subj.Subject})
    for k in REGIONS:
        for h in ("L", "R"):
            if f"{k}_{h}" not in vols:
                vols[f"{k}_{h}"] = np.nan
        vols[f"{k}_bilat"] = vols[f"{k}_L"] + vols[f"{k}_R"]
    if "eTIV" not in vols:
        vols["eTIV"] = np.nan

    csf = nearest_csf(csf_raw, dict(zip(subj.Subject, subj.ScanDate)))
    df = subj.merge(vols, on="Subject", how="left").merge(csf, on="Subject", how="left")

    if a.qc.exists():
        qc = pd.read_csv(a.qc)
        df = df.merge(qc[["Subject", "QC_pass"]].assign(QC_pass=lambda q: pd.to_numeric(q.QC_pass, errors="coerce")),
                      on="Subject", how="left")
    else:
        df["QC_pass"] = np.nan
    if df.QC_pass.isna().all():
        print("WARNING: no segmentation QC recorded yet (QC_pass empty)")
    df["has_volumes"] = df[[f"{k}_bilat" for k in REGIONS]].notna().all(axis=1) & df.eTIV.notna()
    df["has_CSF_sample"] = df.ABETA42.notna() | df.TAU.notna()
    assert df.Subject.is_unique
    df.to_csv(cfg.DERIVED / "analysis_dataset.csv", index=False)

    rep = ["# Merge report (aggregate counts)", "",
           "| Group | In sample | With volumes + eTIV | QC passed | With CSF | Usable (volumes + CSF, not QC-failed) |",
           "|:---|:---:|:---:|:---:|:---:|:---:|"]
    for g in GROUPS:
        x = df[df.Group == g]
        rep.append(f"| {g} | {len(x)} | {int(x.has_volumes.sum())} | {int((x.QC_pass == 1).sum())} | "
                   f"{int(x.has_CSF_sample.sum())} | {int((x.has_volumes & x.has_CSF_sample & (x.QC_pass != 0)).sum())} |")
    gap = df.CSF_gap_days.dropna()
    if len(gap):
        rep += ["", f"CSF-to-MRI gap (days, negative = CSF before MRI): median {gap.median():.0f}, "
                    f"range {gap.min():.0f} to {gap.max():.0f}; within 6 months: {(gap.abs() <= 183).sum()} of {len(gap)}.",
                f"CSF values outside the assay range: {int(df.CSF_out_of_range.fillna(False).astype(bool).sum())}."]
    rep.append("")
    (cfg.OUTPUTS / "merge_report.md").write_text("\n".join(rep))
    print("\n".join(rep))


if __name__ == "__main__":
    main()
