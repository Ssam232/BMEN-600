"""Step 1 - participant selection and age/sex matching.

Inputs (data/raw/, see data/README.md):
  * IDA Advanced Image Search exports for AD, MCI and CN (ADNI3, year-1 visit, MPRAGE)
  * UPenn Roche Elecsys CSF file (used only to know who has CSF data)

What it does:
  1. Collapses the image-level IDA exports to one row per participant
     (age = age at the earliest year-1 scan; all image IDs kept).
  2. Matches each AD participant, then each MCI participant, to CONTROLS_PER_CASE
     cognitively normal (CN) controls of the same sex and within +/- AGE_CALIPER years.
     Matching is solved exactly as an integer linear program (scipy.optimize.milp, HiGHS):
       - first priority: keep as many cases as possible fully matched,
       - second: prefer keeping cases that have CSF data,
       - then: minimise the total age difference, preferring controls with CSF data.
     No control is used twice within a comparison; a control may serve one AD and one MCI
     participant (separate comparisons). For MCI, controls already used for AD are
     mildly preferred so the total control sample stays small.
  3. Optionally keeps listed cases with a single control (config.ALLOW_ONE_CONTROL).
  4. Verifies every pair (same sex, age difference within the caliper, no reuse).
  5. Writes subject-level keep/remove lists (data/derived/, git-ignored) and the
     aggregate demographics table (outputs/, safe to share).

Run:  python code/01_select_and_match.py
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import Bounds, LinearConstraint, milp

import config as cfg


# --------------------------------------------------------------------------- loading
def load_subjects(path, group, csf_ids):
    """Image-level IDA export -> one row per participant."""
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    d = d.assign(Age=d["Age"].astype(int), Acq=pd.to_datetime(d["Acq Date"])).sort_values("Acq")
    s = (d.groupby("Subject")
           .agg(Sex=("Sex", "first"), Age=("Age", "first"), ScanDate=("Acq", "first"),
                Image_IDs=("Image Data ID", "; ".join), N_images=("Image Data ID", "count"))
           .reset_index())
    s.insert(1, "Group", group)
    s["Has_CSF"] = s["Subject"].isin(csf_ids)
    return s


# --------------------------------------------------------------------------- matching
def eligible(case, ctrl):
    return case.Sex == ctrl.Sex and abs(case.Age - ctrl.Age) <= cfg.AGE_CALIPER


def optimal_match(cases, ctrls, extra_cost=None):
    """Return list of (case_index, control_index) pairs from the exact ILP."""
    E = [(a, c) for a in cases.index for c in ctrls.index if eligible(cases.loc[a], ctrls.loc[c])]
    nx, ny = len(E), len(cases)
    cost = np.zeros(nx + ny)
    for k, (a, c) in enumerate(E):
        cost[k] = (abs(cases.Age[a] - ctrls.Age[c]) + (0 if ctrls.Has_CSF[c] else 0.5)
                   + (extra_cost(c) if extra_cost else 0))
    for a in cases.index:  # large reward for each fully matched case, extra if it has CSF
        cost[nx + a] = -(1000 + (50 if cases.Has_CSF[a] else 0))
    rows, lb, ub = [], [], []
    for a in cases.index:  # a kept case gets exactly K controls; a removed case gets none
        r = np.zeros(nx + ny); r[[k for k, (aa, _) in enumerate(E) if aa == a]] = 1; r[nx + a] = -cfg.CONTROLS_PER_CASE
        rows.append(r); lb.append(0); ub.append(0)
    for c in ctrls.index:  # each control used at most once in this comparison
        r = np.zeros(nx + ny); r[[k for k, (_, cc) in enumerate(E) if cc == c]] = 1
        rows.append(r); lb.append(0); ub.append(1)
    res = milp(cost, constraints=LinearConstraint(np.array(rows), lb, ub),
               integrality=np.ones(nx + ny), bounds=Bounds(0, 1))
    if not res.success:
        raise RuntimeError(f"Matching failed: {res.message}")
    x = np.round(res.x).astype(int)
    return [E[k] for k in range(nx) if x[k]]


def add_single_control_cases(cases, ctrls, pairs):
    """Keep config.ALLOW_ONE_CONTROL cases with their closest unused eligible control."""
    used = {c for _, c in pairs}
    matched = {a for a, _ in pairs}
    for sid in cfg.ALLOW_ONE_CONTROL:
        idx = cases.index[cases.Subject == sid]
        if len(idx) == 0 or idx[0] in matched:
            continue
        a = idx[0]
        free = [c for c in ctrls.index if c not in used and eligible(cases.loc[a], ctrls.loc[c])]
        if free:
            c = min(free, key=lambda c: (abs(ctrls.Age[c] - cases.Age[a]), not ctrls.Has_CSF[c], ctrls.Subject[c]))
            pairs.append((a, c)); used.add(c)
    return pairs


def verify(cases, ctrls, pairs):
    ctrl_ids = [c for _, c in pairs]
    assert len(ctrl_ids) == len(set(ctrl_ids)), "a control was reused within a comparison"
    for a, c in pairs:
        assert cases.Sex[a] == ctrls.Sex[c], f"sex mismatch {cases.Subject[a]}"
        assert abs(cases.Age[a] - ctrls.Age[c]) <= cfg.AGE_CALIPER, f"age caliper broken {cases.Subject[a]}"
    n_ctrl = pd.Series([a for a, _ in pairs]).value_counts()
    for a, n in n_ctrl.items():
        assert n == cfg.CONTROLS_PER_CASE or cases.Subject[a] in cfg.ALLOW_ONE_CONTROL


# --------------------------------------------------------------------------- explanations
SEXW = {"M": "male", "F": "female"}


def explain_cases(cases, ctrls, pairs, label):
    used = {c for _, c in pairs}
    out = {}
    for a in cases.index:
        mine = sorted([c for aa, c in pairs if aa == a], key=lambda c: abs(ctrls.Age[c] - cases.Age[a]))
        el = [c for c in ctrls.index if eligible(cases.loc[a], ctrls.loc[c])]
        lo, hi, sx = cases.Age[a] - cfg.AGE_CALIPER, cases.Age[a] + cfg.AGE_CALIPER, SEXW[cases.Sex[a]]
        if len(mine) == cfg.CONTROLS_PER_CASE:
            out[a] = ("KEEP", f"Included: {len(mine)} {sx} controls within +/-{cfg.AGE_CALIPER} y ("
                      + ", ".join(f"{ctrls.Subject[c]} age {ctrls.Age[c]}" for c in mine) + ")")
        elif mine:
            out[a] = ("KEEP", f"Included with 1 control ({ctrls.Subject[mine[0]]}, age {ctrls.Age[mine[0]]}): "
                      f"the only {sx} control aged {lo}-{hi}; kept by team decision")
        elif not el:
            out[a] = ("REMOVE", f"Removed: no {sx} controls aged {lo}-{hi}")
        elif len(el) < cfg.CONTROLS_PER_CASE:
            out[a] = ("REMOVE", f"Removed: only {len(el)} {sx} control aged {lo}-{hi}; {cfg.CONTROLS_PER_CASE} needed")
        else:
            free = [c for c in el if c not in used]
            out[a] = ("REMOVE", f"Removed: {len(el)} eligible {sx} controls aged {lo}-{hi}, but "
                      f"{len(el) - len(free)} were already used by other {label} participants, "
                      f"leaving {len(free)}; {cfg.CONTROLS_PER_CASE} needed")
    return out


# --------------------------------------------------------------------------- demographics
def demographics(ad, ad_ctrl, mci, mci_ctrl, cn_all):
    allp = pd.concat([ad.assign(G="AD"), mci.assign(G="MCI"), cn_all.assign(G="CN")], ignore_index=True)
    same = set(mci_ctrl.Subject) == set(cn_all.Subject)
    groups = [("AD", ad), ("Controls matched to AD", ad_ctrl), ("MCI", mci)]
    if not same:
        groups.append(("Controls matched to MCI", mci_ctrl))
    groups += [("All controls", cn_all), ("Total sample", allp)]
    pct = lambda k, n: f"{k} ({100 * k / n:.1f})"
    rows = []
    def add(label, fn): rows.append([label] + [fn(g) for _, g in groups])
    add("Participants, n", lambda g: str(len(g)))
    add("MRI scans (images), n", lambda g: str(int(g.N_images.sum())))
    add("Sex: Male, n (%)", lambda g: pct((g.Sex == "M").sum(), len(g)))
    add("Sex: Female, n (%)", lambda g: pct((g.Sex == "F").sum(), len(g)))
    add("Age (y), mean (SD)", lambda g: f"{g.Age.mean():.1f} ({g.Age.std():.1f})")
    add("Age (y), median [IQR]", lambda g: f"{g.Age.median():.0f} [{g.Age.quantile(.25):.0f}-{g.Age.quantile(.75):.0f}]")
    add("Age (y), range", lambda g: f"{g.Age.min()}-{g.Age.max()}")
    for lo, hi, lab in [(0, 69, "under 70"), (70, 79, "70-79"), (80, 200, "80+")]:
        add(f"Age group {lab}, n (%)", lambda g, lo=lo, hi=hi: pct(g.Age.between(lo, hi).sum(), len(g)))
    add("Male age (y), mean (SD)", lambda g: f"{g[g.Sex == 'M'].Age.mean():.1f} ({g[g.Sex == 'M'].Age.std():.1f})")
    add("Female age (y), mean (SD)", lambda g: f"{g[g.Sex == 'F'].Age.mean():.1f} ({g[g.Sex == 'F'].Age.std():.1f})")
    add("CSF (Roche Elecsys) available, n (%)", lambda g: pct(int(g.Has_CSF.sum()), len(g)))
    T = pd.DataFrame(rows, columns=["Characteristic"] + [f"{n} (n={len(g)})" for n, g in groups])

    def fisher(a, b, col, val):
        t = [[(a[col] == val).sum(), (a[col] != val).sum()], [(b[col] == val).sum(), (b[col] != val).sum()]]
        return stats.fisher_exact(t)[1]
    P = []
    for lab, a, b in [("AD vs matched controls", ad, ad_ctrl), ("MCI vs matched controls", mci, mci_ctrl)]:
        P.append([lab, stats.ttest_ind(a.Age, b.Age, equal_var=False).pvalue,
                  fisher(a, b, "Sex", "M"), fisher(a, b, "Has_CSF", True)])
    P.append(["AD vs MCI vs all controls", stats.f_oneway(ad.Age, mci.Age, cn_all.Age).pvalue,
              stats.chi2_contingency(pd.crosstab(allp.G, allp.Sex))[1],
              stats.chi2_contingency(pd.crosstab(allp.G, allp.Has_CSF))[1]])
    P = pd.DataFrame(P, columns=["Comparison", "Age p", "Sex p", "CSF availability p"])
    return T, P, same


def table1_markdown(T, P, n_ad, n_mci, n_ad_ctrl, n_cn, same):
    fp = lambda v: "< 0.001" if v < 0.001 else f"{v:.3f}"
    md = ["## Table 1. Demographics of the matched sample", "",
          f"ADNI3, year-1 MPRAGE scans. Matching: each AD and MCI participant has {cfg.CONTROLS_PER_CASE} controls "
          f"of the same sex, aged within +/-{cfg.AGE_CALIPER} years"
          + (f" (one AD participant, {', '.join(cfg.ALLOW_ONE_CONTROL)}, has 1 control because only one exists)."
             if cfg.ALLOW_ONE_CONTROL else "."), "",
          "| " + " | ".join(T.columns) + " |", "|" + "|".join([":---"] + [":---:"] * (len(T.columns) - 1)) + "|"]
    md += ["| " + " | ".join(map(str, r)) + " |" for r in T.itertuples(index=False)]
    md += ["", "### Group comparisons (p-values)", "", "| Comparison | Age | Sex | CSF availability |",
           "|:---|:---:|:---:|:---:|"]
    md += [f"| {r[0]} | {fp(r[1])} | {fp(r[2])} | {fp(r[3])} |" for r in P.itertuples(index=False)]
    overlap = (f"All {n_cn} are matched to MCI participants, and {n_ad_ctrl} of them are also matched to AD "
               "participants, so the control columns overlap." if same else
               "Some controls are matched to both an AD and an MCI participant, so the control columns overlap.")
    md += ["", "**Notes**", "",
           "- Age is age at the year-1 MRI scan, as recorded in IDA. Percentages are of each column's n.",
           f"- \"All controls\" is every control kept. {overlap}",
           f"- Total sample counts each person once ({n_ad} AD + {n_mci} MCI + {n_cn} controls = {n_ad + n_mci + n_cn}).",
           "- \"CSF available\" means the participant has at least one Roche Elecsys CSF result "
           "(UPENNBIOMK_ROCHE_ELECSYS, any visit).",
           "- Tests: age, Welch t-test for pairwise comparisons and one-way ANOVA across the three groups; "
           "sex and CSF availability, Fisher exact test for pairwise comparisons and chi-square across the three groups.",
           "- MCI participants with CSF results were prioritised when controls ran out, so all kept MCI "
           "participants have CSF.", ""]
    return "\n".join(md)


# --------------------------------------------------------------------------- main
def main():
    cfg.DERIVED.mkdir(parents=True, exist_ok=True); cfg.OUTPUTS.mkdir(parents=True, exist_ok=True)
    csf_ids = set(pd.read_csv(cfg.CSF_FILE, usecols=["PTID"]).PTID)
    AD = load_subjects(cfg.IDA_EXPORTS["AD"], "AD", csf_ids)
    MCI = load_subjects(cfg.IDA_EXPORTS["MCI"], "MCI", csf_ids)
    CN = load_subjects(cfg.IDA_EXPORTS["CN"], "CN", csf_ids)
    assert not (set(AD.Subject) | set(MCI.Subject)) & set(CN.Subject), "a participant appears in two groups"

    ad_pairs = add_single_control_cases(AD, CN, optimal_match(AD, CN))
    verify(AD, CN, ad_pairs)
    ad_ctrl_idx = {c for _, c in ad_pairs}
    mci_pairs = optimal_match(MCI, CN, extra_cost=lambda c: 0 if c in ad_ctrl_idx else 0.25)
    verify(MCI, CN, mci_pairs)

    ad_exp, mci_exp = explain_cases(AD, CN, ad_pairs, "AD"), explain_cases(MCI, CN, mci_pairs, "MCI")
    AD["Action"], AD["Explanation"] = zip(*[ad_exp[i] for i in AD.index])
    MCI["Action"], MCI["Explanation"] = zip(*[mci_exp[i] for i in MCI.index])
    AD["Matched_controls"] = [" / ".join(CN.Subject[c] for aa, c in ad_pairs if aa == a) for a in AD.index]
    MCI["Matched_controls"] = [" / ".join(CN.Subject[c] for aa, c in mci_pairs if aa == a) for a in MCI.index]

    to_ad = {c: AD.Subject[a] for a, c in ad_pairs}
    to_mci = {c: MCI.Subject[a] for a, c in mci_pairs}
    CN["Matched_to_AD"] = [to_ad.get(c, "") for c in CN.index]
    CN["Matched_to_MCI"] = [to_mci.get(c, "") for c in CN.index]
    CN["Action"] = np.where((CN.Matched_to_AD != "") | (CN.Matched_to_MCI != ""), "KEEP", "REMOVE")
    kept_cases = pd.concat([AD[AD.Action == "KEEP"], MCI[MCI.Action == "KEEP"]])
    exps = []
    for c in CN.index:
        roles = ([f"control for AD {to_ad[c]}"] if c in to_ad else []) + \
                ([f"control for MCI {to_mci[c]}"] if c in to_mci else [])
        if roles:
            exps.append("Included: " + "; ".join(roles)); continue
        el_any = [a for g in (AD, MCI) for a in g.index if eligible(g.loc[a], CN.loc[c])]
        el_kept = kept_cases[[eligible(r, CN.loc[c]) for r in kept_cases.itertuples()]] if len(kept_cases) else []
        if not el_any:
            exps.append(f"Removed: no AD or MCI participant who is {SEXW[CN.Sex[c]]} and within +/-{cfg.AGE_CALIPER} y")
        elif len(el_kept) == 0:
            exps.append("Removed: its only age/sex matches were AD/MCI participants who could not be fully matched")
        else:
            exps.append("Removed: spare - the AD/MCI participants it could match already have closer-aged controls")
    CN["Explanation"] = exps

    subjects = pd.concat([AD, MCI, CN], ignore_index=True)
    subjects["ScanDate"] = subjects.ScanDate.dt.date
    subjects.to_csv(cfg.DERIVED / "subjects.csv", index=False)
    for label, cases, pairs in [("AD", AD, ad_pairs), ("MCI", MCI, mci_pairs)]:
        rows = []
        for a in cases.index:
            mine = sorted([c for aa, c in pairs if aa == a], key=lambda c: abs(CN.Age[c] - cases.Age[a]))
            r = {"Action": cases.Action[a], "Subject": cases.Subject[a], "Sex": cases.Sex[a], "Age": cases.Age[a],
                 "Has_CSF": cases.Has_CSF[a]}
            for i in range(cfg.CONTROLS_PER_CASE):
                c = mine[i] if i < len(mine) else None
                r[f"CN_{i+1}"] = CN.Subject[c] if c is not None else ""
                r[f"CN_{i+1}_Age"] = CN.Age[c] if c is not None else None
                r[f"CN_{i+1}_Has_CSF"] = CN.Has_CSF[c] if c is not None else None
            r["Explanation"] = cases.Explanation[a]; rows.append(r)
        pd.DataFrame(rows).to_csv(cfg.DERIVED / f"matches_{label}.csv", index=False)

    keep = lambda d: d[d.Action == "KEEP"]
    ad_k, mci_k, cn_k = keep(AD), keep(MCI), keep(CN)
    T, P, same = demographics(ad_k, CN.loc[sorted(ad_ctrl_idx)], mci_k,
                              CN.loc[sorted(c for _, c in mci_pairs)], cn_k)
    T.to_csv(cfg.OUTPUTS / "table1_demographics.csv", index=False)
    P.to_csv(cfg.OUTPUTS / "table1_group_comparisons.csv", index=False)
    (cfg.OUTPUTS / "table1_demographics.md").write_text(table1_markdown(T, P, len(ad_k), len(mci_k), len(ad_ctrl_idx), len(cn_k), same))

    summary = ["# Matching summary", "",
               f"Rules: same sex, age within +/-{cfg.AGE_CALIPER} y, {cfg.CONTROLS_PER_CASE} controls per case, "
               "no control reused within a comparison. All pairs verified automatically.", "",
               "| Group | In IDA export | Kept | Removed | Kept with CSF |", "|:---|:---:|:---:|:---:|:---:|"]
    for g in (AD, MCI, CN):
        summary.append(f"| {g.Group.iloc[0]} | {len(g)} | {(g.Action == 'KEEP').sum()} | "
                       f"{(g.Action == 'REMOVE').sum()} | {int(keep(g).Has_CSF.sum())} |")
    summary.append("")
    (cfg.OUTPUTS / "matching_summary.md").write_text("\n".join(summary))
    print("\n".join(summary))
    print(f"Wrote {cfg.DERIVED / 'subjects.csv'} and outputs/table1_demographics.md")


if __name__ == "__main__":
    main()
