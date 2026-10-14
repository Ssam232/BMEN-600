"""Step 4 - statistical analysis: CSF biomarkers vs hippocampal subfield volumes.

Input:  data/derived/analysis_dataset.csv  (from step 3; one row per participant)
Output: outputs/results/*.csv and outputs/figures/*.png  (aggregate statistics only)

Analyses
  Primary   Within each group (CN, MCI, AD): partial correlation between each CSF marker and
            each subfield volume (left + right), adjusted for age, sex and eTIV.
            95% CIs from 5,000 bootstrap resamples; Benjamini-Hochberg FDR across
            markers x subfields within each group.
  Primary   Pooled model per marker-subfield pair:
              volume ~ marker * group + age + sex + eTIV
            F-test of the marker x group interaction (does the association differ by stage?).
  Compare   Same correlations for whole-hippocampus volume (the measure most prior CSF studies used).
  Secondary Group differences in volume (ANCOVA: volume ~ group + age + sex + eTIV), matched sample.
  Explore   Marker x sex interaction in the pooled model.
  Sensitivity
            (i)   all participants with volumes + CSF, ignoring matching
            (ii)  excluding CSF values outside the assay range
            (iii) left and right hemispheres separately
            (iv)  Spearman (rank-based) partial correlations

Each participant contributes one row; controls matched to both an AD and an MCI
participant appear once. Participants with QC_pass == 0 are excluded.

Run:  python code/04_analysis.py
"""
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

import config as cfg  # (not "C": that name is reserved for categorical terms in model formulas)

GROUPS = ["CN", "MCI", "AD"]
SUB = list(cfg.SUBFIELDS)
WH = list(cfg.WHOLE)
RES = cfg.OUTPUTS / "results"
FIG = cfg.OUTPUTS / "figures"


# --------------------------------------------------------------------------- helpers
def residualize(Y, Z):
    """Residuals of each column of Y after regressing on Z (with intercept)."""
    Z1 = np.column_stack([np.ones(len(Z)), Z])
    beta, *_ = np.linalg.lstsq(Z1, Y, rcond=None)
    return Y - Z1 @ beta


def partial_corr_matrix(X, Y, Z):
    """Partial correlations between every column of X and every column of Y, given Z."""
    rx, ry = residualize(X, Z), residualize(Y, Z)
    rx = (rx - rx.mean(0)) / rx.std(0)
    ry = (ry - ry.mean(0)) / ry.std(0)
    return rx.T @ ry / len(rx)


def partial_p(r, n, k):
    dof = n - 2 - k
    t = r * np.sqrt(dof / np.clip(1 - r ** 2, 1e-12, None))
    return 2 * stats.t.sf(np.abs(t), dof)


def fdr(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def correlations(df, cols, label, rank=False, boot=True):
    """Within-group partial correlations for every marker x column in `cols`."""
    rows = []
    rng = np.random.default_rng(cfg.SEED)
    for g in GROUPS:
        d = df[df.Group == g].dropna(subset=list(cfg.MARKERS) + cols + cfg.COVARIATES)
        n = len(d)
        if n < cfg.MIN_N:
            print(f"  [{label}] {g}: n={n} < {cfg.MIN_N}, skipped")
            continue
        X, Y, Z = (d[list(cfg.MARKERS)].to_numpy(float), d[cols].to_numpy(float), d[cfg.COVARIATES].to_numpy(float))
        if rank:
            X, Y, Z = (np.apply_along_axis(stats.rankdata, 0, A) for A in (X, Y, Z))
        R = partial_corr_matrix(X, Y, Z)
        P = partial_p(R, n, Z.shape[1])
        lo = hi = np.full_like(R, np.nan)
        if boot:
            B = np.empty((cfg.N_BOOT,) + R.shape)
            for b in range(cfg.N_BOOT):
                i = rng.integers(0, n, n)
                with np.errstate(invalid="ignore", divide="ignore"):
                    B[b] = partial_corr_matrix(X[i], Y[i], Z[i])
            lo, hi = np.nanpercentile(B, 2.5, axis=0), np.nanpercentile(B, 97.5, axis=0)
        for i, m in enumerate(cfg.MARKERS):
            for j, v in enumerate(cols):
                rows.append({"analysis": label, "group": g, "marker": m, "region": v, "n": n,
                             "partial_r": R[i, j], "ci_low": lo[i, j], "ci_high": hi[i, j], "p": P[i, j]})
    out = pd.DataFrame(rows)
    if len(out):
        out["q_fdr"] = out.groupby("group")["p"].transform(fdr)
        out["significant_fdr"] = out.q_fdr < cfg.FDR_Q
    return out


def interaction_tests(df, cols, term, extra="", label=""):
    """F-test of marker x `term` interaction in pooled OLS models."""
    rows = []
    for m in cfg.MARKERS:
        for v in cols:
            d = df.dropna(subset=[m, v] + cfg.COVARIATES)
            if d.Group.nunique() < 2:
                continue
            base = f"{v} ~ {m} + {term} + Age + sex_male + eTIV{extra}"
            full = f"{v} ~ {m} * {term} + Age + sex_male + eTIV{extra}"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                red, fit = smf.ols(base, d).fit(), smf.ols(full, d).fit()
                a = anova_lm(red, fit)
            rows.append({"analysis": label, "marker": m, "region": v, "n": len(d),
                         "F": a["F"].iloc[1], "df_num": a["df_diff"].iloc[1], "df_den": fit.df_resid,
                         "p": a["Pr(>F)"].iloc[1]})
    out = pd.DataFrame(rows)
    if len(out):
        out["q_fdr"] = fdr(out.p)
        out["significant_fdr"] = out.q_fdr < cfg.FDR_Q
    return out


def ancova(df, cols):
    rows = []
    for v in cols:
        d = df.dropna(subset=[v] + cfg.COVARIATES)
        fit = smf.ols(f"{v} ~ C(Group, Treatment('CN')) + Age + sex_male + eTIV", d).fit()
        a = anova_lm(fit, typ=2)
        row = {"region": v, "n": len(d), "F_group": a.loc["C(Group, Treatment('CN'))", "F"],
               "p_group": a.loc["C(Group, Treatment('CN'))", "PR(>F)"]}
        for g in ("MCI", "AD"):
            t = f"C(Group, Treatment('CN'))[T.{g}]"
            if t in fit.params:
                ci = fit.conf_int().loc[t]
                row.update({f"{g}_vs_CN_diff_mm3": fit.params[t], f"{g}_ci_low": ci[0], f"{g}_ci_high": ci[1]})
        rows.append(row)
    out = pd.DataFrame(rows)
    out["q_fdr"] = fdr(out.p_group)
    return out


# --------------------------------------------------------------------------- figure
# Diverging scale: blue (negative) - neutral gray - red (positive).
DIVERGING = LinearSegmentedColormap.from_list(
    "bluegrayred", ["#184f95", "#6da7ec", "#f0efec", "#f09a99", "#a8302f"])


def heatmap(res, path):
    fig, axes = plt.subplots(1, len(GROUPS), figsize=(4.2 * len(GROUPS), 3.4), sharey=True)
    for ax, g in zip(np.atleast_1d(axes), GROUPS):
        d = res[res.group == g]
        ax.set_title(f"{g} (n={int(d.n.iloc[0])})" if len(d) else f"{g} (not enough data)", fontsize=10,
                     color="#0b0b0b")
        if not len(d):
            ax.axis("off"); continue
        M = d.pivot(index="marker", columns="region", values="partial_r").loc[list(cfg.MARKERS.values()), SUB]
        Q = d.pivot(index="marker", columns="region", values="q_fdr").loc[list(cfg.MARKERS.values()), SUB]
        im = ax.imshow(M.to_numpy(), cmap=DIVERGING, vmin=-0.8, vmax=0.8, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                r, q = M.iat[i, j], Q.iat[i, j]
                ax.text(j, i, f"{r:.2f}" + ("*" if q < cfg.FDR_Q else ""), ha="center", va="center", fontsize=8,
                        color="white" if abs(r) > 0.55 else "#0b0b0b", fontweight="bold" if q < cfg.FDR_Q else None)
        ax.set_xticks(range(len(SUB)), SUB, rotation=40, ha="right", fontsize=8, color="#52514e")
        ax.set_yticks(range(len(cfg.MARKERS)), list(cfg.MARKERS.values()), fontsize=8, color="#52514e")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("Partial r (age, sex, eTIV adjusted)", fontsize=8, color="#52514e")
    cb.outline.set_visible(False)
    fig.text(0.01, -0.06, f"* FDR q < {cfg.FDR_Q} within group", fontsize=7, color="#52514e")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- main
def main():
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(cfg.DERIVED / "analysis_dataset.csv")
    assert df.Subject.is_unique, "one row per participant expected"
    if df.QC_pass.isna().all():
        print("WARNING: no segmentation QC recorded yet - all participants with volumes are included")
    usable = df.has_volumes & df.has_CSF_sample & (df.QC_pass != 0)
    prim = df[df.in_matched_sample & usable]
    bil = [f"{k}_bilat" for k in SUB]

    print("Primary correlations ...")
    corr = [correlations(prim, bil, "primary")]
    corr.append(correlations(prim, [f"{k}_bilat" for k in WH], "whole_hippocampus"))
    print("Sensitivity analyses ...")
    corr.append(correlations(df[usable], bil, "sens_all_with_csf"))
    corr.append(correlations(prim[~prim.CSF_out_of_range.fillna(False).astype(bool)], bil, "sens_in_range_csf"))
    corr.append(correlations(prim, [f"{k}_L" for k in SUB], "sens_left", boot=False))
    corr.append(correlations(prim, [f"{k}_R" for k in SUB], "sens_right", boot=False))
    corr.append(correlations(prim, bil, "sens_spearman", rank=True, boot=False))
    corr = pd.concat(corr, ignore_index=True)
    for col in ("marker", "region"):
        corr[col] = corr[col].replace({**cfg.MARKERS, **{f"{k}_bilat": k for k in SUB + WH}})
    corr.to_csv(RES / "partial_correlations.csv", index=False)

    print("Interaction models ...")
    interaction_tests(prim, bil + [f"{k}_bilat" for k in WH], "C(Group, Treatment('CN'))",
                      label="marker_x_group").to_csv(RES / "interaction_marker_x_group.csv", index=False)
    interaction_tests(prim, bil, "sex_male", extra=" + C(Group, Treatment('CN'))",
                      label="marker_x_sex").to_csv(RES / "interaction_marker_x_sex_exploratory.csv", index=False)
    ancova(df[df.in_matched_sample & df.has_volumes & (df.QC_pass != 0)],
           bil + [f"{k}_bilat" for k in WH]).to_csv(RES / "ancova_group_differences.csv", index=False)

    heatmap(corr[corr.analysis == "primary"], FIG / "partial_correlation_heatmap.png")
    n = prim.Group.value_counts().reindex(GROUPS).fillna(0).astype(int)
    print("Done. Primary sample:", dict(n), "-> outputs/results/, outputs/figures/")


if __name__ == "__main__":
    main()
