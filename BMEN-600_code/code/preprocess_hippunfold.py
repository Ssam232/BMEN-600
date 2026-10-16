"""Step 2 helpers - preparing scans for HippUnfold and collecting its subfield volumes.

Used by 02_preprocessing.ipynb (run the notebook, not this file). Based on the team's original
BMEN600.ipynb (DICOM -> NIfTI conversion, BIDS restructuring, renaming), generalised to all three
groups, driven by the participant spreadsheet and the paths in config.py.

Workflow
  1. convert_to_bids()      DICOM (or an already-converted NIfTI) -> BIDS_DIR/sub-XXX/anat/sub-XXX_T1w.nii.gz
  2. hippunfold_command()   prints the HippUnfold command to run (outside Python)
  3. collect_volumes()      subfield volumes per hemisphere from HippUnfold's label maps
  4. add_etiv()             adds intracranial volume (eTIV) to the volumes table
  5. make_qc_template()     data/derived/qc.csv for the visual quality check
Source DICOMs are never moved or deleted; finished subjects are skipped on re-runs.
"""
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import config as cfg


# --------------------------------------------------------------------------- participants
def load_participants(path=None):
    """Final study sample -> DataFrame with Subject (ADNI ID), ID (sub-XXX), Group, Sex, Age, ScanDate, Image_ID."""
    p = pd.read_excel(path or cfg.PARTICIPANTS_FILE, dtype=str)
    p = p.rename(columns={"Image Data ID": "Image_ID", "Acq Date": "ScanDate"})
    p["Age"] = p["Age"].astype(int)
    p["ScanDate"] = pd.to_datetime(p["ScanDate"])
    for col in ("Subject", "ID", "Image_ID"):
        assert p[col].is_unique, f"duplicate {col} in participant list"
    return p[["Subject", "ID", "Group", "Sex", "Age", "ScanDate", "Image_ID"]]


# --------------------------------------------------------------------------- 1. DICOM -> BIDS
def _subject_root(row):
    for cand in (cfg.DICOM_DIR / row.Group / row.ID, cfg.DICOM_DIR / row.Group / row.Subject,
                 cfg.DICOM_DIR / row.ID, cfg.DICOM_DIR / row.Subject):
        if cand.is_dir():
            return cand
    return None


def _find_source(row):
    """Return ('nifti', file) or ('dicom', folder) for one participant, or (None, reason)."""
    root = _subject_root(row)
    if root is None:
        return None, f"no folder named {row.ID} or {row.Subject} under {cfg.DICOM_DIR}"
    niftis = sorted(set(root.rglob(f"{row.ID}*T1w.nii.gz")) | set(root.rglob(f"{row.Subject}*T1w.nii.gz")))  # converted earlier
    if niftis:
        return "nifti", niftis[0]
    dirs = sorted({f.parent for f in root.rglob("*") if f.is_file() and f.suffix.lower() == ".dcm"})
    if not dirs:
        return None, f"no .dcm files under {root}"
    if len(dirs) > 1:  # several series: prefer the one for the selected image ID
        match = [d for d in dirs if row.Image_ID in str(d)]
        if len(match) != 1:
            return None, f"{len(dirs)} DICOM series under {root}; could not pick {row.Image_ID}"
        dirs = match
    return "dicom", dirs[0]


def convert_to_bids(participants, overwrite=False):
    """Create BIDS_DIR/sub-XXX/anat/sub-XXX_T1w.nii.gz for every participant. Returns a status table."""
    import dicom2nifti  # imported here so the rest of the module works without it
    cfg.BIDS_DIR.mkdir(parents=True, exist_ok=True)
    status = []
    for row in participants.itertuples():
        target = cfg.BIDS_DIR / row.ID / "anat" / f"{row.ID}_T1w.nii.gz"
        if target.exists() and not overwrite:
            status.append((row.ID, "already done")); continue
        kind, src = _find_source(row)
        if kind is None:
            status.append((row.ID, f"MISSING: {src}")); continue
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            if kind == "nifti":
                shutil.copy2(src, target)
            else:
                with tempfile.TemporaryDirectory() as tmp:
                    dicom2nifti.convert_directory(str(src), tmp, compression=True, reorient=True)
                    out = sorted(Path(tmp).glob("*.nii.gz"))
                    if len(out) != 1:
                        raise RuntimeError(f"expected 1 NIfTI, got {len(out)}")
                    shutil.move(str(out[0]), target)
            status.append((row.ID, f"converted from {kind}"))
        except Exception as e:  # keep going; report at the end
            status.append((row.ID, f"FAILED: {e}"))
    desc = cfg.BIDS_DIR / "dataset_description.json"
    if not desc.exists():
        desc.write_text(json.dumps({"Name": "BMEN600 ADNI3 year-1 T1w", "BIDSVersion": "1.8.0"}, indent=2))
    return pd.DataFrame(status, columns=["ID", "status"])


# --------------------------------------------------------------------------- 2. HippUnfold
def hippunfold_command(participant_ids=None, docker=True):
    """The command to run HippUnfold on the BIDS folder (run it in a terminal, not here)."""
    out = cfg.HIPPUNFOLD_DIRS[0]
    labels = f" --participant-label {' '.join(i.replace('sub-', '') for i in participant_ids)}" if participant_ids else ""
    args = f"participant --modality T1w -p --cores all{labels}"
    if docker:
        return (f'docker run -it --rm -v "{cfg.BIDS_DIR}:/bids" -v "{out}:/output" '
                f"khanlab/hippunfold:latest /bids /output {args}")
    return f'hippunfold "{cfg.BIDS_DIR}" "{out}" {args}'


# --------------------------------------------------------------------------- 3. volumes
def _label_table(dseg_file):
    """Label number -> name, from HippUnfold's *_dseg.tsv next to the map if present, else config defaults."""
    for tsv in (dseg_file.with_name(dseg_file.name.replace(".nii.gz", ".tsv")),
                *sorted(dseg_file.parent.glob("*desc-subfields*dseg.tsv"))):
        if tsv.exists():
            t = pd.read_csv(tsv, sep="\t")
            idx = next((c for c in t.columns if c.lower() in ("index", "label", "value")), None)
            name = next((c for c in t.columns if c.lower() in ("abbreviation", "name", "label_name")), None)
            if idx and name:
                return {int(i): str(n) for i, n in zip(t[idx], t[name])}
    return dict(cfg.HIPPUNFOLD_LABELS)


def _find_dseg(pid, hemi):
    pat = f"{pid}_hemi-{hemi}_space-cropT1w_desc-subfields_*dseg.nii.gz"   # recommended for volumetry
    for root in cfg.HIPPUNFOLD_DIRS:
        hits = sorted(Path(root).rglob(pat)) or sorted(Path(root).rglob(pat.replace("cropT1w", "T1w")))
        if hits:
            return hits[0]
    return None


def subfield_volumes(dseg_file):
    """Volume (mm^3) of each analysis region in one hemisphere's HippUnfold label map."""
    import nibabel as nib
    img = nib.load(str(dseg_file))
    data = np.asarray(img.dataobj).astype(int)
    voxel = float(np.prod(img.header.get_zooms()[:3]))
    names = _label_table(dseg_file)
    by_name = {names[i]: (data == i).sum() * voxel for i in names}
    out = {k: by_name.get(lbl, np.nan) for k, lbl in cfg.SUBFIELDS.items()}
    for k, lbls in cfg.WHOLE.items():
        out[k] = sum(by_name.get(lbl, 0.0) for lbl in lbls)
    missing = [lbl for lbl in cfg.SUBFIELDS.values() if lbl not in by_name]
    if missing:
        raise ValueError(f"{dseg_file.name}: labels {missing} not found; check config.HIPPUNFOLD_LABELS")
    return out


def collect_volumes(participants):
    """data/derived/volumes.csv: Subject, ID, <region>_L, <region>_R (eTIV added by add_etiv)."""
    rows, missing = [], []
    for row in participants.itertuples():
        r = {"Subject": row.Subject, "ID": row.ID}
        for hemi in ("L", "R"):
            f = _find_dseg(row.ID, hemi)
            if f is None:
                missing.append(f"{row.ID} hemi-{hemi}"); continue
            r.update({f"{k}_{hemi}": v for k, v in subfield_volumes(f).items()})
        rows.append(r)
    vols = pd.DataFrame(rows)
    cfg.DERIVED.mkdir(parents=True, exist_ok=True)
    old = cfg.DERIVED / "volumes.csv"
    if old.exists() and "eTIV" in pd.read_csv(old).columns:   # keep eTIV already added
        vols = vols.merge(pd.read_csv(old)[["Subject", "eTIV"]], on="Subject", how="left")
    vols.to_csv(old, index=False)
    return vols, missing


# --------------------------------------------------------------------------- 4. eTIV
def add_etiv(participants, source=None):
    """Add eTIV to volumes.csv.

    source: CSV with Subject,eTIV  (default config.ETIV_FILE), or an ADNI UCSF FreeSurfer table
    (e.g. UCSFFSX7*.csv) with PTID/IMAGEUID/ST10CV columns - matched on the image ID of the scan used.
    """
    src = Path(source or cfg.ETIV_FILE)
    e = pd.read_csv(src)
    if {"Subject", "eTIV"} <= set(e.columns):
        e = e[["Subject", "eTIV"]]
    elif {"PTID", "IMAGEUID", "ST10CV"} <= set(e.columns):
        img = participants.assign(IMAGEUID=participants.Image_ID.str.lstrip("I").astype(int))
        e = e.assign(IMAGEUID=pd.to_numeric(e.IMAGEUID, errors="coerce"))
        e = img.merge(e[["PTID", "IMAGEUID", "ST10CV"]], left_on=["Subject", "IMAGEUID"],
                      right_on=["PTID", "IMAGEUID"], how="inner")[["Subject", "ST10CV"]].rename(columns={"ST10CV": "eTIV"})
    else:
        raise ValueError(f"{src.name}: need columns Subject,eTIV or PTID,IMAGEUID,ST10CV")
    assert e.Subject.is_unique, "more than one eTIV per participant"
    path = cfg.DERIVED / "volumes.csv"
    vols = pd.read_csv(path).drop(columns="eTIV", errors="ignore").merge(e, on="Subject", how="left")
    vols.to_csv(path, index=False)
    return vols


# --------------------------------------------------------------------------- 5. QC
def make_qc_template(participants):
    """data/derived/qc.csv with an empty QC_pass column (1 = pass, 0 = fail). Never overwrites."""
    path = cfg.DERIVED / "qc.csv"
    if path.exists():
        return path, False
    cfg.DERIVED.mkdir(parents=True, exist_ok=True)
    participants[["Subject", "ID", "Group"]].assign(QC_pass="", QC_notes="").to_csv(path, index=False)
    return path, True
