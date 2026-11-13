"""Figure 4c–d reproducible analysis (derived from SimTK Stage4 B021–B023).

Requires the corrected specimen scores and final zonal feature master. No result
is hard-coded. Validates 89 distinct Medial PLM-complete patients.

Figure 4c uses PLM_SZ >= 3; Figure 4d uses full-cohort RF importance for
PLM_Total >= 5. An additional three-model OOF analysis for PLM_Total >= 5
is exported as supplementary numerical source data, without changing either figure.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve

TRAD_COLS = ["True_Total_OARSI_Mean", "True_Total_HHGS_Mean",
             "True_Structure_Mean", "True_SafO_Mean", "True_Cells_Mean"]
BRIGHT_COLS = [
 "Cartilage_Area_mm2", "Median_Thickness_um", "Surface_Roughness",
 "SZ_HE_Basophilia", "MZ_HE_Basophilia", "DZ_HE_Basophilia", "Grad_HE_Baso_DZ_SZ",
 "SZ_HE_Eosinophilia", "DZ_HE_Eosinophilia",
 "SZ_SafO_PG", "MZ_SafO_PG", "DZ_SafO_PG", "Grad_SafO_PG_DZ_SZ",
 "SZ_ChondronDensity_mm2", "MZ_ChondronDensity_mm2", "DZ_ChondronDensity_mm2",
 "SZ_ChondronArea_um2", "MZ_ChondronArea_um2", "DZ_ChondronArea_um2",
 "Grad_ChondronArea_DZ_SZ", "SZ_AspectRatio", "DZ_AspectRatio",
 "SZ_SurfaceAngle_deg", "DZ_SurfaceAngle_deg", "Grad_Angle_DZ_SZ",
 "SZ_CloningRate", "DZ_CloningRate"]

def load_cohort(scores_csv: Path, features_csv: Path) -> pd.DataFrame:
    scores = pd.read_csv(scores_csv)
    zonal = pd.read_csv(features_csv)
    keys = ["Patient", "Region", "SpecimenID"]
    for label, table in (("scores", scores), ("zonal", zonal)):
        if table.duplicated(keys).any():
            raise ValueError(f"Duplicate specimen keys in {label}")
    df = scores.merge(zonal, on=keys, how="inner", suffixes=("_score", ""), validate="one_to_one")
    for c in scores.columns:
        cs = c + "_score"
        if c not in keys and c in df.columns and cs in df.columns:
            if not (df[cs].eq(df[c]) | (df[cs].isna() & df[c].isna())).all():
                raise ValueError(f"Overlapping field mismatch: {c}")
    sub = df[(df["Region"] == "Medial") & df["PLM_Total"].notna()].copy().reset_index(drop=True)
    if len(sub) != 89 or sub["Patient"].nunique() != 89:
        raise ValueError(f"Expected 89 distinct Medial patients; got {len(sub)} rows and {sub['Patient'].nunique()} patients")
    cols = TRAD_COLS + BRIGHT_COLS + ["PLM_SZ", "PLM_Total"]
    missing = sorted(set(cols).difference(sub.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")
    if not sub[cols].notna().all().all():
        raise ValueError("Missing predictor or endpoint values")
    return sub

def calculate_oof_auc(X: np.ndarray, y: np.ndarray):
    """Original Stage4 five-fold recipe and 500-bootstrap OOF confidence interval."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof = np.zeros(len(y), dtype=float)
    for tr, te in cv.split(X, y):
        scaler = StandardScaler()
        train = scaler.fit_transform(X[tr]); test = scaler.transform(X[te])
        lr = LogisticRegression(C=0.2, random_state=42)
        rf = RandomForestClassifier(n_estimators=100, max_depth=4, min_samples_leaf=5, random_state=42)
        lr.fit(train, y[tr]); rf.fit(train, y[tr])
        oof[te] = 0.5 * lr.predict_proba(test)[:, 1] + 0.5 * rf.predict_proba(test)[:, 1]
    auc = roc_auc_score(y, oof)
    rng = np.random.RandomState(42)
    boot = []
    for _ in range(500):
        idx = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[idx])) == 2:
            boot.append(roc_auc_score(y[idx], oof[idx]))
    if not boot:
        raise RuntimeError("Bootstrap resamples lack both endpoint classes")
    low, high = np.percentile(boot, [2.5, 97.5])
    return {"auc": auc, "low": low, "high": high, "oof": oof}

def export_figure(fig, stem: Path):
    for ext in ("pdf", "svg", "tiff"):
        kw = {"bbox_inches": "tight"}
        if ext == "tiff":
            kw.update(dpi=350, pil_kwargs={"compression": "tiff_lzw"})
        fig.savefig(stem.with_suffix("." + ext), **kw)
    plt.close(fig)

def run_three_model_oof(df: pd.DataFrame, y: np.ndarray, models: dict):
    """Return three-model OOF metrics and patient-level probabilities for one endpoint."""
    results = {name: calculate_oof_auc(df[cols].to_numpy(float), y)
               for name, cols in models.items()}
    metrics = pd.DataFrame([
        {"Model": name, "N_patients": len(y), "N_positive": int(y.sum()),
         "Positive_prevalence": float(y.mean()),
         "AUC": d["auc"], "CI_low": d["low"], "CI_high": d["high"]}
        for name, d in results.items()
    ])
    probabilities = df[["Patient", "Region", "SpecimenID"]].copy()
    for name, d in results.items():
        probabilities["OOF_" + name.replace(" ", "_")] = d["oof"]
    return results, metrics, probabilities

def main(args):
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    df = load_cohort(Path(args.scores), Path(args.features))
    models = {"Traditional grading": TRAD_COLS,
              "Brightfield zonal pathomics": BRIGHT_COLS,
              "Combined": TRAD_COLS + BRIGHT_COLS}

    # Figure 4c: superficial PLM endpoint. Keep original figure and filenames.
    y = (df["PLM_SZ"] >= 3).astype(int).to_numpy()
    results, metrics, out_probs = run_three_model_oof(df, y, models)
    metrics[["Model", "N_patients", "AUC", "CI_low", "CI_high"]].to_csv(
        out / "Figure4c_OOF_AUC_Summary.csv", index=False)
    out_probs.insert(3, "Observed_PLMSZ_severe", y)
    out_probs.to_csv(out / "Figure4c_PatientLevel_OOF_Predictions.csv", index=False)
    colors = {"Traditional grading": "#D9A441", "Brightfield zonal pathomics": "#20A39E", "Combined": "#D65A8A"}
    styles = {"Traditional grading": "--", "Brightfield zonal pathomics": "-.", "Combined": "-"}
    fig, ax = plt.subplots(figsize=(5.35, 4.75))
    for name, d in results.items():
        fpr, tpr, _ = roc_curve(y, d["oof"])
        ax.plot(fpr, tpr, color=colors[name], ls=styles[name], lw=2.7,
                label=f"{name}\nAUC {d['auc']:.3f} ({d['low']:.2f}–{d['high']:.2f})")
    ax.plot([0, 1], [0, 1], ls=":", color="#777777", lw=1.3)
    ax.set(xlim=(0, 1), ylim=(0, 1.02), xlabel="False-positive rate (1 − specificity)",
           ylabel="True-positive rate (sensitivity)")
    ax.set_title("Patient-level 5-fold OOF ROC", loc="left")
    ax.legend(loc="lower right", fontsize=8.2); ax.grid(alpha=0.18); fig.tight_layout()
    export_figure(fig, out / "Figure4c_ROC")

    # Supplementary numerical analysis: full-thickness PLM endpoint.
    # The same three model definitions and OOF procedure are used; this does
    # not replace the superficial-endpoint Figure 4c or create a new panel.
    imp_y = (df["PLM_Total"] >= 5).astype(int).to_numpy()
    full_results, full_metrics, full_probs = run_three_model_oof(df, imp_y, models)
    full_metrics.insert(0, "Endpoint", "PLM_Total >= 5")
    full_metrics.to_csv(out / "Figure4_FullThickness_PLMTotal5_OOF_AUC_Summary.csv", index=False)
    full_probs.insert(3, "Observed_PLMTotal_severe", imp_y)
    full_probs.to_csv(out / "Figure4_FullThickness_PLMTotal5_PatientLevel_OOF_Predictions.csv", index=False)

    # Figure 4d: exploratory full-cohort impurity importance for PLM_Total >= 5.
    X = StandardScaler().fit_transform(df[BRIGHT_COLS].to_numpy(float))
    rf = RandomForestClassifier(n_estimators=200, max_depth=5, random_state=42)
    rf.fit(X, imp_y)
    importance = pd.DataFrame({"Feature": BRIGHT_COLS, "Importance": rf.feature_importances_}).sort_values(
        "Importance", ascending=False).reset_index(drop=True)
    importance.insert(0, "Rank", np.arange(1, len(importance) + 1))
    importance.to_csv(out / "Figure4d_All27_FeatureImportance.csv", index=False)
    top10 = importance.head(10)
    cols = ["#6A3D9A" if ("SurfaceAngle" in x or "Grad_Angle" in x) else
            "#D65A8A" if "Basophilia" in x else "#20A39E" for x in top10.Feature]
    fig, ax = plt.subplots(figsize=(6.25, 4.25))
    ypos = np.arange(len(top10)); ax.barh(ypos, top10.Importance, color=cols, edgecolor="white")
    ax.set_yticks(ypos, top10.Feature.str.replace("_", " ", regex=False), fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("Mean decrease in impurity (RF feature importance)")
    ax.set_title("Random-forest feature importance", loc="left")
    ax.grid(axis="x", alpha=0.18)
    vmax = top10.Importance.max(); ax.set_xlim(0, vmax * 1.24)
    for i, v in enumerate(top10.Importance):
        ax.text(v + vmax * .018, i, f"{v:.3f}", va="center", fontsize=9)
    fig.tight_layout(); export_figure(fig, out / "Figure4d_RF_Importance")
    print("Figure 4c — PLM_SZ >= 3")
    print(metrics.to_string(index=False))
    print("\nSupplementary full-thickness analysis — PLM_Total >= 5")
    print(full_metrics.to_string(index=False))
    return metrics

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, help="Corrected OA_Specimen_Scores_177.csv")
    parser.add_argument("--features", required=True, help="Validated final OA_Zonal_Pathomics_Master_177_20261101.csv")
    parser.add_argument("--output-dir", default="outputs/simtk")
    main(parser.parse_args())
    
