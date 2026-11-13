
"""Synthetic functional tests for Figure 4 predictive models.

These tests check computational behaviour only. They do not reproduce
the manuscript AUC values or validate the real patient cohort.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "simtk"
    / "figure4_models.py"
)


@pytest.fixture(scope="module")
def figure4():
    """Load the Figure 4 analysis module directly from its file."""
    assert MODULE_PATH.is_file(), f"Missing file: {MODULE_PATH}"

    spec = importlib.util.spec_from_file_location(
        "figure4_models_under_test", MODULE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def synthetic_cohort(figure4):
    """Create 89 independent patients with two distinct PLM endpoints."""
    rng = np.random.RandomState(42)
    n = 89

    data = {
        "Patient": [f"P{i:03d}" for i in range(n)],
        "Region": ["Medial"] * n,
        "SpecimenID": [f"S{i:03d}" for i in range(n)],
    }

    for feature in dict.fromkeys(
        figure4.TRAD_COLS + figure4.BRIGHT_COLS
    ):
        data[feature] = rng.normal(size=n)

    # Different labels ensure the two endpoints are not conflated.
    superficial = np.array([1] * 48 + [0] * 41)
    full_thickness = np.array([1] * 47 + [0] * 42)

    rng.shuffle(superficial)
    rng.shuffle(full_thickness)

    data["PLM_SZ"] = np.where(superficial == 1, 3, 2)
    data["PLM_Total"] = np.where(full_thickness == 1, 5, 4)

    return pd.DataFrame(data)


def test_three_model_oof_predictions(figure4, synthetic_cohort):
    """Each model must provide one valid OOF probability per patient."""
    models = {
        "Traditional grading": figure4.TRAD_COLS,
        "Brightfield zonal pathomics": figure4.BRIGHT_COLS,
        "Combined": figure4.TRAD_COLS + figure4.BRIGHT_COLS,
    }

    y = (synthetic_cohort["PLM_SZ"] >= 3).astype(int).to_numpy()

    results, metrics, probabilities = figure4.run_three_model_oof(
        synthetic_cohort, y, models
    )

    assert len(results) == 3
    assert len(metrics) == 3
    assert len(probabilities) == 89
    assert probabilities["Patient"].is_unique

    for name in models:
        oof = results[name]["oof"]

        assert oof.shape == (89,)
        assert np.isfinite(oof).all()
        assert ((oof >= 0) & (oof <= 1)).all()

        assert np.isclose(
            results[name]["auc"],
            roc_auc_score(y, oof),
        )

        assert (
            0 <= results[name]["low"]
            <= results[name]["high"] <= 1
        )

        column = "OOF_" + name.replace(" ", "_")
        assert column in probabilities.columns
        np.testing.assert_allclose(
            probabilities[column].to_numpy(), oof
        )


def test_full_thickness_endpoint_is_distinct(
    figure4, synthetic_cohort
):
    """PLM_Total >= 5 must not reuse superficial endpoint labels."""
    superficial = (
        synthetic_cohort["PLM_SZ"] >= 3
    ).astype(int).to_numpy()

    full_thickness = (
        synthetic_cohort["PLM_Total"] >= 5
    ).astype(int).to_numpy()

    assert len(superficial) == len(full_thickness) == 89
    assert superficial.sum() == 48
    assert full_thickness.sum() == 47
    assert not np.array_equal(superficial, full_thickness)

    models = {
        "Traditional grading": figure4.TRAD_COLS,
        "Brightfield zonal pathomics": figure4.BRIGHT_COLS,
        "Combined": figure4.TRAD_COLS + figure4.BRIGHT_COLS,
    }

    _, metrics, probabilities = figure4.run_three_model_oof(
        synthetic_cohort, full_thickness, models
    )

    assert len(metrics) == 3
    assert (metrics["N_patients"] == 89).all()
    assert (metrics["N_positive"] == 47).all()
    assert np.allclose(
        metrics["Positive_prevalence"],
        47 / 89,
    )
    assert len(probabilities) == 89


def test_oof_reproducibility(figure4):
    """Fixed random seeds should reproduce OOF predictions."""
    rng = np.random.RandomState(123)
    X = rng.normal(size=(30, 5))
    y = np.array([0, 1] * 15)

    first = figure4.calculate_oof_auc(X, y)
    second = figure4.calculate_oof_auc(X, y)

    np.testing.assert_allclose(first["oof"], second["oof"])
    assert first["auc"] == second["auc"]
    assert first["low"] == second["low"]
    assert first["high"] == second["high"]


def test_model_feature_definitions(figure4):
    """Confirm the intended traditional and brightfield predictors."""
    assert len(figure4.TRAD_COLS) == 5
    assert len(figure4.BRIGHT_COLS) == 27

    assert len(set(figure4.TRAD_COLS)) == 5
    assert len(set(figure4.BRIGHT_COLS)) == 27

    assert not (
        set(figure4.TRAD_COLS)
        & set(figure4.BRIGHT_COLS)
    )

    assert "PLM_SZ" not in figure4.TRAD_COLS
    assert "PLM_SZ" not in figure4.BRIGHT_COLS
    assert "PLM_Total" not in figure4.TRAD_COLS
    assert "PLM_Total" not in figure4.BRIGHT_COLS
