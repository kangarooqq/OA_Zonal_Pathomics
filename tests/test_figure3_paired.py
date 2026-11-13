
"""Synthetic functional tests for the Figure 3 paired analysis.

These tests validate the code using artificial patient-level data.
They do not reproduce the manuscript's results from real patient data.
"""

from pathlib import Path
import importlib.util

import numpy as np
import pandas as pd
import pytest
from scipy.stats import wilcoxon


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "simtk"
    / "figure3_paired.py"
)


def load_figure3_module():
    assert MODULE_PATH.is_file(), f"Missing source file: {MODULE_PATH}"

    spec = importlib.util.spec_from_file_location(
        "figure3_paired_under_test", MODULE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def figure3():
    return load_figure3_module()


@pytest.fixture
def paired_data(figure3):
    """Construct 87 complete CL/CM pairs for all seven outcomes."""
    records = []

    for patient in range(1, 88):
        for region in ("Lateral", "Medial"):
            row = {
                "Patient": f"P{patient:03d}",
                "Region": region,
            }

            for index, feature in enumerate(
                figure3.FEATURES + figure3.SCORE_COLS
            ):
                baseline = 10.0 + patient * 0.1 + index
                difference = 0.2 + (patient % 7) * 0.03

                row[feature] = (
                    baseline
                    if region == "Lateral"
                    else baseline + difference
                )

            records.append(row)

    return pd.DataFrame(records)


def test_87_complete_pairs_and_output(figure3, paired_data):
    """All seven outcomes must contain exactly 87 complete pairs."""
    summary, source_data = figure3.paired_summary(paired_data)

    expected_features = figure3.FEATURES + figure3.SCORE_COLS

    assert len(summary) == 7
    assert set(summary["Feature"]) == set(expected_features)
    assert (summary["N_pairs"] == 87).all()

    assert len(source_data) == 7 * 87
    assert {"Feature", "Patient", "CL", "CM"}.issubset(
        source_data.columns
    )

    for feature in expected_features:
        subset = source_data[source_data["Feature"] == feature]
        assert len(subset) == 87
        assert subset["Patient"].nunique() == 87

    assert summary["Wilcoxon_P"].notna().all()


def test_86_pairs_rejected_for_every_outcome(
    figure3, paired_data
):
    """Missing one side of any outcome must trigger an error."""
    for feature in figure3.FEATURES + figure3.SCORE_COLS:
        modified = paired_data.copy()

        mask = (
            (modified["Patient"] == "P001")
            & (modified["Region"] == "Medial")
        )
        modified.loc[mask, feature] = np.nan

        with pytest.raises(
            ValueError,
            match="expected 87 complete CL/CM pairs",
        ):
            figure3.paired_summary(modified)


def test_duplicate_patient_region_rejected(
    figure3, paired_data
):
    """Duplicate patient-region records must not be accepted."""
    duplicate = pd.concat(
        [paired_data, paired_data.iloc[[0]]],
        ignore_index=True,
    )

    with pytest.raises(
        ValueError,
        match="Duplicate patient-region records",
    ):
        figure3.paired_summary(duplicate)


def test_wilcoxon_matches_independent_calculation(
    figure3, paired_data
):
    """Reported paired Wilcoxon statistics must match SciPy."""
    summary, source_data = figure3.paired_summary(paired_data)

    for feature in figure3.FEATURES + figure3.SCORE_COLS:
        subset = source_data[
            source_data["Feature"] == feature
        ]

        expected_stat, expected_p = wilcoxon(
            subset["CL"].to_numpy(dtype=float),
            subset["CM"].to_numpy(dtype=float),
            alternative="two-sided",
            zero_method="wilcox",
        )

        reported = summary.loc[
            summary["Feature"] == feature
        ].iloc[0]

        assert np.isclose(
            reported["Wilcoxon_stat"], expected_stat
        )
        assert np.isclose(
            reported["Wilcoxon_P"], expected_p
        )
        assert np.isclose(
            reported["Lateral_mean"], subset["CL"].mean()
        )
        assert np.isclose(
            reported["Medial_mean"], subset["CM"].mean()
        )
