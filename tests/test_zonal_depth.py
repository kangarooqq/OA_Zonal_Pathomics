"""Synthetic functional tests for simtk/zonal_depth.py.

These tests verify numerical and structural behaviour on an idealised rectangular
cartilage/bone geometry. They do NOT validate segmentation on real histology.
"""
from pathlib import Path
import importlib.util

import numpy as np
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "simtk" / "zonal_depth.py"


def load_depth_function():
    assert MODULE_PATH.is_file(), f"Missing source file: {MODULE_PATH}"
    spec = importlib.util.spec_from_file_location("zonal_depth_under_test", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_zonal_depth_field_v2


@pytest.fixture
def synthetic_masks():
    """Horizontal cartilage band immediately above a substantial bone band."""
    cartilage = np.zeros((240, 300), dtype=np.uint8)
    bone = np.zeros_like(cartilage)
    cartilage[45:165, 45:255] = 1
    bone[165:220, 35:265] = 1
    return cartilage, bone


def test_output_structure_and_background(synthetic_masks):
    cartilage, bone = synthetic_masks
    surface, tide, depth, zones, thickness = load_depth_function()(cartilage, bone)

    assert surface.shape == tide.shape == depth.shape == zones.shape == cartilage.shape
    assert surface.dtype == tide.dtype == zones.dtype == np.uint8
    assert depth.dtype == np.float32
    assert np.any(surface == 1), "Surface boundary was not detected"
    assert np.any(tide == 1), "Tidemark boundary was not detected"
    assert set(np.unique(surface)).issubset({0, 1})
    assert set(np.unique(tide)).issubset({0, 1})
    assert np.all(zones[cartilage == 0] == 0)
    assert np.all(np.isnan(depth[cartilage == 0]))
    assert np.all(np.isfinite(depth[cartilage == 1]))
    assert np.all((depth[cartilage == 1] >= 0) & (depth[cartilage == 1] <= 1))
    assert np.isfinite(thickness) and thickness > 0


def test_three_zones_and_thresholds(synthetic_masks):
    cartilage, bone = synthetic_masks
    _, _, depth, zones, _ = load_depth_function()(cartilage, bone)
    valid = cartilage == 1

    assert {1, 2, 3}.issubset(set(np.unique(zones[valid])))
    expected = np.zeros_like(zones)
    expected[valid & (depth <= 0.15)] = 1
    expected[valid & (depth > 0.15) & (depth <= 0.50)] = 2
    expected[valid & (depth > 0.50)] = 3
    np.testing.assert_array_equal(zones, expected)

    # In the central column, surface -> middle -> deep is top -> bottom.
    assert zones[45, 150] == 1
    assert zones[100, 150] == 2
    assert zones[164, 150] == 3
    assert depth[45, 150] < depth[100, 150] < depth[164, 150]


def test_mpp_scales_thickness_not_normalized_depth(synthetic_masks):
    cartilage, bone = synthetic_masks
    fn = load_depth_function()
    result_a = fn(cartilage, bone, mpp=6.88)
    result_b = fn(cartilage, bone, mpp=3.44)

    np.testing.assert_array_equal(result_a[0], result_b[0])
    np.testing.assert_array_equal(result_a[1], result_b[1])
    np.testing.assert_allclose(result_a[2], result_b[2], atol=1e-6, equal_nan=True)
    np.testing.assert_array_equal(result_a[3], result_b[3])
    assert result_b[4] == pytest.approx(result_a[4] / 2)
