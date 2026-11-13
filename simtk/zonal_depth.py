"""Stage2 B064: surface/tidemark boundary and three-zone depth map.

The original algorithm is preserved; a caller must supply validated binary cartilage
and bone masks at the source image resolution. mpp defaults to historical 6.88.
This module has not been validated against supplied microscopy images in this package.
"""
import numpy as np
import cv2
from scipy.ndimage import distance_transform_edt

def build_zonal_depth_field_v2(cart_mask, subbone_raw, mpp=6.88):
    """Separate cartilage surface/tidemark after morphological removal of bone-shell artefacts."""
    # Historical annotation removed during English-language code curation.
    sub_eroded = cv2.morphologyEx(
        subbone_raw,
        cv2.MORPH_OPEN,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (19, 19)),
    )
    num_b, lbl_b, stats_b, _ = cv2.connectedComponentsWithStats(sub_eroded, 8)
    true_bone = (
        lbl_b == (1 + int(np.argmax(stats_b[1:, cv2.CC_STAT_AREA])))
    ).astype(np.uint8)

    # Historical annotation removed during English-language code curation.
    cnts, _ = cv2.findContours(
        cart_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
    )
    main_cnt = max(cnts, key=cv2.contourArea)[:, 0, :]  # shape: (N, 2) -> (x, y)
    cx_arr, cy_arr = main_cnt[:, 0], main_cnt[:, 1]

    # Historical annotation removed during English-language code curation.
    dist_to_bone = distance_transform_edt(1 - true_bone)
    is_tide = dist_to_bone[cy_arr, cx_arr] <= 18

    tidemark_mask = np.zeros_like(cart_mask)
    tidemark_mask[cy_arr[is_tide], cx_arr[is_tide]] = 1

    # Historical annotation removed during English-language code curation.
    # Historical annotation removed during English-language code curation.
    tide_indices = np.where(is_tide)[0]
    # Historical annotation removed during English-language code curation.
    d_from_tide = distance_transform_edt(1 - tidemark_mask)[cy_arr, cx_arr]
    apex_idx = int(np.argmax(d_from_tide))
    n_pts = len(main_cnt)
    shift = (n_pts // 2) - apex_idx
    rolled_cnt = np.roll(main_cnt, shift, axis=0)
    rolled_tide = np.roll(is_tide, shift)

    # Historical annotation removed during English-language code curation.
    non_tide_idx = np.where(~rolled_tide)[0]
    # Historical annotation removed during English-language code curation.
    k_start = int(len(non_tide_idx) * 0.24)
    k_end = int(len(non_tide_idx) * 0.76)
    surf_indices = non_tide_idx[k_start:k_end]

    surface_mask = np.zeros_like(cart_mask)
    sx_pts = rolled_cnt[surf_indices, 0]
    sy_pts = rolled_cnt[surf_indices, 1]
    surface_mask[sy_pts, sx_pts] = 1

    # Historical annotation removed during English-language code curation.
    d_surf = distance_transform_edt(1 - surface_mask) * mpp
    d_tide = distance_transform_edt(1 - tidemark_mask) * mpp
    depth_norm = np.full(cart_mask.shape, np.nan, dtype=np.float32)
    valid = cart_mask == 1
    depth_norm[valid] = d_surf[valid] / (d_surf[valid] + d_tide[valid] + 1e-6)

    zone_map = np.zeros(cart_mask.shape, dtype=np.uint8)
    zone_map[valid & (depth_norm <= 0.15)] = 1
    zone_map[valid & (depth_norm > 0.15) & (depth_norm <= 0.50)] = 2
    zone_map[valid & (depth_norm > 0.50)] = 3

    med_thick_um = float(np.median(d_surf[tidemark_mask == 1]))
    return surface_mask, tidemark_mask, depth_norm, zone_map, med_thick_um
