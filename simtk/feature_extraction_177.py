"""SimTK Stage3 B025: complete original 177-specimen feature extraction candidate.

This is the latest *complete* feature extraction candidate found in the provided
archives, not a verified final publication run. Requires real TIFFs, classifier
training annotations, source CSVs and their correct paths; see docs/OPEN_ISSUES.md.
The original feature computations have not been tuned to target published numbers.
"""
# ==============================================================================
# Historical annotation removed during English-language code curation.
# Historical annotation removed during English-language code curation.
# Historical annotation removed during English-language code curation.
# ==============================================================================

import json
import os
from pathlib import Path
import time
import warnings
import cv2
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.ndimage import distance_transform_edt
from skimage import measure
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from sklearn.ensemble import RandomForestClassifier
import tifffile

warnings.filterwarnings("ignore")
t0 = time.time()

DATA_DIR = Path(os.environ["OA_SIMTK_DATA_DIR"])
df_img = pd.read_csv(DATA_DIR / "OA_Image_Properties_791.csv")
df_spec = pd.read_csv(DATA_DIR / "OA_Specimen_Scores_177.csv")

M_he = np.array(
    [[0.650, 0.704, 0.286], [0.072, 0.990, 0.105], [0.268, 0.570, 0.776]],
    dtype=np.float32,
)
M_he_inv = np.linalg.inv(M_he)


def white_balance_slide(img_rgb_in):
    img_f = img_rgb_in.astype(np.float32)
    gray = cv2.cvtColor(img_rgb_in, cv2.COLOR_RGB2GRAY)
    bg_mask = (gray >= np.percentile(gray, 80)) & (
        gray <= np.percentile(gray, 99)
    )
    bg_rgb = np.median(img_f[bg_mask], axis=0)
    return np.clip(img_f / (bg_rgb + 1e-5) * 252.0, 0, 255).astype(np.uint8)


def extract_multiscale_features(img_rgb_s):
    img_f = img_rgb_s.astype(np.float32) / 255.0
    od = -np.log10(np.clip(img_f, 1.0 / 255.0, 1.0))
    stains = np.dot(od, M_he_inv)
    h_ch = np.clip(stains[:, :, 0], 0, 3.0)
    e_ch = np.clip(stains[:, :, 1], 0, 3.0)
    base_channels = [
        img_f[:, :, 0],
        img_f[:, :, 1],
        img_f[:, :, 2],
        h_ch,
        e_ch,
        h_ch - e_ch,
    ]
    feat_list = []
    for sigma in [1.0, 3.0, 7.0, 15.0]:
        ksize = int(sigma * 4) | 1
        for ch in base_channels:
            feat_list.append(cv2.GaussianBlur(ch, (ksize, ksize), sigma))
        for ch in [img_f[:, :, 1], h_ch, e_ch]:
            mc = cv2.GaussianBlur(ch, (ksize, ksize), sigma)
            sq = cv2.GaussianBlur(ch * ch, (ksize, ksize), sigma)
            feat_list.append(np.sqrt(np.maximum(sq - mc * mc, 0)))
    return np.stack(feat_list, axis=-1)


# Historical annotation removed during English-language code curation.
if "rf_clf" not in globals():
    print("Kernel 36 rf_clf...")
    p60_raw_path = Path(
        df_img[
            (df_img["Patient"] == "P060")
            & (df_img["Region"] == "Medial")
            & (df_img["Stain"] == "HE")
            & (df_img["Section"] == "S1")
        ].iloc[0]["ImagePath"]
    )
    p60_img = tifffile.imread(p60_raw_path)
    with open(
        DATA_DIR / "P060_Medial_HE_S1_Cartilage_Manual_20261101.geojson",
        "r",
        encoding="utf-8",
    ) as f:
        geo_d = json.load(f)
    feats_g = (
        geo_d.get("features", [geo_d]) if isinstance(geo_d, dict) else geo_d
    )
    coords_g = feats_g[0]["geometry"]["coordinates"]
    polys_g = (
        coords_g
        if feats_g[0]["geometry"]["type"] == "MultiPolygon"
        else [coords_g]
    )
    gt_p60_full = np.zeros(p60_img.shape[:2], dtype=np.uint8)
    for poly in polys_g:
        cv2.fillPoly(
            gt_p60_full,
            [np.round(np.array(poly[0], dtype=np.float32)).astype(np.int32)],
            1,
        )
    gt_p60_full = ndimage.binary_fill_holes(gt_p60_full).astype(np.uint8)

    p60_1500 = cv2.resize(p60_img, (1500, 1500), interpolation=cv2.INTER_AREA)
    gt_p60_1500 = cv2.resize(
        gt_p60_full, (1500, 1500), interpolation=cv2.INTER_NEAREST
    )
    pred_v2_1500 = cv2.resize(
        tifffile.imread(DATA_DIR / "P060_HE_RTrees_4Class_v2_20261101.ome.tif"),
        (1500, 1500),
        interpolation=cv2.INTER_NEAREST,
    )
    y_tr_map = pred_v2_1500.copy()
    y_tr_map[gt_p60_1500 == 1] = 2
    y_tr_map[(gt_p60_1500 == 0) & (y_tr_map == 2)] = 1

    X_p60 = extract_multiscale_features(p60_1500)
    np.random.seed(42)
    X_s, y_s = [], []
    for cid in [0, 1, 2, 3]:
        yc, xc = np.where(y_tr_map == cid)
        idx_c = np.random.choice(
            len(yc), size=min(15000, len(yc)), replace=False
        )
        X_s.append(X_p60[yc[idx_c], xc[idx_c]])
        y_s.append(np.full(len(idx_c), cid, dtype=np.uint8))

    p50_p = Path(
        df_img[
            (df_img["Patient"] == "P050")
            & (df_img["Region"] == "Medial")
            & (df_img["Stain"] == "HE")
            & (df_img["Section"] == "S1")
        ].iloc[0]["ImagePath"]
    )
    X_p50 = extract_multiscale_features(
        cv2.resize(
            tifffile.imread(p50_p), (1500, 1500), interpolation=cv2.INTER_AREA
        )
    )
    yb, xb = np.where((X_p50[:, :, 5] > 0.04) & (X_p50[:, :, 1] < 0.85))
    idx_b = np.random.choice(len(yb), size=min(10000, len(yb)), replace=False)
    X_s.append(X_p50[yb[idx_b], xb[idx_b]])
    y_s.append(np.full(len(idx_b), 2, dtype=np.uint8))

    rf_clf = RandomForestClassifier(
        n_estimators=100,
        max_depth=16,
        min_samples_leaf=10,
        n_jobs=-1,
        random_state=42,
    )
    rf_clf.fit(np.vstack(X_s), np.concatenate(y_s))


def segment_cartilage_wb(img_full, clf, target_size=1500):
    im_raw = cv2.resize(
        img_full, (target_size, target_size), interpolation=cv2.INTER_AREA
    )
    im_s = white_balance_slide(im_raw)
    gray = cv2.cvtColor(im_s, cv2.COLOR_RGB2GRAY)
    tissue_raw = cv2.morphologyEx(
        (gray < 232).astype(np.uint8),
        cv2.MORPH_OPEN,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)),
    )
    tissue_soft = (gray < 243).astype(np.uint8)
    spec_solid = ndimage.binary_fill_holes(
        cv2.morphologyEx(
            tissue_raw,
            cv2.MORPH_CLOSE,
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (65, 65)),
        )
    ).astype(np.uint8)
    num_s, lbl_s, stats_s, _ = cv2.connectedComponentsWithStats(spec_solid, 8)
    spec_solid = (
        lbl_s == (1 + int(np.argmax(stats_s[1:, cv2.CC_STAT_AREA])))
    ).astype(np.uint8)

    inner_spec = cv2.erode(
        spec_solid, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (19, 19))
    )
    raw_holes = cv2.morphologyEx(
        ((gray >= 230) & (inner_spec == 1)).astype(np.uint8),
        cv2.MORPH_OPEN,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7)),
    )
    num_h, lbl_h, stats_h, _ = cv2.connectedComponentsWithStats(raw_holes, 8)
    true_marrow = np.zeros_like(raw_holes)
    for hid in range(1, num_h):
        if stats_h[hid, cv2.CC_STAT_AREA] >= 200:
            true_marrow[lbl_h == hid] = 1

    feats = extract_multiscale_features(im_s)
    probs = clf.predict_proba(feats.reshape(-1, feats.shape[-1])).reshape(
        target_size, target_size, 4
    )
    marrow_hull = ndimage.binary_fill_holes(
        cv2.morphologyEx(
            cv2.dilate(
                true_marrow,
                cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (45, 45)),
            ),
            cv2.MORPH_CLOSE,
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (75, 75)),
        )
    ).astype(np.uint8)
    near_marrow = cv2.dilate(
        marrow_hull, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35, 35))
    )
    subbone_plate = (
        (near_marrow == 1)
        & ((probs[:, :, 1] + probs[:, :, 3]) > probs[:, :, 2] * 1.1)
    ).astype(np.uint8)
    subbone_zone = ((marrow_hull == 1) | (subbone_plate == 1)) & (
        spec_solid == 1
    )

    cart_cand = cv2.morphologyEx(
        ((tissue_soft == 1) & (spec_solid == 1) & (subbone_zone == 0)).astype(
            np.uint8
        ),
        cv2.MORPH_OPEN,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11)),
    )
    num_c, lbl_c, stats_c, _ = cv2.connectedComponentsWithStats(cart_cand, 8)
    clean_c = ndimage.binary_fill_holes(
        cv2.morphologyEx(
            (lbl_c == (1 + int(np.argmax(stats_c[1:, cv2.CC_STAT_AREA])))).astype(
                np.uint8
            ),
            cv2.MORPH_CLOSE,
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (19, 19)),
        )
    ).astype(np.uint8)
    clean_b = ((spec_solid == 1) & (clean_c == 0)).astype(np.uint8)
    return im_s, clean_c, clean_b


def build_zonal_depth_field_v2(cart_mask, subbone_raw, mpp=6.88):
    sub_eroded = cv2.morphologyEx(
        subbone_raw,
        cv2.MORPH_OPEN,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (19, 19)),
    )
    num_b, lbl_b, stats_b, _ = cv2.connectedComponentsWithStats(sub_eroded, 8)
    true_bone = (
        lbl_b == (1 + int(np.argmax(stats_b[1:, cv2.CC_STAT_AREA])))
    ).astype(np.uint8)

    cnts, _ = cv2.findContours(
        cart_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
    )
    main_cnt = max(cnts, key=cv2.contourArea)[:, 0, :]
    cx_arr, cy_arr = main_cnt[:, 0], main_cnt[:, 1]
    is_tide = distance_transform_edt(1 - true_bone)[cy_arr, cx_arr] <= 18

    tidemark_mask = np.zeros_like(cart_mask)
    tidemark_mask[cy_arr[is_tide], cx_arr[is_tide]] = 1

    d_from_tide = distance_transform_edt(1 - tidemark_mask)[cy_arr, cx_arr]
    apex_idx = int(np.argmax(d_from_tide))
    shift = (len(main_cnt) // 2) - apex_idx
    rolled_cnt = np.roll(main_cnt, shift, axis=0)
    rolled_tide = np.roll(is_tide, shift)

    non_tide_idx = np.where(~rolled_tide)[0]
    surf_indices = non_tide_idx[
        int(len(non_tide_idx) * 0.24) : int(len(non_tide_idx) * 0.76)
    ]
    surface_mask = np.zeros_like(cart_mask)
    surface_mask[rolled_cnt[surf_indices, 1], rolled_cnt[surf_indices, 0]] = 1

    d_surf = distance_transform_edt(1 - surface_mask) * mpp
    d_tide = distance_transform_edt(1 - tidemark_mask) * mpp
    depth_norm = np.full(cart_mask.shape, np.nan, dtype=np.float32)
    valid = cart_mask == 1
    depth_norm[valid] = d_surf[valid] / (d_surf[valid] + d_tide[valid] + 1e-6)

    zone_map = np.zeros(cart_mask.shape, dtype=np.uint8)
    zone_map[valid & (depth_norm <= 0.15)] = 1
    zone_map[valid & (depth_norm > 0.15) & (depth_norm <= 0.50)] = 2
    zone_map[valid & (depth_norm > 0.50)] = 3
    return (
        surface_mask,
        tidemark_mask,
        depth_norm,
        zone_map,
        float(np.median(d_surf[tidemark_mask == 1])),
    )


def extract_specimen_master_features(pat, reg, sub_imgs, clf):
    he_rows = sub_imgs[sub_imgs["Stain"] == "HE"].sort_values("Section")
    sf_rows = sub_imgs[sub_imgs["Stain"] == "SafO"].sort_values("Section")
    pl_rows = sub_imgs[sub_imgs["Stain"] == "PLM"]
    if len(he_rows) == 0:
        return None

    p_he = Path(he_rows.iloc[0]["ImagePath"])
    im_he_full = tifffile.imread(p_he)
    mpp_1500 = im_he_full.shape[0] / 1500.0 * 1.72

    im_s, c_mask, b_mask = segment_cartilage_wb(
        im_he_full, clf, target_size=1500
    )
    s_mask, t_mask, d_field, z_map, thick_um = build_zonal_depth_field_v2(
        c_mask, b_mask, mpp=mpp_1500
    )

    valid_c = c_mask == 1
    area_mm2 = valid_c.sum() * (mpp_1500**2) / 1e6

    sy_pts, sx_pts = np.where(s_mask > 0)
    if len(sy_pts) > 10:
        chord = np.hypot(
            sy_pts.max() - sy_pts.min(), sx_pts.max() - sx_pts.min()
        )
        surf_roughness = float(len(sy_pts)) / (chord + 1e-5)
    else:
        surf_roughness = 1.0

    od_he = -np.log10(np.clip(im_s.astype(np.float32), 1.0, 255.0) / 255.0)
    stains_he = np.dot(od_he, M_he_inv)
    h_pure = np.clip(stains_he[:, :, 0], 0, 3.0)
    e_pure = np.clip(stains_he[:, :, 1], 0, 3.0)

    ys, xs = np.where(valid_c)
    y0, y1, x0, x1 = ys.min(), ys.max(), xs.min(), xs.max()
    scale_roi = 2.0
    mpp_roi = mpp_1500 / scale_roi
    h_r, w_r = int((y1 - y0) * scale_roi), int((x1 - x0) * scale_roi)

    fy0 = int(y0 * im_he_full.shape[0] / 1500.0)
    fy1 = int(y1 * im_he_full.shape[0] / 1500.0)
    fx0 = int(x0 * im_he_full.shape[1] / 1500.0)
    fx1 = int(x1 * im_he_full.shape[1] / 1500.0)

    roi_rgb = white_balance_slide(
        cv2.resize(
            im_he_full[fy0:fy1, fx0:fx1],
            (w_r, h_r),
            interpolation=cv2.INTER_AREA,
        )
    )
    roi_z = cv2.resize(
        z_map[y0:y1, x0:x1], (w_r, h_r), interpolation=cv2.INTER_NEAREST
    )
    roi_d = cv2.resize(
        d_field[y0:y1, x0:x1], (w_r, h_r), interpolation=cv2.INTER_LINEAR
    )
    roi_m = (roi_z > 0).astype(np.uint8)

    od_r = -np.log10(np.clip(roi_rgb.astype(np.float32), 1.0, 255.0) / 255.0)
    h_r_ch = np.clip(np.dot(od_r, M_he_inv)[:, :, 0], 0, 3.0)
    gray_r = cv2.cvtColor(roi_rgb, cv2.COLOR_RGB2GRAY).astype(np.float32)

    cell_dark = (
        np.clip(cv2.GaussianBlur(gray_r, (31, 31), 7.5) - gray_r, 0, None)
        / 255.0
    )
    pcm_halo = np.clip(
        cv2.GaussianBlur(h_r_ch, (5, 5), 1.2)
        - cv2.GaussianBlur(h_r_ch, (43, 43), 11.0),
        0,
        None,
    )
    comb_sig = (
        cv2.GaussianBlur(pcm_halo + cell_dark * 0.45, (5, 5), 1.2) * roi_m
    )

    in_m = cv2.erode(roi_m, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    ch_bw = np.zeros_like(roi_m, dtype=np.uint8)
    for zid, pct_th, mr in [(1, 88.0, 2), (2, 86.0, 3), (3, 85.0, 3)]:
        zp = (roi_z == zid) & (in_m == 1)
        if zp.sum() > 50:
            bw_z = ((comb_sig > np.percentile(comb_sig[zp], pct_th)) & zp).astype(
                np.uint8
            )
            bw_z = cv2.morphologyEx(
                bw_z,
                cv2.MORPH_CLOSE,
                cv2.getStructuringElement(
                    cv2.MORPH_ELLIPSE, (mr * 2 + 1, mr * 2 + 1)
                ),
            )
            bw_z = ndimage.binary_fill_holes(bw_z).astype(np.uint8)
            ch_bw[bw_z == 1] = 1

    dist_c = distance_transform_edt(ch_bw)
    pk_c = peak_local_max(
        cv2.GaussianBlur(dist_c, (5, 5), 1.2), min_distance=7, labels=ch_bw
    )
    mk_c = np.zeros(dist_c.shape, dtype=bool)
    if len(pk_c) > 0:
        mk_c[tuple(pk_c.T)] = True
    mk_c, _ = ndimage.label(mk_c)
    lbl_ch = watershed(-dist_c, mk_c, mask=ch_bw)

    nuc_dk = (
        cv2.GaussianBlur(gray_r, (15, 15), 4.0)
        - cv2.GaussianBlur(gray_r, (3, 3), 0.5)
    ) * ch_bw
    pk_n = peak_local_max(nuc_dk, min_distance=3, threshold_abs=9.5)
    nuc_map = np.zeros_like(roi_m, dtype=np.uint8)
    if len(pk_n) > 0:
        nuc_map[tuple(pk_n.T)] = 1

    gy_r = cv2.Sobel(np.nan_to_num(roi_d, nan=0.0), cv2.CV_32F, 0, 1, ksize=11)
    gx_r = cv2.Sobel(np.nan_to_num(roi_d, nan=0.0), cv2.CV_32F, 1, 0, ksize=11)

    ch_list = []
    for p in measure.regionprops(lbl_ch, intensity_image=nuc_map):
        a_um2 = p.area * (mpp_roi**2)
        if not (30.0 <= a_um2 <= 2400.0):
            continue
        cy, cx = int(round(p.centroid[0])), int(round(p.centroid[1]))
        zid = int(roi_z[cy, cx])
        if zid == 0:
            continue
        ar = max(p.axis_major_length, 1.0) / max(p.axis_minor_length, 1.0)
        cvy, cvx = np.cos(p.orientation), -np.sin(p.orientation)
        lny, lnx = gy_r[cy, cx], gx_r[cy, cx]
        nm = np.hypot(lny, lnx) + 1e-6
        lny, lnx = lny / nm, lnx / nm
        cos_norm = np.clip(abs(cvy * lny + cvx * lnx), 0.0, 1.0)
        surf_ang = 90.0 - np.degrees(np.arcsin(cos_norm))
        n_nuc = max(1, int(round(p.intensity_mean * p.area)))
        ch_list.append((zid, a_um2, ar, surf_ang, n_nuc, int(n_nuc >= 2)))

    df_c = pd.DataFrame(
        ch_list, columns=["ZoneID", "Area", "AR", "SurfAngle", "NucCount", "Clone"]
    )

    sf_pg = {1: np.nan, 2: np.nan, 3: np.nan}
    if len(sf_rows) > 0:
        im_sf = white_balance_slide(
            cv2.resize(
                tifffile.imread(Path(sf_rows.iloc[0]["ImagePath"])),
                (1500, 1500),
                interpolation=cv2.INTER_AREA,
            )
        )
        od_sf = -np.log10(np.clip(im_sf.astype(np.float32), 1.0, 255.0) / 255.0)
        pg_map = np.clip(od_sf[:, :, 1] - 0.5 * od_sf[:, :, 0], 0, 3.0)
        t_sf_m = (cv2.cvtColor(im_sf, cv2.COLOR_RGB2GRAY) < 232).astype(np.uint8)
        for zid in [1, 2, 3]:
            mz = (z_map == zid) & (t_sf_m == 1)
            if mz.sum() > 20:
                sf_pg[zid] = float(np.median(pg_map[mz]))

    plm_int = {1: np.nan, 2: np.nan, 3: np.nan}
    plm_coh = {1: np.nan, 2: np.nan, 3: np.nan}
    if len(pl_rows) > 0:
        im_pl = tifffile.imread(Path(pl_rows.iloc[0]["ImagePath"]))
        pl_g = cv2.cvtColor(im_pl, cv2.COLOR_RGB2GRAY).astype(np.float32)
        pl_top = np.clip(pl_g - cv2.GaussianBlur(pl_g, (201, 201), 0) - 2.5, 0, None)
        gxp = cv2.Sobel(pl_g, cv2.CV_32F, 1, 0, ksize=5)
        gyp = cv2.Sobel(pl_g, cv2.CV_32F, 0, 1, ksize=5)
        Jxx = cv2.GaussianBlur(gxp * gxp, (15, 15), 0)
        Jyy = cv2.GaussianBlur(gyp * gyp, (15, 15), 0)
        Jxy = cv2.GaussianBlur(gxp * gyp, (15, 15), 0)
        pl_c = np.sqrt((Jxx - Jyy) ** 2 + 4 * Jxy**2) / (Jxx + Jyy + 1e-5)

        z_pl = cv2.resize(z_map, (938, 938), interpolation=cv2.INTER_NEAREST)
        yz, xz = np.where(z_pl > 0)
        if len(yz) > 50:
            std_p = cv2.GaussianBlur(np.hypot(gxp, gyp), (25, 25), 0)
            yp, xp = np.where(std_p > np.percentile(std_p, 82))
            M_p = np.float32(
                [
                    [1, 0, int(np.percentile(xp, 25) - np.median(xz))],
                    [0, 1, int(np.percentile(yp, 25) - np.median(yz))],
                ]
            )
            z_on_p = cv2.warpAffine(
                z_pl, M_p, (pl_g.shape[1], pl_g.shape[0]), flags=cv2.INTER_NEAREST
            )
            for zid in [1, 2, 3]:
                mp = z_on_p == zid
                if mp.sum() > 20:
                    plm_int[zid] = float(np.mean(pl_top[mp]))
                    plm_coh[zid] = float(np.median(pl_c[mp]))

    rec = {
        "Patient": pat,
        "Region": reg,
        "Cartilage_Area_mm2": round(area_mm2, 3),
        "Median_Thickness_um": round(thick_um, 1),
        "Surface_Roughness": round(surf_roughness, 3),
        "Total_Chondrons": len(df_c),
        "Total_Nuclei": int(df_c["NucCount"].sum()) if len(df_c) > 0 else 0,
    }

    for zid, zname in [(1, "SZ"), (2, "MZ"), (3, "DZ")]:
        zm = z_map == zid
        z_area = max(zm.sum() * (mpp_1500**2) / 1e6, 1e-4)
        sc = df_c[df_c["ZoneID"] == zid]

        rec[f"{zname}_HE_Basophilia"] = round(
            float(np.median(h_pure[zm])) if zm.any() else np.nan, 4
        )
        rec[f"{zname}_HE_Eosinophilia"] = round(
            float(np.median(e_pure[zm])) if zm.any() else np.nan, 4
        )
        rec[f"{zname}_SafO_PG"] = (
            round(sf_pg[zid], 4) if pd.notna(sf_pg[zid]) else np.nan
        )
        rec[f"{zname}_ChondronDensity_mm2"] = round(len(sc) / z_area, 1)
        rec[f"{zname}_ChondronArea_um2"] = (
            round(float(sc["Area"].mean()), 1) if len(sc) > 0 else np.nan
        )
        rec[f"{zname}_AspectRatio"] = (
            round(float(sc["AR"].mean()), 3) if len(sc) > 0 else np.nan
        )
        rec[f"{zname}_SurfaceAngle_deg"] = (
            round(float(sc["SurfAngle"].mean()), 2) if len(sc) > 0 else np.nan
        )
        rec[f"{zname}_NucleiPerChondron"] = (
            round(float(sc["NucCount"].mean()), 3) if len(sc) > 0 else np.nan
        )
        rec[f"{zname}_CloningRate"] = (
            round(float(sc["Clone"].mean()), 3) if len(sc) > 0 else np.nan
        )
        rec[f"{zname}_PLM_Intensity"] = (
            round(plm_int[zid], 3) if pd.notna(plm_int[zid]) else np.nan
        )
        rec[f"{zname}_PLM_Coherence"] = (
            round(plm_coh[zid], 3) if pd.notna(plm_coh[zid]) else np.nan
        )

    rec["Grad_HE_Baso_DZ_SZ"] = round(
        rec["DZ_HE_Basophilia"] - rec["SZ_HE_Basophilia"], 4
    )
    rec["Grad_SafO_PG_DZ_SZ"] = (
        round(rec["DZ_SafO_PG"] - rec["SZ_SafO_PG"], 4)
        if pd.notna(rec["DZ_SafO_PG"]) and pd.notna(rec["SZ_SafO_PG"])
        else np.nan
    )
    rec["Grad_ChondronArea_DZ_SZ"] = (
        round(rec["DZ_ChondronArea_um2"] - rec["SZ_ChondronArea_um2"], 1)
        if pd.notna(rec["DZ_ChondronArea_um2"])
        and pd.notna(rec["SZ_ChondronArea_um2"])
        else np.nan
    )
    rec["Grad_Angle_DZ_SZ"] = (
        round(rec["DZ_SurfaceAngle_deg"] - rec["SZ_SurfaceAngle_deg"], 2)
        if pd.notna(rec["DZ_SurfaceAngle_deg"])
        and pd.notna(rec["SZ_SurfaceAngle_deg"])
        else np.nan
    )
    return rec


print(" 177 ...")
all_feats = []
total_n = len(df_spec)
for idx, s_row in df_spec.iterrows():
    pat, reg = s_row["Patient"], s_row["Region"]
    sub_imgs = df_img[
        (df_img["Patient"] == pat) & (df_img["Region"] == reg)
    ]
    try:
        f_rec = extract_specimen_master_features(pat, reg, sub_imgs, rf_clf)
        if f_rec is not None:
            all_feats.append(f_rec)
    except Exception as e:
        print(f"   [Warning] {pat} {reg} : {e}")

    if (idx + 1) % 25 == 0 or (idx + 1) == total_n:
        print(
            f"   Processed: {idx + 1} / {total_n} specimens (elapsed: {time.time() - t0:.1f}s)"
        )

df_zonal = pd.DataFrame(all_feats)
df_master = pd.merge(df_spec, df_zonal, on=["Patient", "Region"], how="inner")
out_master = Path(os.environ.get("OA_OUTPUT_DIR", "outputs/simtk")) / "OA_Zonal_Pathomics_Master_177.csv"
out_master.parent.mkdir(parents=True, exist_ok=True)
df_master.to_csv(out_master, index=False, encoding="utf-8-sig")

print(
    f"\n1. Cohort extraction completed; processed specimen count: {len(df_master)} / {total_n}, total elapsed: {time.time() - t0:.1f} seconds"
)
print(f"2. : {out_master.name} (: {df_master.shape[1]})")
print("\n3.  (n=177)  (SZ vs MZ vs DZ) ±")
preview_cols = [
    "SZ_HE_Basophilia",
    "MZ_HE_Basophilia",
    "DZ_HE_Basophilia",
    "SZ_SafO_PG",
    "MZ_SafO_PG",
    "DZ_SafO_PG",
    "SZ_ChondronArea_um2",
    "MZ_ChondronArea_um2",
    "DZ_ChondronArea_um2",
    "SZ_SurfaceAngle_deg",
    "MZ_SurfaceAngle_deg",
    "DZ_SurfaceAngle_deg",
    "SZ_NucleiPerChondron",
    "MZ_NucleiPerChondron",
    "DZ_NucleiPerChondron",
]
print(
    df_master[preview_cols]
    .describe()
    .loc[["count", "mean", "std", "50%"]]
    .round(3)
    .to_string()
)
