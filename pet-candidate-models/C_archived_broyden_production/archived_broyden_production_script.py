"""
MEMI/PET production pipeline for CR1000 multi-station study.

What this script does
- Loads station datasets (Stations 1,3,4,5,6,7,8).
- Applies Tg QC and (for Station 4) Tg recon handling using precomputed recon columns.
- Loads SVF metrics (hemispheric + directional 12-sector SVF).
- Computes sun azimuth for each timestamp (default Melbourne CBD lat/lon; editable).
- Applies directional SVF to shortwave forcing in MEMI: Kg_eff = Kg * SVF_dir(sun sector).
- Solves full 2-node MEMI using inverse-Jacobian Broyden (Option A / full solution).
- Computes PET as Ta_ref in a standard indoor reference env giving identical Ts (PET definition).
- Exports: Excel workbook + per-station CSV extracts.

Key references (Harvard)
- Broyden, C.G. (1965) A class of methods for solving nonlinear simultaneous equations.
- Höppe, P. (1984; 1999) MEMI and thermophysiological modeling foundations used in PET.
- Matzarakis, A., Mayer, H. and Iziomon, M.G. (1999) PET methodology and outdoor human-biometeorology.
- Bröde, P. et al. (2012) PET/thermophysiological modeling and radiation terms.
- ISO 7726 (1998) Globe thermometer and mean radiant temperature measurement principles.
- Fanger, P.O. (1970) Thermal comfort fundamentals (respiratory loss parameterisations commonly adopted).

Notes on “SVF wiring”
- Directional SVF is used to attenuate local shortwave load (via sun-azimuth sector).
- Hemispheric SVF is retained as a station descriptor and can be used for longwave partitioning
  if you later extend the model (sky vs surrounding surfaces); here the MEMI longwave uses MRT from globe.
"""

from __future__ import annotations

import math
import os
import re
import numpy as np
import pandas as pd
from dataclasses import dataclass

# ----------------------------
# USER SETTINGS
# ----------------------------

BASE_DIR = "."  # folder containing your station files + SVF excel files

# Solar position anchor (editable).
# Default: Melbourne CBD approximate. Update to your exact site if desired.
SITE_LAT = -37.8136
SITE_LON = 144.9631
SITE_TZ  = "Australia/Melbourne"  # used to interpret timestamps if you localize

HOT_DAYS = ["2018-01-06", "2018-01-18", "2018-01-19"]

# Input files (update only if your filenames differ)
STATION_FILES = {
    1: "CR1000_1_Table2 copy.dat.csv",
    3: "CR1000_3_Table2 copy.dat.csv",
    4: "CR1000_4_Table2_QC_TgReconstructed.csv",
    5: "CR1000_5_Table2 copy.dat.csv",
    6: "CR1000_6_Table2 copy.dat.csv",
    7: "CR1000_7_Table2_reconstructed_full_fixed.xlsx",
    8: "CR1000_8_Table2_reconstructed.xlsx",
}

SVF_DIR_FILE    = "SVF_directional_Stations1-8.xlsx"
SVF_TOTAL_FILE  = "SVF_results_Stations1-8.xlsx"  # optional; directional file often includes totals too

# Image-to-station mapping (two images per station)
IMG_TO_STATION = {
    "IMG_2258.JPG": 1, "IMG_6274.JPG": 1,
    "IMG_2295.JPG": 3, "IMG_6258.JPG": 3,
    "IMG_2280.JPG": 4, "IMG_6250.JPG": 4,
    "IMG_2270.JPG": 5, "IMG_6243.JPG": 5,
    "IMG_2305.JPG": 6, "IMG_6260.JPG": 6,
    "IMG_2264.JPG": 7, "IMG_6270.JPG": 7,
    "IMG_2320.JPG": 8, "IMG_6266.JPG": 8,
}

# ----------------------------
# CONSTANTS AND MODEL SETTINGS
# ----------------------------

SIGMA = 5.670374419e-8
P_ATM = 101325.0
L_V   = 2.42e6

# Tuned parameter set (as per your workflow)
# (Höppe, 1999; Matzarakis et al., 1999; Bröde et al., 2012; Fanger, 1970)
@dataclass
class Person2Node:
    M: float = 80.0
    emiss: float = 0.97
    alpha_sw: float = 0.60
    f_p: float = 0.70
    k_core_skin: float = 8.0

PERSON = Person2Node()

# ----------------------------
# IO HELPERS
# ----------------------------

TOA5_COLS = [
    "TIMESTAMP","RECORD","BattV_Avg","SlrW_Avg","WS_ms_Avg","AirTC_Avg","RH_Max","BB_mV_Avg","BB_temp_Avg"
]

def read_station_file(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()

    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path)
        return df

    # For your “copy.dat.csv” files: rows are wrapped in quotes as one field.
    # We skip the 4 metadata/header rows and split the quoted line ourselves.
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()

    # Find first data row after header block: assume first 4 lines are metadata/headers.
    data_lines = lines[4:]
    cleaned = []
    for ln in data_lines:
        ln = ln.strip().strip('"')
        if not ln:
            continue
        cleaned.append(ln.split(","))

    df = pd.DataFrame(cleaned, columns=TOA5_COLS[:len(cleaned[0])])

    return df

def coerce_station_types(df: pd.DataFrame) -> pd.DataFrame:
    # Standardize timestamp and numeric columns
    df = df.copy()
    df["TIMESTAMP"] = pd.to_datetime(df["TIMESTAMP"], errors="coerce")

    for c in ["RECORD","BattV_Avg","SlrW_Avg","WS_ms_Avg","AirTC_Avg","RH_Max","BB_mV_Avg","BB_temp_Avg"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # keep RECORD as metadata; primary index is TIMESTAMP
    df = df.dropna(subset=["TIMESTAMP"]).sort_values("TIMESTAMP")
    df = df.set_index("TIMESTAMP")
    return df

# ----------------------------
# SVF PROCESSING
# ----------------------------

def load_station_svf() -> pd.DataFrame:
    """
    Returns station-level SVF table with:
    - SVF_total (if present)
    - directional SVF sectors: SVF_000_030 ... SVF_330_360
    Aggregated across the two images per station (mean).
    """
    svf = pd.read_excel(os.path.join(BASE_DIR, SVF_DIR_FILE))
    # Expect a column containing image filename; try common names
    img_col = None
    for cand in ["image","Image","filename","Filename","file","File","img","IMG"]:
        if cand in svf.columns:
            img_col = cand
            break
    if img_col is None:
        # fallback: first column often image id
        img_col = svf.columns[0]

    svf["station"] = svf[img_col].map(IMG_TO_STATION)
    if svf["station"].isna().any():
        missing = svf.loc[svf["station"].isna(), img_col].unique()
        raise ValueError(f"SVF file contains images not mapped to stations: {missing}")

    # Identify directional columns
    dir_cols = [c for c in svf.columns if re.match(r"SVF_\d{3}_\d{3}", str(c))]
    keep_cols = ["station"] + dir_cols
    if "SVF_total" in svf.columns:
        keep_cols.append("SVF_total")

    out = svf[keep_cols].groupby("station").mean(numeric_only=True)
    return out

def sector_from_azimuth_deg(az_deg: float) -> str:
    """
    Returns directional SVF column name for 12-sector binning (30° sectors).
    az_deg: 0..360 from North clockwise.
    """
    az = az_deg % 360.0
    start = int(math.floor(az / 30.0) * 30)
    end = start + 30
    if end == 360:
        return f"SVF_{start:03d}_360"
    return f"SVF_{start:03d}_{end:03d}"

# ----------------------------
# SOLAR POSITION (azimuth)
# ----------------------------

def solar_azimuth_noaa(dt_utc: pd.Timestamp, lat: float, lon: float) -> float:
    """
    Compact NOAA-style solar azimuth approximation.

    Inputs:
    - dt_utc must be UTC.
    Output:
    - azimuth degrees from North clockwise.

    This is sufficient for sector-binning; for publication-grade solar geometry,
    replace with pysolar/astral and document the implementation.
    """
    # Julian day
    t = dt_utc.to_pydatetime()
    # Convert to day of year and fractional hour
    doy = t.timetuple().tm_yday
    frac_hour = t.hour + t.minute/60 + t.second/3600

    # fractional year gamma
    gamma = 2*math.pi/365 * (doy - 1 + (frac_hour - 12)/24)

    # equation of time (min)
    eot = 229.18 * (
        0.000075
        + 0.001868*math.cos(gamma)
        - 0.032077*math.sin(gamma)
        - 0.014615*math.cos(2*gamma)
        - 0.040849*math.sin(2*gamma)
    )

    # declination (rad)
    decl = (
        0.006918
        - 0.399912*math.cos(gamma)
        + 0.070257*math.sin(gamma)
        - 0.006758*math.cos(2*gamma)
        + 0.000907*math.sin(2*gamma)
        - 0.002697*math.cos(3*gamma)
        + 0.00148*math.sin(3*gamma)
    )

    # time offset (min)
    time_offset = eot + 4*lon  # UTC => no timezone correction

    # true solar time (min)
    tst = (t.hour*60 + t.minute + t.second/60 + time_offset) % 1440

    # hour angle (deg)
    ha = (tst / 4) - 180
    ha_rad = math.radians(ha)

    lat_rad = math.radians(lat)

    # solar zenith
    cos_zen = math.sin(lat_rad)*math.sin(decl) + math.cos(lat_rad)*math.cos(decl)*math.cos(ha_rad)
    cos_zen = min(1.0, max(-1.0, cos_zen))
    zen = math.acos(cos_zen)

    # azimuth (NOAA convention)
    sin_az = -math.sin(ha_rad)*math.cos(decl) / math.sin(zen) if math.sin(zen) != 0 else 0.0
    cos_az = (math.sin(decl) - math.sin(lat_rad)*math.cos(zen)) / (math.cos(lat_rad)*math.sin(zen)) if math.sin(zen) != 0 else 0.0
    az = math.degrees(math.atan2(sin_az, cos_az))
    az = (az + 360) % 360
    return az

# ----------------------------
# THERMOPHYSIOLOGY (MEMI + PET)
# ----------------------------

def sat_vp(Tc: float) -> float:
    # Saturation vapour pressure [Pa] (common Tetens form; document if you swap)
    return 610.78 * math.exp(17.27 * Tc / (Tc + 237.3))

def mrt_from_globe(Tg: float, Ta: float, v: float,
                   D_globe: float = 0.15, emiss_globe: float = 0.95) -> float:
    # Globe-based MRT (ISO 7726-style approach as used in many PET workflows)
    TgK = Tg + 273.15
    TaK = Ta + 273.15
    v_eff = max(v, 0.1)
    factor = 1.1e8 * (v_eff**0.6) / (emiss_globe * (D_globe**0.4))
    TmrtK4 = TgK**4 + factor*(TgK - TaK)
    return (TmrtK4**0.25) - 273.15

def memi2_residual(x: np.ndarray, Ta: float, RH: float, v: float, Kg: float, Tg: float, p: Person2Node) -> np.ndarray:
    """
    2-node MEMI residual vector F=[F_core, F_skin]=0.
    (Höppe, 1984; Höppe, 1999; Matzarakis et al., 1999; Bröde et al., 2012)
    """
    Tc, Ts = float(x[0]), float(x[1])
    Tmrt = mrt_from_globe(Tg, Ta, v)

    TsK   = Ts + 273.15
    TmrtK = Tmrt + 273.15

    # core-skin conduction
    Qc = p.k_core_skin * (Tc - Ts)

    # convection (typical PET-style power law; document if you swap)
    v_eff = max(v, 0.1)
    h_c = 8.3 * (v_eff**0.6)
    C_env = h_c * (Ts - Ta)

    # radiation: shortwave absorbed + longwave net
    K_abs  = p.f_p * p.alpha_sw * Kg
    LW_net = p.emiss * SIGMA * (TsK**4 - TmrtK**4)
    R_env  = LW_net - K_abs

    # evaporation (simple mass transfer coupling to h_c; keep consistent with your tuned workflow)
    p_air  = (RH/100.0) * sat_vp(Ta)
    p_skin = sat_vp(Ts)
    k_e = h_c / (P_ATM * 0.016)
    E_sw = L_V * k_e * (p_skin - p_air) / P_ATM

    # respiratory losses (classical PET/MEMI implementations)
    p_air_kPa = p_air / 1000.0
    C_res = 0.0014 * p.M * (34.0 - Ta)
    E_res = 0.0173 * p.M * (5.87 - p_air_kPa)
    Q_res = C_res + E_res

    F_core = p.M - Qc - Q_res
    F_skin = Qc - (C_env + R_env + E_sw)

    return np.array([F_core, F_skin], dtype=float)

def broyden_inverse(F, x0, args=(), tol=1e-4, max_iter=50) -> tuple[np.ndarray, bool, int]:
    """
    Inverse-Jacobian Broyden solver (Broyden, 1965).
    """
    x = np.array(x0, dtype=float)
    fx = F(x, *args)
    H = np.eye(x.size, dtype=float)

    for it in range(1, max_iter+1):
        if np.linalg.norm(fx, 2) < tol:
            return x, True, it

        s = -H.dot(fx)
        x_new = x + s
        fx_new = F(x_new, *args)

        y = fx_new - fx
        denom = float(y @ y)
        if denom == 0.0 or not np.isfinite(denom):
            break

        H = H + np.outer(s - H.dot(y), y) / denom
        x, fx = x_new, fx_new

    return x, False, max_iter

def solve_body_state_full(Ta: float, RH: float, v: float, Kg: float, Tg: float, p: Person2Node) -> tuple[float, float, bool]:
    x, ok, _ = broyden_inverse(memi2_residual, [37.0, 34.0], args=(Ta, RH, v, Kg, Tg, p))
    if not ok:
        return np.nan, np.nan, False
    Tc, Ts = float(x[0]), float(x[1])
    # physiologic plausibility bounds (document in methods)
    if not (30.0 <= Tc <= 43.0 and 20.0 <= Ts <= 45.0):
        return np.nan, np.nan, False
    return Tc, Ts, True

def pet_residual_full(Ta_ref: float, Ts_target: float, p: Person2Node) -> float:
    # indoor reference env: RH=50%, v=0.1 m/s, Kg=0, Tg=Ta_ref (common PET reference form)
    Tc_ref, Ts_ref, ok = solve_body_state_full(Ta_ref, 50.0, 0.1, 0.0, Ta_ref, p)
    if not ok:
        return np.nan
    return Ts_ref - Ts_target

def secant_pet_full(Ts_target: float, Ta_guess: float, p: Person2Node, tol: float = 0.3, max_iter: int = 30) -> float:
    if not np.isfinite(Ts_target):
        return np.nan
    a = Ta_guess - 10.0
    b = Ta_guess + 10.0
    fa = pet_residual_full(a, Ts_target, p)
    fb = pet_residual_full(b, Ts_target, p)

    for _ in range(max_iter):
        if not np.isfinite(fa) or not np.isfinite(fb):
            return np.nan
        if abs(fb) < tol:
            return float(b)
        denom = fb - fa
        if denom == 0.0:
            break
        b_new = b - fb * (b - a) / denom
        a, fa, b = b, fb, b_new
        fb = pet_residual_full(b, Ts_target, p)

    return float(b) if np.isfinite(b) else np.nan

# ----------------------------
# Tg QC + Station-specific handling
# ----------------------------

def apply_tg_qc(df: pd.DataFrame, tg_col: str = "BB_temp_Avg") -> pd.DataFrame:
    df = df.copy()
    df["Tg_flag_70_75"] = (df[tg_col] >= 70.0) & (df[tg_col] <= 75.0)
    df["Tg_QC"] = df[tg_col].where(df[tg_col] <= 75.0, np.nan)
    return df

def choose_tg_for_station(df: pd.DataFrame, station: int) -> pd.DataFrame:
    """
    Creates Tg_used and Tg_used_alt columns.

    - Stations != 4:
        Tg_used = Tg_QC (no recon available)
    - Station 4:
        Prefer QC-passed Tg; else fill with recon columns (two methods retained).
        This assumes your Station 4 file already includes:
        - Tg_rec_A
        - Tg_rec_D
    """
    df = apply_tg_qc(df, "BB_temp_Avg")

    if station != 4:
        df["Tg_used"] = df["Tg_QC"]
        df["Tg_used_alt"] = df["Tg_QC"]
        return df

    # Station 4: use QC where possible; otherwise recon
    if "Tg_rec_A" not in df.columns or "Tg_rec_D" not in df.columns:
        raise ValueError("Station 4 file missing Tg_rec_A and/or Tg_rec_D columns.")

    df["Tg_used_A"] = df["Tg_QC"].fillna(df["Tg_rec_A"])
    df["Tg_used_D"] = df["Tg_QC"].fillna(df["Tg_rec_D"])

    # Primary choice: method D (as per your prior workflow preference); keep A as alt
    df["Tg_used"] = df["Tg_used_D"]
    df["Tg_used_alt"] = df["Tg_used_A"]
    return df

# ----------------------------
# MAIN COMPUTATION
# ----------------------------

def compute_station_outputs(station: int, df: pd.DataFrame, svf_station: pd.Series) -> pd.DataFrame:
    """
    Adds:
    - sun azimuth
    - svf_sun (directional SVF at sun azimuth)
    - Kg_eff (shortwave attenuated)
    - MEMI outputs (Tc, Ts, MRT, PET)
    """
    out = df.copy()

    # solar azimuth: interpret timestamps as local naive; convert to UTC by localization if needed
    # Here we assume timestamps are local time already (common CR1000 field practice).
    # If your timestamps are UTC, remove the UTC conversion below.
    idx = out.index

    # Use pandas timezone localization if available
    try:
        idx_local = idx.tz_localize(SITE_TZ, nonexistent="shift_forward", ambiguous="NaT")
        idx_utc = idx_local.tz_convert("UTC")
    except Exception:
        # fallback: treat as UTC
        idx_utc = idx.tz_localize("UTC") if idx.tz is None else idx.tz_convert("UTC")

    az_list = [solar_azimuth_noaa(ts, SITE_LAT, SITE_LON) for ts in idx_utc]
    out["sun_azimuth_deg"] = az_list

    # directional SVF for sun sector
    svf_sun = []
    svf_col_used = []
    for az in out["sun_azimuth_deg"].values:
        col = sector_from_azimuth_deg(float(az))
        svf_col_used.append(col)
        svf_val = svf_station[col] if col in svf_station.index else np.nan
        svf_sun.append(float(svf_val) if np.isfinite(svf_val) else np.nan)

    out["svf_dir_col"] = svf_col_used
    out["SVF_sun_dir"] = svf_sun

    # hemispheric SVF as metadata if present
    if "SVF_total" in svf_station.index:
        out["SVF_total"] = float(svf_station["SVF_total"])
    else:
        out["SVF_total"] = np.nan

    # shortwave attenuation
    out["SlrW_eff"] = out["SlrW_Avg"] * out["SVF_sun_dir"]

    # MEMI per timestep
    Tc_list, Ts_list, MRT_list, PET_list, ok_list = [], [], [], [], []

    for _, r in out.iterrows():
        Ta = float(r["AirTC_Avg"]) if np.isfinite(r["AirTC_Avg"]) else np.nan
        RH = float(r["RH_Max"]) if np.isfinite(r["RH_Max"]) else np.nan
        v  = float(r["WS_ms_Avg"]) if np.isfinite(r["WS_ms_Avg"]) else np.nan
        Kg = float(r["SlrW_eff"]) if np.isfinite(r["SlrW_eff"]) else np.nan
        Tg = float(r["Tg_used"]) if np.isfinite(r["Tg_used"]) else np.nan

        if not np.all(np.isfinite([Ta, RH, v, Kg, Tg])):
            Tc_list.append(np.nan); Ts_list.append(np.nan); MRT_list.append(np.nan); PET_list.append(np.nan); ok_list.append(False)
            continue

        Tc, Ts, ok = solve_body_state_full(Ta, RH, v, Kg, Tg, PERSON)
        if not ok:
            Tc_list.append(np.nan); Ts_list.append(np.nan); MRT_list.append(np.nan); PET_list.append(np.nan); ok_list.append(False)
            continue

        MRT = mrt_from_globe(Tg, Ta, v)
        PET = secant_pet_full(Ts, Ta, PERSON)
        if (not np.isfinite(PET)) or PET < -20 or PET > 80:
            PET = np.nan

        Tc_list.append(Tc); Ts_list.append(Ts); MRT_list.append(MRT); PET_list.append(PET); ok_list.append(True)

    out["Tc_full"] = Tc_list
    out["Ts_full"] = Ts_list
    out["MRT"]     = MRT_list
    out["PET"]     = PET_list
    out["MEMI_ok"] = ok_list

    return out

def summarise_station_day(station: int, df_day: pd.DataFrame) -> dict:
    return {
        "station": station,
        "date": str(df_day.index.min().date()) if len(df_day) else "",
        "n_total": int(len(df_day)),
        "n_pet_valid": int(df_day["PET"].notna().sum()),
        "PET_mean": float(df_day["PET"].mean(skipna=True)) if df_day["PET"].notna().any() else np.nan,
        "PET_max": float(df_day["PET"].max(skipna=True)) if df_day["PET"].notna().any() else np.nan,
        "PET_p95": float(df_day["PET"].quantile(0.95)) if df_day["PET"].notna().any() else np.nan,
        "Ta_max": float(df_day["AirTC_Avg"].max(skipna=True)) if "AirTC_Avg" in df_day else np.nan,
        "MRT_max": float(df_day["MRT"].max(skipna=True)) if df_day["MRT"].notna().any() else np.nan,
        "SVF_total": float(df_day["SVF_total"].iloc[0]) if "SVF_total" in df_day.columns and len(df_day) else np.nan,
    }

def main():
    # Load station-level SVF
    svf_station = load_station_svf()

    all_rows = []
    summaries = []

    for st, fn in STATION_FILES.items():
        path = os.path.join(BASE_DIR, fn)
        df = read_station_file(path)
        df = coerce_station_types(df)
        df = choose_tg_for_station(df, st)

        # filter to hot days
        df_hot = pd.concat([df.loc[d] for d in HOT_DAYS if d in df.index.strftime("%Y-%m-%d")], axis=0)
        # alternative robust filter:
        df_hot = df[df.index.strftime("%Y-%m-%d").isin(HOT_DAYS)].copy()

        out = compute_station_outputs(st, df_hot, svf_station.loc[st])
        out["station"] = st
        all_rows.append(out)

        for d in HOT_DAYS:
            day_df = out[out.index.strftime("%Y-%m-%d") == d]
            summaries.append(summarise_station_day(st, day_df))

    out_all = pd.concat(all_rows).sort_values(["station"]).sort_index()

    summary_df = pd.DataFrame(summaries).sort_values(["date","station"])

    # Export workbook
    out_xlsx = "MEMI_PET_SVF_outputs_Stations1-8_Jan2018.xlsx"
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        out_all.to_excel(xw, sheet_name="timestep_outputs")
        summary_df.to_excel(xw, sheet_name="daily_summary", index=False)
        svf_station.to_excel(xw, sheet_name="SVF_station")

        # A minimal references sheet (expand if your thesis style requires full fields)
        refs = pd.DataFrame({
            "Reference (Harvard)": [
                "Broyden, C.G. (1965) ‘A class of methods for solving nonlinear simultaneous equations’, Mathematics of Computation, 19(92), pp. 577–593.",
                "Höppe, P. (1999) ‘The physiological equivalent temperature – a universal index for the biometeorological assessment of the thermal environment’, International Journal of Biometeorology, 43, pp. 71–75.",
                "Matzarakis, A., Mayer, H. and Iziomon, M.G. (1999) ‘Applications of a universal thermal index: physiological equivalent temperature’, International Journal of Biometeorology, 43, pp. 76–84.",
                "Bröde, P. et al. (2012) ‘Deriving the operational procedure for the universal thermal climate index (UTCI)’, International Journal of Biometeorology, 56, pp. 481–494. (Background for thermophysiological/radiative treatment conventions).",
                "ISO (1998) ISO 7726: Ergonomics of the thermal environment — Instruments for measuring physical quantities."
            ]
        })
        refs.to_excel(xw, sheet_name="References", index=False)

    # Export per-station CSVs
    for st in sorted(STATION_FILES.keys()):
        df_st = out_all[out_all["station"] == st].copy()
        df_st.to_csv(f"Station{st}_MEMI_PET_SVF_Jan2018.csv")

    print(f"Written: {out_xlsx}")
    print("Written: Station*_MEMI_PET_SVF_Jan2018.csv (per station)")

if __name__ == "__main__":
    main()
