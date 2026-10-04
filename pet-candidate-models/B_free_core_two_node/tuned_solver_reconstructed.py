import numpy as np
import pandas as pd

SIGMA = 5.670374419e-8
P_ATM = 101325.0
L_V   = 2.42e6

def sat_vp(T):
    """Saturation vapour pressure [Pa] at T [°C]."""
    return 610.78 * np.exp(17.27 * T / (T + 237.3))

def mrt_from_globe(Tg, Ta, v, D_globe=0.15, emiss_globe=0.95):
    """MRT from 0.15 m black globe (ISO / Hoeppe-style)."""
    TgK = Tg + 273.15
    TaK = Ta + 273.15
    v_eff = max(v, 0.1)
    factor = 1.1e8 * v_eff**0.6 / (emiss_globe * D_globe**0.4)
    TmrtK4 = TgK**4 + factor * (TgK - TaK)
    return TmrtK4**0.25 - 273.15

class Person2Node:
    def __init__(self, M=70.0, emiss=0.97, alpha_sw=0.50, f_p=0.70, k_core_skin=8.0):
        self.M = M
        self.emiss = emiss
        self.alpha_sw = alpha_sw
        self.f_p = f_p
        self.k_core_skin = k_core_skin

def memi2_residual(x, Ta, RH, v, Kg, Tg, person):
    """2-node MEMI residuals: F_core, F_skin = 0 (Hoeppe, 1984; Hoeppe, 1999)."""
    Tc, Ts = x
    Tmrt = mrt_from_globe(Tg, Ta, v)
    TsK   = Ts   + 273.15
    TmrtK = Tmrt + 273.15

    Qc = person.k_core_skin * (Tc - Ts)

    v_eff = max(v, 0.1)
    h_c = 8.3 * v_eff**0.6

    C_env = h_c * (Ts - Ta)

    K_abs  = person.f_p * person.alpha_sw * Kg
    LW_net = person.emiss * SIGMA * (TsK**4 - TmrtK**4)
    R_env  = LW_net - K_abs

    p_air  = RH / 100.0 * sat_vp(Ta)
    p_skin = sat_vp(Ts)
    k_e = h_c / (P_ATM * 0.016)
    E_sw = L_V * k_e * (p_skin - p_air) / P_ATM

    p_air_kPa = p_air / 1000.0
    C_res = 0.0014 * person.M * (34.0 - Ta)
    E_res = 0.0173 * person.M * (5.87 - p_air_kPa)
    Q_res = C_res + E_res

    F_core = person.M - Qc - Q_res
    F_skin = Qc - (C_env + R_env + E_sw)

    return np.array([F_core, F_skin])

def broyden_inverse(F, x0, args=(), tol=1e-4, max_iter=40):
    x = np.array(x0, dtype=float)
    fx = F(x, *args)
    n = x.size
    H = np.eye(n)
    for _ in range(max_iter):
        if np.linalg.norm(fx, 2) < tol:
            return x
        s = -H.dot(fx)
        x_new = x + s
        fx_new = F(x_new, *args)
        y = fx_new - fx
        denom = float(y @ y)
        if denom == 0.0:
            break
        H = H + np.outer(s - H.dot(y), y) / denom
        x, fx = x_new, fx_new
    return x

def solve_body_state_full(Ta, RH, v, Kg, Tg, person,
                           Tc_init=37.0, Ts_init=34.0, tol=1e-4, max_iter=40):
    x0 = np.array([Tc_init, Ts_init])
    x = broyden_inverse(memi2_residual, x0, args=(Ta, RH, v, Kg, Tg, person),
                         tol=tol, max_iter=max_iter)
    Tc, Ts = x
    if not (30.0 <= Tc <= 43.0 and 20.0 <= Ts <= 45.0):
        return np.nan, np.nan
    return Tc, Ts

def pet_residual_full(Ta_ref, Ts_target, person):
    Tc_ref, Ts_ref = solve_body_state_full(Ta_ref, 50.0, 0.1, 0.0, Ta_ref, person)
    if np.isnan(Ts_ref):
        return np.nan
    return Ts_ref - Ts_target

def secant_pet_full(Ts_target, Ta_act, person, tol=0.3, max_iter=30):
    if np.isnan(Ts_target):
        return np.nan
    a = Ta_act - 10.0
    b = Ta_act + 10.0
    fa = pet_residual_full(a, Ts_target, person)
    fb = pet_residual_full(b, Ts_target, person)
    for _ in range(max_iter):
        if np.isnan(fa) or np.isnan(fb):
            return np.nan
        if abs(fb) < tol:
            return b
        denom = fb - fa
        if denom == 0.0:
            break
        b_new = b - fb * (b - a) / denom
        a, fa, b, fb = b, fb, b_new, pet_residual_full(b_new, Ts_target, person)
    return b

def compute_pet_row(Ta, RH, v, Kg, Tg, person):
    if np.any(np.isnan([Ta, RH, v, Kg, Tg])):
        return np.nan, np.nan, np.nan, np.nan
    Tc, Ts = solve_body_state_full(Ta, RH, v, Kg, Tg, person)
    if np.isnan(Tc) or np.isnan(Ts):
        return np.nan, np.nan, np.nan, np.nan
    MRT = mrt_from_globe(Tg, Ta, v)
    PET = secant_pet_full(Ts, Ta, person)
    if PET < -20 or PET > 80:
        PET = np.nan
    return Tc, Ts, MRT, PET
