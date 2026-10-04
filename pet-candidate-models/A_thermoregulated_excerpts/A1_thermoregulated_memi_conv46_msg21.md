# Candidate thermoregulated MEMI — excerpt A1

Source: ChatGPT archive file 1784085500983_conversations-001.json, conversation 46 ("Broyden equation explanation"), assistant message 21.
Link: https://chatgpt.com/c/691bf9c6-ee64-8322-908d-00cb787fc2e8

Extracted verbatim. Not yet tested against the production file.

---

Short answer: I can get you *very close* to a full MEMI-style multi-node model with thermo-regulation, but I can’t honestly claim a mathematically exact replica of Höppe’s original dissertation, because the full set of equations is not publicly visible and sits behind paywalls. What I *can* do is:

- Use the published MEMI energy balance structure (core–skin–clothing, M+W+R+C+ED+ERe+ESw+S=0) citeturn2view0turn1search16  
- Use the widely documented Gagge-type control laws for vasomotion and sweating, which MEMI builds on citeturn2view0turn4view0  
- Follow the same 3-node layout and regulatory ideas used in modern open implementations of PET (e.g. Walther & Goestchel’s correction, as implemented in `pythermalcomfort`) citeturn6view1turn7view0  

So below is a MEMI-style residual function that:

- Has **three nodes**: core, skin, clothing  
- Includes **blood-flow regulation (vasomotion)** and **sweating control**  
- Implements **radiation, convection, evaporation, and respiratory** terms explicitly  
- Is structured so you can plug it straight into your existing Broyden solver

I’ll keep parameters explicit so you can later match them to specific citations in your thesis.

---

## 1. MEMI-style data structures

```python
from dataclasses import dataclass
import numpy as np

@dataclass
class EnvState:
    """Environmental / meteorological state for MEMI/PET."""
    T_air: float          # Air temperature [°C]
    T_mrt: float          # Mean radiant temperature [°C]
    v_air: float          # Wind speed at body level [m/s]
    rh: float             # Relative humidity [%]
    p_atm: float = 1013.25  # Atmospheric pressure [hPa]

@dataclass
class PersonParams:
    """Person and activity parameters for MEMI."""
    height: float = 1.75      # [m]
    weight: float = 75.0      # [kg]
    age: int = 35             # [years]
    sex: str = "male"         # "male" or "female"
    met: float = 1.4          # Activity level [met]
    w_ext: float = 0.0        # External work [met]
    clo: float = 0.9          # Clothing insulation [clo]
    posture: str = "sitting"  # "sitting", "standing", "standing_forced"
```

---

## 2. Helper functions: saturation vapour pressure, body surface area

These are standard and not specific to MEMI.

```python
def p_sat_Tetens(T_c):
    """
    Saturation vapour pressure over water [hPa] using a Tetens-type formula.
    Valid ~0–50 °C. (e.g. as used in many biometeorological models)
    """
    return 6.1078 * np.exp(17.2694 * T_c / (T_c + 237.3))


def body_surface_area_dubois(weight, height):
    """
    DuBois body surface area [m²].
    """
    return 0.20247 * (height**0.725) * (weight**0.425)
```

---

## 3. Thermoregulatory control: vasomotion & sweating

Here we implement a **Gagge-type 2-node controller** (on which MEMI is based) citeturn2view0turn4view0  

The constants (set-point temperatures, gains, limits) are taken from the thermoregulation literature and are conceptually the same type as in MEMI / PET implementations citeturn1search16turn6view1.

```python
def vasomotion_blood_flow(T_core, T_skin):
    """
    Blood flow and mass distribution between core and skin.
    Approximates the MEMI / Gagge-style vasomotor control.

    Returns
    -------
    m_blood : float
        Blood perfusion [kg m⁻² h⁻¹]
    alpha   : float
        Mass fraction at skin temperature [-]
    """
    # Set-point temperatures (core and skin) [°C]
    Tc_set = 36.6
    Tsk_set = 34.0

    # Deviation signals (only warm core and cool skin drive vasodilation)
    sig_core = max(0.0, T_core - Tc_set)
    sig_skin = max(0.0, Tsk_set - T_skin)

    # Baseline & gain for blood flow (Gagge-type, as used in PET codes)
    m_blood = (6.3 + 75.0 * sig_core) / (1.0 + 0.5 * sig_skin)  # [L m⁻² h⁻¹ approx]

    # Upper limit for blood flow
    m_blood = min(m_blood, 90.0)

    # Convert to kg m⁻² h⁻¹ with blood density ~1 kg/L
    m_blood_kg = m_blood  # density ≈ 1 kg/L

    # Fraction of body mass at skin temperature (alpha)
    # (empirical fit used in PET-type models)
    alpha = 0.0417737 + 0.7451833 / (m_blood + 0.585417)

    return m_blood_kg, alpha


def sweating_rate(T_core, T_skin):
    """
    Regulatory sweating as a function of weighted body temperature.

    Returns
    -------
    m_rsw : float
        Regulatory sweat rate [g m⁻² h⁻¹].
    """
    Tc_set = 36.6
    Tsk_set = 34.0

    # Reference body temperature as weighted core/skin
    T_body_set = 0.9 * Tc_set + 0.1 * Tsk_set

    # Current weighted body temperature
    T_body = 0.9 * T_core + 0.1 * T_skin

    sig_body = max(0.0, T_body - T_body_set)

    # Gagge-type proportional control (widely cited)
    m_rsw = 304.94 * sig_body   # [g m⁻² h⁻¹]

    # Upper limit (heavy sweating)
    m_rsw = min(m_rsw, 500.0)

    return m_rsw
```

---

## 4. MEMI 3-node residual function with regulation

This is the heart of what you asked for: a **multi-node MEMI-style residual** with regulation. The structure mirrors published descriptions of MEMI and its PET applications: core, skin and clothing balances plus full heat-budget terms (M, W, R, C, ED, ERe, ESw) citeturn2view0turn1search16turn6view1  

```python
def memi_residual(T_nodes, env: EnvState, person: PersonParams):
    """
    MEMI-style 3-node steady-state energy balance with thermoregulation.

    Parameters
    ----------
    T_nodes : array-like
        [T_core, T_skin, T_clo] in °C.
    env : EnvState
        Environmental state (Ta, Tmrt, v, RH, p_atm).
    person : PersonParams
        Anthropometry, clothing, activity etc.

    Returns
    -------
    res : ndarray, shape (3,)
        Residuals of energy balance for core, skin, clothing [W m⁻²].
        At solution, all three should be ~0 for steady state.
    """
    T_core, T_skin, T_clo = T_nodes

    # --- CONSTANTS & BASIC PHYSICS ---------------------------------------
    sigma = 5.67e-8         # Stefan–Boltzmann [W m⁻² K⁻⁴]
    eps_skin = 0.99         # emissivity skin [-]
    eps_clo = 0.95          # emissivity clothing [-]
    L_v = 2.42e6            # latent heat of evaporation [J kg⁻¹]
    c_blood = 3640.0        # specific heat blood [J kg⁻¹ K⁻¹]
    Lewis = 1.67            # Lewis relation [K hPa⁻¹] ~ 16.7e-1
    i_m = 0.38              # Woodcock vapour transfer ratio [-]

    Ta = env.T_air
    Tmrt = env.T_mrt
    v = max(env.v_air, 0.01)  # avoid 0 in v^0.67
    rh = env.rh
    p_atm = env.p_atm

    # Person
    h = person.height
    w = person.weight
    age = person.age
    sex = person.sex.lower()
    met = person.met
    w_ext = person.w_ext
    clo = person.clo
    posture = person.posture

    # --- BODY SURFACE AREA ----------------------------------------------
    A_D = body_surface_area_dubois(w, h)  # [m²]

    # --- METABOLIC HEAT PRODUCTION & BASE METABOLISM --------------------
    # Convert met units to W/m²
    met_Wm2 = met * 58.2                 # 1 met ≈ 58.2 W/m²
    w_ext_Wm2 = w_ext * 58.2

    # Sex-specific basal metabolism as in PET / MEMI implementations
    if sex == "male":
        M_bas = 3.45 * w**0.75 * (1.0 + 0.004*(30 - age) +
                                  0.01*(h*100.0 / (w**(1.0/3.0)) - 43.4))
    else:
        M_bas = 3.19 * w**0.75 * (1.0 + 0.004*(30 - age) +
                                  0.018*(h*100.0 / (w**(1.0/3.0)) - 42.1))

    M_tot = (met_Wm2 * A_D + M_bas) / A_D      # total internal heat production [W m⁻²]
    W_ext = w_ext_Wm2                           # external mechanical work [W m⁻²]
    H_int = M_tot * (1.0 - W_ext / max(M_tot, 1e-6))  # effective metabolic heat [W m⁻²]

    # --- CLOTHING GEOMETRY & AREAS --------------------------------------
    # Clothing area factor f_cl (Burton relation)
    f_cl = 1.0 + 0.31 * clo

    # Fraction of body covered by clothing (empirical)
    if clo >= 2.0:
        y = 1.0
    elif clo > 0.6:
        y = (h - 0.2)/h
    elif clo > 0.3:
        y = 0.5
    elif clo > 0.0:
        y = 0.1
    else:
        y = 0.0

    # Fractional covered area (polynomial approx used in PET-type models)
    f_a_cl = (173.51*clo - 2.36 - 100.76*clo**2 + 19.28*clo**3) / 100.0
    f_a_cl = min(max(f_a_cl, 0.0), 1.0)

    A_clo = A_D * (f_a_cl + (f_cl - 1.0))  # clothed surface
    f_eff = 0.696 if posture == "standing" else 0.725  # effective radiative factor
    A_r_eff = A_D * f_eff                   # effective radiative area

    # --- CONVECTION COEFFICIENT -----------------------------------------
    if posture == "sitting":
        h_c = 2.67 + 6.5 * v**0.67
    elif posture == "standing":
        h_c = 2.26 + 7.42 * v**0.67
    else:  # standing forced convection
        h_c = 8.6 * v**0.513

    # Natural convection floor at low pressure (Hoeppe-style)
    h_cc = 3.0 * (p_atm/1013.25)**0.53
    h_c = max(h_c, h_cc)
    h_c *= (p_atm/1013.25)**0.55

    # --- RESPIRATORY HEAT LOSSES (E_Re = C_res + Q_res) -----------------
    # Expired air temperature
    T_exp = 0.47 * Ta + 21.0
    # Pulmonary ventilation (per area)
    D_vent = H_int * 1.44e-6          # [m³ s⁻¹ m⁻²] approx
    # Sensible respiratory loss
    C_res = 1010.0 * (Ta - T_exp) * D_vent    # [W m⁻²]
    # Latent respiratory loss
    vpa = rh/100.0 * p_sat_Tetens(Ta) / 100.0 # [hPa]
    vp_exp = p_sat_Tetens(T_exp) / 100.0      # [hPa]
    Q_res = 0.623 * L_v / p_atm * (vpa - vp_exp) * D_vent
    E_re = C_res + Q_res

    # --- THERMOREGULATION: BLOOD FLOW & SWEAT ---------------------------
    m_blood, alpha = vasomotion_blood_flow(T_core, T_skin)
    T_body = alpha * T_skin + (1.0 - alpha) * T_core
    m_rsw = sweating_rate(T_core, T_skin)   # [g m⁻² h⁻¹]

    # Evaporative loss due to sweat
    E_sw = (L_v / 1000.0) * (m_rsw / 3600.0)   # [W m⁻²] (g→kg and h→s)

    # --- DIFFUSIVE EVAPORATION THROUGH SKIN / CLOTHING (E_D) -----------
    # Clothing resistance to sensible heat
    R_cl = clo / 6.45                # [m² K W⁻¹]

    # Saturation vapour pressure at skin temperature
    p_vs_skin = p_sat_Tetens(T_skin) / 100.0  # [hPa]

    # Convective mass transfer coefficient via Lewis relation
    h_e_air = h_c * Lewis         # [W m⁻² hPa⁻¹ approx]
    # Effective vapour transfer efficiency of clothing
    f_e_cl = 1.0 / (1.0 + 0.92 * h_c * R_cl)
    # Max possible skin evaporation (diffusion-limited)
    E_max = h_e_air * f_e_cl * (p_vs_skin - vpa)  # [W m⁻²]

    if E_max <= 0.0:
        E_max = 1e-3

    # Skin wettedness
    w_skin = E_sw / E_max
    if w_skin > 1.0:
        w_skin = 1.0
        # cap E_sw at E_max if necessary
        E_sw = min(E_sw, E_max)

    # Clothing vapour resistance (Woodcock)
    R_e_cl = (1.0 / (f_cl * h_c) + R_cl) / (Lewis * i_m)
    E_diff = (1.0 - w_skin) * (p_vs_skin - vpa) / R_e_cl  # [W m⁻²]

    # Total evaporative loss from skin surface
    E_tot = -(E_diff + E_sw)   # negative sign = loss from body

    # --- RADIATIVE EXCHANGE ---------------------------------------------
    # Convert to Kelvin
    Tk_air = Ta + 273.15
    Tk_mrt = Tmrt + 273.15
    Tk_skin = T_skin + 273.15
    Tk_clo = T_clo + 273.15

    # Bare skin area (fraction not covered)
    A_bare = A_D * (1.0 - f_a_cl)

    # Radiation from bare skin
    R_bare = (A_r_eff * (1.0 - f_a_cl) * eps_skin * sigma *
              (Tk_mrt**4 - Tk_skin**4) / A_D)

    # Radiation from clothing
    R_clo = (f_eff * A_clo * eps_clo * sigma *
             (Tk_mrt**4 - Tk_clo**4) / A_D)

    R_tot = R_bare + R_clo   # [W m⁻²]

    # --- CONVECTIVE EXCHANGE --------------------------------------------
    C_bare = h_c * (Ta - T_skin) * A_bare / A_D
    C_clo = h_c * (Ta - T_clo) * A_clo / A_D
    C_tot = C_bare + C_clo   # [W m⁻²]

    # --- CONDUCTION CORE–SKIN & SKIN–CLOTHING --------------------------
    # Clothing geometry: internal/external radii (cylindrical shell)
    if y > 0 and clo > 0:
        r2 = A_D * (f_cl - 1.0 + f_a_cl) / (2.0 * np.pi * h * y)
        r1 = f_a_cl * A_D / (2.0 * np.pi * h * y)
        di = max(r2 - r1, 1e-4)
        # Equivalent thermal conductance for tissues/clothing
        h_t_cl = (2.0 * np.pi * h * y * di) / (R_cl * np.log(r2 / r1) * A_clo)
    else:
        h_t_cl = 5.28  # fallback conductance [W m⁻² K⁻¹] (matches dry transfer constant)

    # Heat transfer core↔skin via blood and static conduction
    F_cs = (m_blood / 3600.0 * c_blood + 5.28) * (T_core - T_skin)
    # Heat transfer skin↔clothing
    F_sc = h_t_cl * (T_skin - T_clo)

    # --- MEMI-STYLE BALANCES -------------------------------------------
    # 1) Core node: metabolic heat – blood transfer to skin – storage (=0 at steady state)
    res_core = H_int + E_re - F_cs         # [W m⁻²]

    # 2) Skin node: radiation + convection + evaporation + blood inflow – conduction to clothing
    res_skin = R_bare + C_bare + E_tot + F_cs - F_sc

    # 3) Clothing node: radiation + convection + conduction in from skin
    res_clo = R_clo + C_clo + F_sc

    return np.array([res_core, res_skin, res_clo], dtype=float)
```

---

## 5. How this maps to the “original MEMI” concepts

In terms of **structure**, this does what you asked:

- **Multi-node**  
  - Core node (`res_core`)  
  - Skin node (`res_skin`)  
  - Clothing node (`res_clo`)  

- **MEMI energy balance terms** (using the notation in Höppe citeturn2view0turn1search16):  
  - \(M\): `M_tot` (metabolic heat)  
  - \(W\): `W_ext` (external work, subtracted in `H_int`)  
  - \(R\): `R_tot = R_bare + R_clo`  
  - \(C\): `C_tot = C_bare + C_clo`  
  - \(E_{Re}\): `E_re = C_res + Q_res`  
  - \(E_{Sw}\): `E_sw` (regulatory sweating)  
  - \(E_D\): `E_diff` (diffusion through skin/clothing)  
  - \(S\): set to **zero** here because we are solving for steady-state node temps; if you wanted transient MEMI, you’d add node heat capacities and write ODEs for \(dT/dt\) instead.

- **Thermo-regulation**  
  - **Vasomotor control** through `vasomotion_blood_flow`, modifying `m_blood` and thus `F_cs`, following Gagge-style control with core and skin set points citeturn4view0turn6view1  
  - **Sweating control** via `sweating_rate`, again Gagge-type body-temperature control feeding into `E_sw` citeturn4view0turn6view1  

- **Radiation, convection, evaporation & clothing**  
  - Radiative exchange uses Stefan–Boltzmann with separate bare-skin and clothed areas, consistent with MEMI and PET descriptions citeturn1search16turn1search19  
  - Convective exchange uses empirical correlations vs. wind speed and posture, as in contemporary PET implementations citeturn6view1turn3search18  
  - Vapour diffusion through clothing uses a Woodcock-type resistance and Lewis relation, which is exactly the conceptual framework described for MEMI-based PET citeturn6view1turn1search7  

Because the *full original* MEMI equations reside primarily in Höppe’s 1984 dissertation and are not fully exposed in open sources, I cannot prove that every coefficient and branching logic here is identical to that document. What this code represents is:

> A **3-node MEMI-style steady-state model with regulation** that is compatible with the way PET is described and implemented in modern literature, and that follows the same physical and physiological assumptions.

If you like, next step I can:

- Show how to hook this `memi_residual` into your existing **Broyden inverse-Jacobian solver**; then  
- Re-derive PET as “indoor Ta such that MEMI gives identical core and skin temperatures” for your thesis, and then  
- Draft the **methods section text** that explains this precise model with full inline references in Harvard style.