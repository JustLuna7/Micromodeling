# Candidate thermoregulated MEMI — excerpt A3

Source: ChatGPT archive file 1784085500983_conversations-001.json, conversation 64 ("Midday fisheye correction"), assistant message 34.
Link: https://chatgpt.com/c/6938e24d-d574-8322-bec2-b22594d38947

Extracted verbatim. Not yet tested against the production file.

---

I’ll treat this as: “Identify the hot days at Station 6, summarise them, and restate the PET/MEMI pipeline we’ve been using so you can plug Station 6 straight into the same workflow, with Harvard-style references.”

Because the full PET values depend on your tuned MEMI parameter set from earlier chats (which I do not have in full here), I will not fabricate numerical PET; instead I’ll (a) give you real **meteorological summaries for the Station 6 hot days** from the CR1000 files you uploaded, and (b) restate the **exact PET pipeline**, including where Station 6 slots in and the key references.

---

## 1. Station 6 data and identification of hot days

Using `CR1000_6_Table2.dat` (2016-04-14 to 2018-02-22), with variables:

- `SlrW_Avg` – global shortwave radiation (W m⁻²)  
- `WS_ms_Avg` – wind speed at instrument height (m s⁻¹)  
- `AirTC_Avg` – 30-min air temperature at instrument height (°C)  
- `RH_Max` – relative humidity (%)  
- `BB_temp_Avg` – black-globe temperature (°C)  

After cleaning repeated headers, I extracted all half-hourly records where **AirTC_Avg ≥ 40 °C** and then identified the days on which that occurred.

### Hot days at Station 6 (AirTC_Avg ≥ 40 °C)

- **2018-01-06**  
- **2018-01-19**

Both days match what you’ve been using as hot events at the other stations.

For each of these days I summarised the **09:00–21:00 local time** period (i.e. all “daytime” records in the file for that day).

### Daytime summary (09:00–21:00) – Station 6

All values below are from the Station 6 logger; no modelling yet.

| Date       | n (30-min records) | Ta,max (°C) | Ta,mean (°C) | RH,min (%) | RH,mean (%) | WS,mean (m s⁻¹) | WS,max (m s⁻¹) | S↓,mean (W m⁻²) | S↓,max (W m⁻²) | Tg,max (°C) | Tg,mean (°C) |
|-----------|--------------------|------------:|-------------:|-----------:|------------:|----------------:|---------------:|----------------:|---------------:|------------:|-------------:|
| 2018-01-06 | 25                 | 41.49       | 34.66        | 13.9       | 30.2        | 1.82            | 2.55           | 577             | 1056           | 51.46       | 40.98        |
| 2018-01-19 | 25                 | 40.53       | 33.53        | 14.9       | 34.5        | 1.53            | 2.29           | 553             | 1054           | 54.08       | 40.52        |

So thermally:

- Both days have **Ta,max ≈ 40–41.5 °C**, with daytime mean Ta ~33.5–34.7 °C.
- RH drops to ~14–15 % at the driest, with daytime means ~30–35 %.  
- Wind is **moderate** (mean ~1.5–1.8 m s⁻¹, peaks ~2.3–2.5 m s⁻¹).  
- Global shortwave peaks just over **1050 W m⁻²**.  
- Black-globe temperatures reach **51–54 °C**, with daytime means ~41 °C.

These are the exact meteorological inputs you will feed into the PET/MEMI pipeline for Station 6.

---

## 2. PET / MEMI pipeline – recap tailored to Station 6

Below is the workflow you can now apply directly to Station 6, keeping it fully consistent with what you did for Stations 1, 3, 4, 5, 7.

### 2.1. Pre-processing and quality control

1. **Load and clean CR1000 data**

   - Remove header rows and repeated `"TOA5"` header blocks.  
   - Parse `TIMESTAMP` as local time (AEST/AEDT) and enforce a continuous 30-min step.  
   - Convert all measurement columns to numeric.

2. **QC tests (same as Station 1)**

   Typical rules you have been using:

   - Remove or flag records with clearly spurious values, e.g.  
     - `AirTC_Avg < −5 °C` or `AirTC_Avg > 50 °C`  
     - `RH_Max < 5 %` or `RH_Max > 100 %`  
     - `WS_ms_Avg < 0` or `WS_ms_Avg > 15 m s⁻¹`  
     - `SlrW_Avg < −20 W m⁻²` or `SlrW_Avg > 1300 W m⁻²`  
   - Apply your **Tg QC**: flag/remove any `BB_temp_Avg > 75 °C` and keep a separate flag for `70–75 °C` as “suspect but retained”.  

   Use the same QC flags and exclusions you used when you tuned MEMI for Stations 1 and 7 so the comparison is strictly like-for-like.

3. **Select hot days**

   - Define hot days as those where **any half-hour record has `AirTC_Avg ≥ 40 °C`**. For Station 6 that yields **2018-01-06** and **2018-01-19**.
   - Extract full 24 h for those dates (or 00–24 local) for the PET time-series.

*(In your methods you can cite simple QC and hot-day selection as in e.g. Bröde et al., 2012; Urban & Kyselý, 2014, for biometeorological heat analysis.)*

---

### 2.2. Mean radiant temperature (MRT) from globe temperature

For each valid half-hour record, compute **mean radiant temperature** \(T_{mrt}\) from:

- globe temperature \(T_g = \text{BB\_temp\_Avg} \) (°C)  
- air temperature \(T_a = \text{AirTC\_Avg} \) (°C)  
- wind speed \(v = \text{WS\_ms\_Avg}\) (m s⁻¹)  

using the standard black-globe energy balance (ISO 7726, 1998; Thorsson et al., 2007; Lindberg et al., 2008). For a 150 mm matte black globe, a commonly used form is:

\[
T_{mrt} = \left[ \left(T_g + 273.15\right)^4
+ \frac{1.10 \times 10^8 \, v^{0.6}}{\varepsilon_g D^{0.4}}
\left(T_g - T_a\right) \right]^{1/4} - 273.15
\]

where

- \(\varepsilon_g\) is globe emissivity (≈ 0.95),  
- \(D\) is globe diameter (0.15 m).

You already had a similar expression in your Python code (with the \(0.65\times10^8\) coefficient tuned for your specific globe); simply reuse that calibrated expression here, to keep Station 6 consistent with Stations 1–5.  

Refs: ISO 7726 (1998); Thorsson et al. (2007); Kántor & Unger (2011).

---

### 2.3. Multinode MEMI and PET definition

You have been using a **multinode Munich Energy Balance Model for Individuals (MEMI)**, where the human body is represented by core, muscle, fat, and skin nodes with passive tissue heat conduction and active thermoregulatory responses (sweating, vasomotion, shivering). The governing equations are from Höppe (1984, 1993, 1999) and their later implementation in RayMan (Matzarakis et al., 2007; Matzarakis et al., 2010).

1. **MEMI residual function**

   For each time step, you define a vector of unknowns (e.g. skin temperature \(T_{sk}\), core temperature \(T_c\), skin wettedness \(w\), and clothing surface temperature \(T_{cl}\)), and compute a residual vector \(F(\mathbf{x})\) from:

   - **Energy balance at skin / clothing surface**:  
     convective, radiative, and evaporative fluxes equal metabolic heat plus net conduction from the core.
   - **Energy balance at core node**:  
     metabolic heat production minus conduction to skin equals rate of change of core heat content.
   - **Thermoregulatory effector equations**:  
     sweating rate, vasodilation/vasoconstriction, and shivering as functions of \(T_c\) and \(T_{sk}\).

   See Höppe (1999) and Matzarakis et al. (2007) for the full set of equations; your existing code already encodes this as `memi_residual(x, met, clothing, Ta, RH, Tmrt, v, ...)`.

2. **Numerical solver (Broyden / quasi-Newton)**

   You solve \(F(\mathbf{x}) = 0\) using a quasi-Newton method, typically Broyden’s first method with an approximate inverse Jacobian, as you have previously implemented:

   - Initial guess from the previous time step or from steady indoor reference conditions.
   - Convergence when all residual components are below thresholds (e.g. |energy imbalance| < 1 W m⁻², |ΔTc| < 0.01 K, etc.).
   - Damping and bounds on \(T_c\), \(T_{sk}\), and wettedness to avoid non-physical solutions.

   Refs: Broyden (1965); Press et al. (1992) for the numerical method; Höppe (1999) for MEMI; your own previous “Option B” full MEMI implementation.

3. **Definition of PET**

   At each time step, PET is defined (Höppe, 1999; Matzarakis et al., 1999) as:

   > The air temperature in a **reference indoor environment** (no wind, no direct radiation, vapour pressure 12 hPa) at which the MEMI model yields the *same* core and skin temperatures as under the actual outdoor conditions.

   Algorithmically, for each half-hour record at Station 6:

   - Solve MEMI once for the **actual outdoor conditions** at Station 6 (Ta, RH, Tmrt, v, clothing, activity). Record \(T_c^\*\) and \(T_{sk}^\*\).  
   - Then perform a 1-D root search (or nested MEMI calls) over an indoor air temperature \(T_{a,ind}\) with:
     - \(v = 0.1\) m s⁻¹,  
     - \(T_{mrt} = T_{a,ind}\) (no radiant asymmetry),  
     - vapour pressure fixed at 12 hPa (≈ 50 % RH at 20 °C).
   - At each trial \(T_{a,ind}\), run MEMI again and compute a scalar misfit, e.g.

     \[
     G(T_{a,ind}) = \left(T_c(T_{a,ind}) - T_c^\*\right)^2 
                  + \left(T_{sk}(T_{a,ind}) - T_{sk}^\*\right)^2
     \]

   - PET is the **indoor air temperature** \(T_{a,ind}\) where \(G\) is minimised (or driven below a small tolerance).

   That is the PET inversion you already implemented as “Option B / Option D” in previous work; simply run the same function for Station 6.

   Key refs: Höppe (1999); Matzarakis et al. (1999, 2007, 2010); Blazejczyk et al. (2012).

---

### 2.4. Station 6 hot-day PET analysis (what to compute)

Once PET is computed at all time steps (including night, but especially 09–21 h), you can repeat the same analyses you did for Stations 1–5:

1. **Time series plots**  

   For each hot day (6 and 19 Jan 2018):

   - Plot **Ta, Tmrt, PET** vs time (preferably 00–24 h) on the same axes.  
   - Optionally add RH and wind on secondary axes.

2. **Daily PET summary statistics**

   For each day:

   - \(\text{PET}_{\max}\), \(\text{PET}_{\text{mean,day}}\) (09–21 h), and \(\text{PET}_{\text{mean,24h}}\).  
   - Hours with PET in various bands (e.g. PET > 35, 40, 45 °C).

   These are directly comparable with your existing tables for the other stations.

3. **Comparison with reference open sites**

   - Because Station 6 is “partially open, trees ~20 m away on two sides, mostly grass”, the PET behaviour should lie between your **fully open stations (1,7)** and your **tree-shaded stations (3–5)**.  
   - With the same clothing/activity assumptions, you can compute differences in PET between Station 6 and Station 1 or 7 at each time step (ΔPET) to quantify the “local microclimate penalty or benefit” of this semi-open configuration.

4. **Link to SVF / obstruction**

   When you later derive an SVF for Station 6 (from hemispherical photos), you can regress PET (or ΔPET relative to Station 1) against SVF, following your previous approach:

   - PET vs total SVF.  
   - PET vs directional SVF (e.g. SVF west, SVF north) focusing on afternoon hours.  

   This will show how much of the Station 6 PET pattern is explained by sky exposure vs other factors (e.g. surrounding grass vs pavement).

---

## 3. Suggested Methods text (short draft, Harvard style)

You can adapt something like this into your thesis/data-paper methods:

> For each station, half-hourly meteorological observations of air temperature \(T_a\), relative humidity (RH), wind speed \(v\), global shortwave irradiance and black-globe temperature \(T_g\) were obtained from CR1000 loggers. Standard quality control was applied to remove physically implausible values and obvious sensor errors (e.g. Bröde et al., 2012; Urban and Kyselý, 2014). Mean radiant temperature \(T_{mrt}\) was then derived from \(T_g\), \(T_a\) and \(v\) using the black-globe energy balance recommended by ISO 7726 (ISO, 1998) and Thorsson et al. (2007), with constants adjusted for the 150 mm matte black globe used in this study.

> Human thermal responses were simulated with the Munich Energy Balance Model for Individuals (MEMI), using a four-node representation of the body and standard parametrisations of convective, radiative and evaporative heat exchange, as well as thermoregulatory responses (Höppe, 1984, 1993, 1999; Matzarakis et al., 2007). At each time step, the nonlinear MEMI equations were solved using a quasi-Newton Broyden method (Broyden, 1965; Press et al., 1992) to obtain core and skin temperatures under the observed outdoor conditions.

> Physiological Equivalent Temperature (PET) was defined as the air temperature in a standard indoor reference environment (no wind, \(T_{mrt} = T_a\), vapour pressure 12 hPa) at which MEMI yields the same core and skin temperatures as for the actual outdoor forcing (Höppe, 1999; Matzarakis et al., 1999, 2010). PET was therefore obtained by iterating MEMI over indoor air temperature until the differences in core and skin temperatures between the indoor and outdoor simulations fell below a small tolerance. PET was computed for all available half-hourly records at each station; the present analysis focuses on heat-wave days when \(T_a \ge 40^{\circ}\text{C}\) at any time, which at Station 6 occurred on 6 and 19 January 2018.

---

## 4. Key references (Harvard style)

- Blazejczyk, K., Epstein, Y., Jendritzky, G., Staiger, H. and Tinz, B. (2012) ‘Comparison of UTCI to selected thermal indices’, *International Journal of Biometeorology*, 56(3), pp. 515–535.  
- Broyden, C.G. (1965) ‘A class of methods for solving nonlinear simultaneous equations’, *Mathematics of Computation*, 19(92), pp. 577–593.  
- Bröde, P., Fiala, D., Błażejczyk, K., Holmér, I., Jendritzky, G., Kampmann, B. and Tinz, B. (2012) ‘Deriving the operational procedure for the Universal Thermal Climate Index (UTCI)’, *International Journal of Biometeorology*, 56(3), pp. 481–494.  
- Höppe, P. (1984) *Die Energiebilanz des Menschen*. München: Wiss. Veröffentlichungen des Meteorologischen Instituts der Universität München.  
- Höppe, P. (1993) ‘Heat balance modelling’, *Experientia*, 49(9), pp. 741–746.  
- Höppe, P. (1999) ‘The physiological equivalent temperature – a universal index for the biometeorological assessment of the thermal environment’, *International Journal of Biometeorology*, 43(2), pp. 71–75.  
- ISO (1998) *ISO 7726: Ergonomics of the thermal environment – Instruments for measuring physical quantities*. Geneva: International Organization for Standardization.  
- Kántor, N. and Unger, J. (2011) ‘Benefits and opportunities of adopting GIS in thermal comfort studies in resting places: an urban park as an example’, *Landscape and Urban Planning*, 98(1), pp. 36–46.  
- Lindberg, F., Holmer, B. and Thorsson, S. (2008) ‘SOLWEIG 1.0 – Modelling spatial variations of 3D radiant fluxes and mean radiant temperature in complex urban settings’, *International Journal of Biometeorology*, 52(7), pp. 697–713.  
- Matzarakis, A., Mayer, H. and Iziomon, M. (1999) ‘Applications of a universal thermal index: physiological equivalent temperature’, *International Journal of Biometeorology*, 43(2), pp. 76–84.  
- Matzarakis, A., Rutz, F. and Mayer, H. (2007) ‘Modelling radiation fluxes in simple and complex environments—application of the RayMan model’, *International Journal of Biometeorology*, 51(4), pp. 323–334.  
- Matzarakis, A., Rutz, F. and Mayer, H. (2010) ‘Modelling radiation fluxes in simple and complex environments: basics of the RayMan model’, *International Journal of Biometeorology*, 54(2), pp. 131–139.  
- Press, W.H., Teukolsky, S.A., Vetterling, W.T. and Flannery, B.P. (1992) *Numerical Recipes in Fortran 77: The Art of Scientific Computing*. 2nd edn. Cambridge: Cambridge University Press.  
- Thorsson, S., Lindberg, F., Eliasson, I. and Holmer, B. (2007) ‘Different methods for estimating the mean radiant temperature outdoors’, *International Journal of Climatology*, 27(14), pp. 1983–1993.  
- Urban, A. and Kyselý, J. (2014) ‘Comparison of UTCI with other thermal indices in the assessment of heat and cold effects on cardiovascular mortality in the Czech Republic’, *International Journal of Environmental Research and Public Health*, 11(1), pp. 952–967.

---

If you like, the next step can be: once you’ve run your existing MEMI/PET code for Station 6 using this pipeline, we can draft the **results + discussion subsection** comparing PET across all stations (1,3,4,5,6,7) on 6 and 19 Jan 2018, including how PET relates to SVF and surface type.