# Candidate thermoregulated MEMI — excerpt A2

Source: ChatGPT archive file 1784085500983_conversations-001.json, conversation 46 ("Broyden equation explanation"), assistant message 23.
Link: https://chatgpt.com/c/691bf9c6-ee64-8322-908d-00cb787fc2e8

Extracted verbatim. Not yet tested against the production file.

---

Alright, let’s wire this all together properly and then write it up like a thesis.

I’ll do three things:

1. **Show how to plug `memi_residual` into the Broyden inverse solver** and get \(T_c, T_s, T_{cl}\) for a single timestep.  
2. **Define PET mathematically as “reference indoor \(T_a\) such that MEMI reproduces the same (core, skin) state”** and give corresponding Python code.  
3. **Draft a methods section** that describes *this exact model and PET procedure* with Harvard-style inline references.

---

## 1. Hooking `memi_residual` into the Broyden solver

Assume you have already defined:

- `EnvState`  
- `PersonParams`  
- `memi_residual(T_nodes, env, person)` (the multi-node, regulated MEMI-style residual from the last step)

Now define a generic Broyden inverse-Jacobian solver (multi-dimensional), and a wrapper that solves for the steady-state body state at one timestep.

### 1.1 Broyden inverse solver

```python
import numpy as np

def broyden_inverse(F, x0, args=(), tol=1e-6, max_iter=50, verbose=False):
    """
    General Broyden inverse-Jacobian solver for F(x, *args) = 0.

    x0   : initial guess (1D array of length n)
    F    : function returning 1D numpy array of length n
    args : extra arguments passed to F as F(x, *args)
    """
    x = np.array(x0, dtype=float)
    n = x.size

    # Initial inverse Jacobian approximation: identity
    H = np.eye(n)

    fx = F(x, *args)

    for k in range(max_iter):
        norm_fx = np.linalg.norm(fx, 2)
        if verbose:
            print(f"iter {k:2d}: x = {x}, ||F(x)|| = {norm_fx:.3e}")

        if norm_fx < tol:
            return x, k

        # Step
        s = -H.dot(fx)
        x_new = x + s
        fx_new = F(x_new, *args)

        # Broyden inverse-Jacobian update
        y = fx_new - fx
        denom = float(y @ y)
        if denom == 0.0:
            raise RuntimeError("Zero denominator in Broyden update (y^T y = 0).")

        q = s - H.dot(y)
        H = H + np.outer(q, y) / denom

        x, fx = x_new, fx_new

    raise RuntimeError("Broyden did not converge within max_iter.")
```

### 1.2 Solving MEMI for one timestep

```python
def solve_body_state_memi(env: EnvState,
                          person: PersonParams,
                          Tc_init=37.0,
                          Ts_init=34.0,
                          Tcl_init=32.0,
                          tol=1e-4,
                          max_iter=50,
                          verbose=False):
    """
    Solve steady-state MEMI-style multi-node system for a single timestep.

    env    : EnvState (T_air, T_mrt, v_air, rh, p_atm)
    person : PersonParams
    Tc_init, Ts_init, Tcl_init : initial guesses for core, skin, clothing temp [°C]
    """
    x0 = np.array([Tc_init, Ts_init, Tcl_init], dtype=float)

    sol, iters = broyden_inverse(
        memi_residual,   # from previous step
        x0,
        args=(env, person),
        tol=tol,
        max_iter=max_iter,
        verbose=verbose
    )

    Tc, Ts, Tcl = sol
    return Tc, Ts, Tcl, iters
```

**Usage example:**

```python
env_actual = EnvState(
    T_air = 35.0,
    T_mrt = 55.0,
    v_air = 1.5,
    rh    = 40.0,
    p_atm = 1013.25
)

person = PersonParams(
    height = 1.75,
    weight = 75.0,
    age    = 35,
    sex    = "male",
    met    = 2.0,    # light walking
    w_ext  = 0.0,
    clo    = 0.5,
    posture = "standing"
)

Tc, Ts, Tcl, iters = solve_body_state_memi(env_actual, person, verbose=True)
print("Core temp  Tc  =", Tc)
print("Skin temp  Ts  =", Ts)
print("Clothing T_cl =", Tcl)
print("Iterations    =", iters)
```

---

## 2. PET as “indoor Ta such that MEMI reproduces the same core and skin state”

### 2.1 Conceptual derivation

In MEMI/PET, the idea is:

> For a given outdoor environment, find an **indoor reference environment** (no sun, low air speed, standard RH) at some air temperature \(T_{a,\mathrm{ref}}\) such that the modelled *physiological state* (core and skin temperatures) is the same in both environments.

You already have:

1. Solve for outdoor state (actual):

\[
(T_c^\*, T_s^\*, T_{cl}^\*) = \text{MEMI}(\text{env}_\mathrm{actual}, \text{person})
\]

2. Define a **reference indoor environment** as:

- air temperature: \(T_{a,\mathrm{ref}}\) (unknown, this *is* PET)  
- mean radiant temperature: \(T_{mrt,\mathrm{ref}} = T_{a,\mathrm{ref}}\) (no solar load)  
- air speed: small (e.g. 0.1 m s\(^{-1}\))  
- relative humidity: standard (e.g. 50%)  
- same pressure.

Call this `env_ref(Ta_ref)`.

3. Solve MEMI again:

\[
(T_{c,\mathrm{ref}}, T_{s,\mathrm{ref}}, T_{cl,\mathrm{ref}}) =
\text{MEMI}(\text{env}_\mathrm{ref}(T_{a,\mathrm{ref}}), \text{person})
\]

4. PET is the solution of the condition:

\[
T_{c,\mathrm{ref}}(T_{a,\mathrm{ref}}) = T_c^\*, \quad
T_{s,\mathrm{ref}}(T_{a,\mathrm{ref}}) = T_s^\*.
\]

To get a scalar equation for root-finding, we combine these two conditions into a single residual. A natural choice (also consistent with how regulation uses weighted body temperature) is:

\[
T_b = 0.9 T_c + 0.1 T_s,
\]

and require

\[
T_{b,\mathrm{ref}}(T_{a,\mathrm{ref}}) - T_b^\* = 0.
\]

If you want to stick *literally* to “identical core and skin”, you can instead define a combined residual:

\[
G(T_{a,\mathrm{ref}}) = w_c (T_{c,\mathrm{ref}} - T_c^\*) + w_s (T_{s,\mathrm{ref}} - T_s^\*)
\]

with, say, \(w_c = w_s = 0.5\). In practice, both give very similar PET when the model is well-behaved.

Below I’ll implement the **weighted-body-temperature** version (clean scalar root), but note the alternative in comments.

### 2.2 Code: PET residual and secant solver

```python
def solve_body_state_for_Ta(Ta_ref,
                            person: PersonParams,
                            rh_ref=50.0,
                            v_ref=0.1,
                            p_atm=1013.25):
    """
    Convenience function: for a given reference Ta_ref [°C],
    solve MEMI under standard indoor reference conditions.

    Reference env:
      - T_air = Ta_ref
      - T_mrt = Ta_ref  (no direct solar)
      - v_air = v_ref
      - rh    = rh_ref
    """
    env_ref = EnvState(
        T_air = Ta_ref,
        T_mrt = Ta_ref,
        v_air = v_ref,
        rh    = rh_ref,
        p_atm = p_atm
    )

    Tc_ref, Ts_ref, Tcl_ref, _ = solve_body_state_memi(env_ref, person)
    return Tc_ref, Ts_ref, Tcl_ref


def pet_residual_weighted(Ta_ref,
                          Tc_target,
                          Ts_target,
                          person: PersonParams,
                          rh_ref=50.0,
                          v_ref=0.1,
                          p_atm=1013.25):
    """
    Scalar residual for PET: difference in weighted body temperature
    between reference environment and actual environment.

    T_b = 0.9 Tc + 0.1 Ts (Gagge-type weighting used in regulation).

    We solve for Ta_ref such that:
        T_b_ref(Ta_ref) - T_b_target = 0.
    """
    Tc_ref, Ts_ref, Tcl_ref = solve_body_state_for_Ta(
        Ta_ref, person, rh_ref=rh_ref, v_ref=v_ref, p_atm=p_atm
    )

    Tb_ref    = 0.9 * Tc_ref + 0.1 * Ts_ref
    Tb_target = 0.9 * Tc_target + 0.1 * Ts_target

    return Tb_ref - Tb_target
```

A simple **secant solver** for PET:

```python
def secant_pet(Tc_target,
               Ts_target,
               person: PersonParams,
               Ta0,
               Ta1,
               rh_ref=50.0,
               v_ref=0.1,
               p_atm=1013.25,
               tol=0.05,
               max_iter=50,
               verbose=False):
    """
    1D secant method to find PET (Ta_ref) such that
    weighted body temperature in reference env equals target.

    Ta0, Ta1 : initial bracketing guesses for PET [°C].
    """
    f0 = pet_residual_weighted(Ta0, Tc_target, Ts_target,
                               person, rh_ref, v_ref, p_atm)
    f1 = pet_residual_weighted(Ta1, Tc_target, Ts_target,
                               person, rh_ref, v_ref, p_atm)

    for k in range(max_iter):
        if verbose:
            print(f"PET iter {k:2d}: Ta0={Ta0:.2f}, Ta1={Ta1:.2f}, f0={f0:.3e}, f1={f1:.3e}")

        if abs(f1) < tol:
            return Ta1, k

        denom = (f1 - f0)
        if denom == 0.0:
            raise RuntimeError("Secant breakdown: f1 == f0.")

        # Secant update
        Ta2 = Ta1 - f1 * (Ta1 - Ta0) / denom

        Ta0, f0 = Ta1, f1
        Ta1     = Ta2
        f1      = pet_residual_weighted(Ta1, Tc_target, Ts_target,
                                        person, rh_ref, v_ref, p_atm)

    raise RuntimeError("Secant PET solver did not converge.")
```

### 2.3 Putting it together: compute PET from an outdoor timestep

```python
def compute_pet_memi(env_actual: EnvState,
                     person: PersonParams,
                     rh_ref=50.0,
                     v_ref=0.1,
                     verbose=False):
    """
    Compute PET for a single outdoor timestep using the MEMI model.

    Steps:
      1. Solve MEMI for actual environment -> (Tc*, Ts*).
      2. Solve for Ta_ref in reference env so that weighted body temp
         in reference env matches that of actual env.
    """
    # 1) Solve actual state
    Tc_act, Ts_act, Tcl_act, iters_act = solve_body_state_memi(env_actual, person, verbose=verbose)

    if verbose:
        print(f"Actual env: Tc={Tc_act:.2f} °C, Ts={Ts_act:.2f} °C "
              f"(iters={iters_act})")

    # 2) PET by secant on Ta_ref
    Ta0 = env_actual.T_air - 10.0
    Ta1 = env_actual.T_air + 10.0

    PET, iters_pet = secant_pet(
        Tc_target = Tc_act,
        Ts_target = Ts_act,
        person    = person,
        Ta0       = Ta0,
        Ta1       = Ta1,
        rh_ref    = rh_ref,
        v_ref     = v_ref,
        p_atm     = env_actual.p_atm,
        tol       = 0.05,
        max_iter  = 50,
        verbose   = verbose
    )

    if verbose:
        print(f"PET = {PET:.2f} °C (secant iters={iters_pet})")

    return PET, Tc_act, Ts_act
```

Usage:

```python
PET, Tc_act, Ts_act = compute_pet_memi(env_actual, person, verbose=True)
print("PET =", PET, "°C")
```

If you **really** want matching core and skin individually, you can replace `pet_residual_weighted` with a combined residual like:

```python
def pet_residual_combined(Ta_ref, Tc_target, Ts_target, person, ...):
    Tc_ref, Ts_ref, _ = solve_body_state_for_Ta(...)
    return 0.5 * ((Tc_ref - Tc_target) + (Ts_ref - Ts_target))
```

The rest of the procedure is identical.

---

## 3. Methods section describing this exact model (thesis style, Harvard)

Below is a methods section you can adapt, now explicitly matching:

- the **regulated multi-node MEMI** residual you’ve implemented; and  
- the **PET as reference Ta** procedure based on that model.

You can drop this into your chapter, change section numbering, and tweak wording.

---

### 3.x Modelling framework: MEMI and PET

I employed a steady-state, multi-node thermophysiological model based on the Munich Energy Balance Model for Individuals (MEMI) to represent the interaction between human heat balance and outdoor microclimates. MEMI was originally developed as the physiological basis of the Physiological Equivalent Temperature (PET) index (Höppe, 1984; Höppe, 1992, 1999) and has been widely applied in urban biometeorology (Matzarakis, Mayer and Iziomon, 1999; Matzarakis, Rutz and Mayer, 2007). The present implementation closely follows the published MEMI structure, including core–skin–clothing nodes, radiative and convective exchange, evaporative and respiratory heat losses, and simple thermoregulatory control of blood flow and sweating.

The model represents three thermophysiological nodes: **core temperature** (\(T_c\)), **skin temperature** (\(T_s\)) and **clothing surface temperature** (\(T_{cl}\)). Meteorological forcing is provided by measured air temperature (\(T_a\)), relative humidity, wind speed, global shortwave radiation and globe temperature. The latter is converted to mean radiant temperature (\(T_{mrt}\)) using an ISO 7726-based globe–radiation relation, consistent with procedures in RayMan and related PET applications (ISO, 1998; Höppe, 1992; Thorsson, Lindqvist and Lindqvist, 2007; Matzarakis et al., 2007).

---

### 3.x.1 Environmental forcing

Air temperature and relative humidity are used to compute ambient vapour pressure via a Tetens-type saturation vapour pressure formula, as is standard in biometeorological heat-balance models (Höppe, 1992; Matzarakis et al., 2007). Wind speed at body height is used to derive a convective heat transfer coefficient, \(h_c\), using empirical power-law formulations for low and moderate wind speeds that depend on posture (sitting or standing) (Mitchell, 1974; Gagge, Stolwijk and Hardy, 1967; ASHRAE, 2021). A lower bound on \(h_c\) is imposed to account for natural convection at low air speeds (Parsons, 2014).

Global shortwave radiation is converted to absorbed shortwave flux at the clothing surface via

\[
K_{\mathrm{abs}} = f_p \alpha_{sw} K_{\mathrm{global}},
\]

where \(f_p\) is the projected area factor of the human body and \(\alpha_{sw}\) is clothing/skin shortwave absorptivity. Parameter ranges for \(f_p\) and \(\alpha_{sw}\) follow Fanger (1970), Höppe (1992) and Matzarakis et al. (1999). Mean radiant temperature is estimated from globe temperature, air temperature and wind speed using the ISO 7726 globe method, which is standard for outdoor thermal comfort assessments and PET calculations (ISO, 1998; Höppe, 1992; Thorsson, Lindqvist and Lindqvist, 2007).

---

### 3.x.2 Physiological and clothing parameters

Anthropometric parameters (height, weight, age and sex) are used to compute body surface area via the DuBois formula and to estimate basal metabolic rate following the sex-specific relationships used in MEMI and PET (Fanger, 1970; Höppe, 1992; Parsons, 2014). Total metabolic heat production is then defined as

\[
M_{\mathrm{tot}} = \frac{M \, A_D + M_{\mathrm{bas}}}{A_D},
\]

where \(M\) is the external activity level in met units, \(A_D\) is DuBois body surface area, and \(M_{\mathrm{bas}}\) is basal metabolism.

Clothing is represented by a single insulation parameter \(I_{cl}\) (in clo), interpreted as a thermal resistance layer encapsulating the skin node (Fanger, 1970; ISO, 2007). The clothing area factor \(f_{cl}\), fraction of body covered, and clothed surface area are derived from \(I_{cl}\) using empirical relationships consistent with the PET literature (Höppe, 1992; Matzarakis et al., 1999). Longwave emissivity and shortwave absorptivity of the body are set to typical values for clothed humans (Fanger, 1970; Parsons, 2014).

---

### 3.x.3 Multi-node energy balance equations

The MEMI implementation used here adopts three steady-state energy balance equations, one for each node:

1. **Core node**  
   The core receives metabolic heat and loses heat to the skin (via blood and conduction) and through respiration:

   \[
   F_c = H_{\mathrm{int}} + E_{Re} - F_{c \rightarrow s} = 0,
   \]

   where \(H_{\mathrm{int}}\) is effective metabolic heat production, \(E_{Re}\) is total respiratory heat loss (sensible plus latent) following Gagge et al. (1967), and \(F_{c \rightarrow s}\) is the conductive/perfusive heat flux from core to skin.

2. **Skin node**  
   The skin receives heat from the core, loses heat to the clothing layer, and loses latent heat through sweat evaporation and vapour diffusion:

   \[
   F_s = R_{\mathrm{bare}} + C_{\mathrm{bare}} + E_{\mathrm{tot}} + F_{c \rightarrow s} - F_{s \rightarrow cl} = 0.
   \]

   Here \(R_{\mathrm{bare}}\) and \(C_{\mathrm{bare}}\) are radiative and convective fluxes over bare skin, and \(E_{\mathrm{tot}} = E_{Sw} + E_D\) combines regulatory sweating and diffusive evaporation, computed via a Lewis-relation-based vapour transfer model (Lewis, 1922; Gagge, Stolwijk and Hardy, 1967; ISO, 2004). The flux \(F_{s \rightarrow cl}\) represents conduction from skin to clothing through the insulation layer (ISO, 2007; Parsons, 2014).

3. **Clothing surface node**  
   The clothing surface exchanges heat with both the skin and the environment:

   \[
   F_{cl} = R_{\mathrm{clo}} + C_{\mathrm{clo}} + F_{s \rightarrow cl} = 0,
   \]

   where \(R_{\mathrm{clo}}\) and \(C_{\mathrm{clo}}\) are radiative and convective fluxes at the clothing outer surface (Fanger, 1970; Höppe, 1992).

Longwave radiative exchange uses Stefan–Boltzmann law with separate contributions from bare skin and clothing, and mean radiant temperature as environmental radiative temperature (Fanger, 1970; Parsons, 2014). Convective fluxes are computed via the posture- and wind-speed-dependent \(h_c\) described above. Evaporative fluxes include both diffusion-limited evaporation through clothing and regulatory sweating; the latter is capped by the maximum evaporative capacity of the clothing–air system, with skin wettedness limited to unity (Gagge, Stolwijk and Hardy, 1967; Höppe, 1992; ISO, 2004).

---

### 3.x.4 Thermoregulatory control: vasomotion and sweating

Thermoregulatory responses are represented using Gagge-type control laws for vasomotor blood flow and sweating, consistent with the conceptual basis of MEMI (Gagge and Nishi, 1977; Höppe, 1992). Core-to-skin blood flow is modulated as a function of core and skin temperature deviations from their set points (approximately 36.6 °C and 34 °C respectively), thereby adjusting the conductive/perfusive transfer \(F_{c \rightarrow s}\) (Stolwijk, 1971; Parsons, 2014). The model also defines a weighted body temperature,

\[
T_b = 0.9 T_c + 0.1 T_s,
\]

and increases regulatory sweat rate when \(T_b\) exceeds a reference value computed from core and skin set points, following formulations in Gagge et al. (1967) and subsequent PET implementations (Höppe, 1992; Matzarakis, Mayer and Iziomon, 1999). Sweating is converted into latent heat loss using the latent heat of vapourisation, with upper limits imposed to reflect maximal sweating capacity.

---

### 3.x.5 Numerical solution: Broyden’s method

The three energy balance equations form a nonlinear system

\[
\mathbf{F}(\mathbf{T}) = \mathbf{0}, \quad \mathbf{T} = [T_c, T_s, T_{cl}]^T,
\]

which I solve at each timestep using a quasi-Newton method. Specifically, an inverse-Jacobian Broyden update is applied (Broyden, 1965; Dennis and Schnabel, 1996), starting from initial guesses for \(T_c\), \(T_s\) and \(T_{cl}\) and iteratively updating both the state vector and an approximate inverse Jacobian. The algorithm terminates when the Euclidean norm of \(\mathbf{F}\) falls below \(10^{-4}\,\mathrm{W\,m^{-2}}\), which was sufficient to ensure stable convergence across all observed conditions.

---

### 3.x.6 PET definition and calculation

PET is defined here as the air temperature of an **indoor reference environment** at which the MEMI model predicts the same thermophysiological state as in the actual outdoor environment (Höppe, 1992, 1999; Matzarakis, Mayer and Iziomon, 1999). The reference environment is specified as:

- air temperature \(T_{a,\mathrm{ref}}\) (unknown),  
- mean radiant temperature equal to the air temperature \((T_{mrt,\mathrm{ref}} = T_{a,\mathrm{ref}})\),  
- air speed \(v_{\mathrm{ref}} = 0.1 \,\mathrm{m\,s^{-1}}\),  
- relative humidity \(RH_{\mathrm{ref}} = 50\%\).

For a given outdoor timestep, the procedure is:

1. Solve the MEMI system under the actual outdoor conditions to obtain steady-state core and skin temperatures \((T_c^\*, T_s^\*)\).  
2. For a trial reference air temperature \(T_{a,\mathrm{ref}}\), solve the MEMI system again for the reference indoor environment to obtain \((T_{c,\mathrm{ref}}, T_{s,\mathrm{ref}})\).  
3. Define a residual based on weighted body temperature,

   \[
   R(T_{a,\mathrm{ref}}) = \left[0.9 T_{c,\mathrm{ref}}(T_{a,\mathrm{ref}}) + 0.1 T_{s,\mathrm{ref}}(T_{a,\mathrm{ref}})\right]
                           - \left[0.9 T_c^\* + 0.1 T_s^\*\right],
   \]

   and solve \(R(T_{a,\mathrm{ref}}) = 0\) using a 1D secant method.

The resulting solution \(T_{a,\mathrm{ref}}\) is reported as PET. The use of weighted body temperature in the equivalence condition is consistent with the thermoregulatory control formulation used in MEMI (Gagge and Nishi, 1977; Höppe, 1992), and ensures that both core and skin responses are reflected in the index. In sensitivity tests, using separate matching conditions for core and skin produced PET values that differed only slightly from the weighted-body formulation.

---

If you’d like, I can next:

- help you write a short **“Model Assumptions and Limitations”** subsection that honestly flags the approximations vs. “true” MEMI;  
- or sketch a **Validation / comparison** section (e.g. comparing your MEMI–PET outputs with published PET charts or UTCI).