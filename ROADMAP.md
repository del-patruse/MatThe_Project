# Roadmap — Perfect Plasticity with Hill's Anisotropic Yield Criterion

**Exam 2 · Introduction to Theory of Materials · TU Dortmund, Institute of Mechanics**

Goal: implement, from zero, a small-strain rate-independent **perfectly plastic** material model with **Hill's orthotropic yield criterion**, driven by a strain-driven constitutive driver, under **equi-biaxial** loading (ramp + cyclic), with a **numerical tangent**, and verify it by recovering **von Mises** as a special case.

> Deadlines: code → **August 28th, 9:00 am** to tim.furlan@tu-dortmund.de · slides → day before the exam, 11:59 pm · exam → **September 3rd**.

---

## Phase 0 — Groundwork (before writing any model equation)

### 0.1 Fix your notation and storage convention

Everything downstream depends on one consistent convention. Decide it now and write it down at the top of your derivation notes:

- Small strains, additive split: $\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}^e + \boldsymbol{\varepsilon}^p$.
- Stress and strain as symmetric 3×3 tensors internally, plus a **6-component matrix (Voigt-like) representation** for the yield function, exactly as the assignment defines it:

$$\boldsymbol{\sigma}^m_{\text{dev}} = \begin{bmatrix} \sigma_{\text{dev}\,11} & \sigma_{\text{dev}\,22} & \sigma_{\text{dev}\,33} & \sigma_{\text{dev}\,12} & \sigma_{\text{dev}\,23} & \sigma_{\text{dev}\,13} \end{bmatrix}^T$$

- **Watch the ordering!** The assignment uses (11, 22, 33, **12, 23, 13**) — this is *not* the standard Voigt order (which puts 23 first among the shears). Every mapping tensor↔vector must respect the assignment's order, and the $\mathbf{G}$ matrix below is written for exactly this order.
- **No factor-2 on shear strains here.** The quadratic form uses *stress* components directly; the factors of 2 for shear are already built into $\mathbf{G}$ (the $2L, 2M, 2N$ entries). Do not additionally double shear components when converting.

### 0.2 Understand the constitutive driver from the tutorials

The driver imposes a mixed strain/stress history: some strain components are *prescribed* over time, the remaining components are found by a global Newton iteration enforcing the corresponding stress components to be zero (stress-free directions). Key ingredients you will reuse:

- `partition()` — decides which of the 6 components are strain-controlled vs. stress-controlled. **You must modify it** so that *both* $\varepsilon_{11}$ *and* $\varepsilon_{22}$ are prescribed (equi-biaxial), with all other stress components driven to zero.
- `numericalTangent()` — forward-difference approximation of $\partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}$. The assignment **explicitly forbids** the analytical algorithmic tangent — you use the numerical one.
- Consequence for your code structure (per the assignment hint): implement **two separate functions** — one returning the stress for a given strain and history, one returning the (numerical) tangent by repeatedly calling the stress function with perturbed strains.

### 0.3 Collect the material parameters (Table 1 of the assignment)

| Parameter | $E$ | $\nu$ | $\sigma_y^{11}$ | $\sigma_y^{22}$ | $\sigma_y^{33}$ | $\sigma_y^{12}$ | $\sigma_y^{23}$ | $\sigma_y^{13}$ |
|---|---|---|---|---|---|---|---|---|
| Value | 200 000 | 0.3 | 580 | 510 | 460 | 590 | 530 | 550 |

(Stress units MPa; consistent with the lecture's plate-with-hole example, which uses the same yield stresses with $G = 76\,923$ MPa, $K = 166\,667$ MPa — a useful cross-check for your elastic constants.)

---

## Phase 1 — Theory: derive the model on paper (Task 1)

Do this completely on paper **before** coding. This derivation is also the core of your presentation and oral exam.

### 1.1 Elasticity

Linear isotropic Hooke elasticity (lecture §5.4):

$$\boldsymbol{\sigma} = \mathbb{E}^e : \boldsymbol{\varepsilon}^e, \qquad \mathbb{E}^e = 2G\,\mathbb{I}^{sym}_{dev} + K\,\boldsymbol{I}\otimes\boldsymbol{I}$$

with $G = \dfrac{E}{2(1+\nu)} \approx 76\,923$ MPa and $K = \dfrac{E}{3(1-2\nu)} \approx 166\,667$ MPa. Equivalently split:

$$\boldsymbol{\sigma}_{dev} = 2G\,\boldsymbol{\varepsilon}^e_{dev}, \qquad \sigma_m = K\,\varepsilon^e_{vol}$$

### 1.2 Hill's yield function

Tensor form (lecture eq. 6.51 / assignment eq. 2):

$$\Phi = \boldsymbol{\sigma}_{dev} : \mathbb{G} : \boldsymbol{\sigma}_{dev} - 1 \le 0$$

Matrix form (assignment eq. 3–5), with the 6×6 structural matrix in the assignment's component order:

$$\mathbf{G} = \begin{bmatrix}
F+G & -F & -G & 0 & 0 & 0\\
-F & F+H & -H & 0 & 0 & 0\\
-G & -H & G+H & 0 & 0 & 0\\
0 & 0 & 0 & 2L & 0 & 0\\
0 & 0 & 0 & 0 & 2M & 0\\
0 & 0 & 0 & 0 & 0 & 2N
\end{bmatrix}$$

Expanded (classical Hill 1950 format, assignment eq. 6):

$$\Phi = F[\sigma_{dev11}-\sigma_{dev22}]^2 + G[\sigma_{dev11}-\sigma_{dev33}]^2 + H[\sigma_{dev22}-\sigma_{dev33}]^2 + 2L\sigma_{dev12}^2 + 2M\sigma_{dev23}^2 + 2N\sigma_{dev13}^2 - 1$$

**Important observations to derive and state:**

- Differences of deviatoric components equal differences of total components: $\sigma_{dev11}-\sigma_{dev22} = \sigma_{11}-\sigma_{22}$, etc. So $\Phi$ is automatically **mean-stress independent** ($\mathbb{G}:\boldsymbol{I} = \boldsymbol{0}$, lecture eq. 6.69) — evaluating it with $\boldsymbol{\sigma}_{dev}$ or $\boldsymbol{\sigma}$ gives the same value. Still follow the assignment and feed it the deviator.
- $\mathbb{G}$ has major and minor symmetry and is constant (perfect plasticity, no hardening: the yield surface never changes).

### 1.3 Parameter identification (assignment eqs. 7–10 / lecture eq. 6.80–6.81)

From uniaxial yield stresses along the orthotropy axes and shear yield stresses:

$$F = \tfrac{1}{2}\Big[\tfrac{1}{(\sigma_y^{11})^2} + \tfrac{1}{(\sigma_y^{22})^2} - \tfrac{1}{(\sigma_y^{33})^2}\Big], \quad
G = \tfrac{1}{2}\Big[\tfrac{1}{(\sigma_y^{11})^2} + \tfrac{1}{(\sigma_y^{33})^2} - \tfrac{1}{(\sigma_y^{22})^2}\Big], \quad
H = \tfrac{1}{2}\Big[\tfrac{1}{(\sigma_y^{22})^2} + \tfrac{1}{(\sigma_y^{33})^2} - \tfrac{1}{(\sigma_y^{11})^2}\Big]$$

$$L = \tfrac{1}{2(\sigma_y^{12})^2}, \qquad M = \tfrac{1}{2(\sigma_y^{23})^2}, \qquad N = \tfrac{1}{2(\sigma_y^{13})^2}$$

Derive these yourself by evaluating $\Phi = 0$ for uniaxial stress along each axis (three equations, lecture eq. 6.79) and pure shear in each plane — this derivation is a likely oral-exam question.

**Convexity check** (lecture eq. 6.90–6.94, good presentation material): $L, M, N > 0$ automatically; $F, G, H$ must satisfy $F+G+H > 0$ and $FG + FH + GH \ge 0$. Verify numerically with the given parameters. The eigenvalues of $\mathbf{G}_{nn}$ are $\mu_1 = 0$ (hydrostatic direction) and $\mu_{2,3} = F{+}G{+}H \pm \sqrt{F^2{+}G^2{+}H^2 - (FG{+}FH{+}GH)}$.

**Von Mises as special case** (assignment eqs. 11–12): with
$$F = G = H = \frac{1}{2\sigma_y^2}, \qquad L = M = N = \frac{3}{2\sigma_y^2}$$
Hill's criterion collapses to $\Phi = \frac{3}{2}\frac{\|\boldsymbol{\sigma}_{dev}\|^2}{\sigma_y^2} - 1$, i.e. the von Mises criterion $\sigma_e \le \sigma_y$ in "squared" form. Derive this equivalence on paper (substitute into eq. 6 and compare with $\sigma_e = \sqrt{\tfrac{3}{2}}\|\boldsymbol{\sigma}_{dev}\|$). This becomes your main verification tool later. Note the two forms are *equivalent as yield surfaces* ($\Phi = 0$ coincides) but the yield-function *values* and flow-rule magnitudes differ ($\lambda$ scales differently) — the stress–strain response must nevertheless coincide, because the product $\lambda\,\boldsymbol{\nu}$ is what enters the update.

### 1.4 Constitutive framework — perfect plasticity (lecture §5.2)

Assemble the full set of continuum equations:

1. **Additive split:** $\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}^e + \boldsymbol{\varepsilon}^p$; internal variable: $\boldsymbol{\varepsilon}^p$ only (no hardening variables).
2. **Elastic law:** $\boldsymbol{\sigma} = \mathbb{E}^e : [\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^p]$.
3. **Yield condition:** $\Phi(\boldsymbol{\sigma}) \le 0$ with $\Phi$ from §1.2 — note $\Phi$ depends on stress only (perfect plasticity: no $\kappa$, no $\boldsymbol{\alpha}$).
4. **Associative flow rule** (normality, from the principle of maximum dissipation, lecture §5.2/eq. 6.100):
$$\dot{\boldsymbol{\varepsilon}}^p = \lambda\,\boldsymbol{\nu}, \qquad \boldsymbol{\nu} = \frac{\partial\Phi}{\partial\boldsymbol{\sigma}} = 2\,\mathbb{G} : \boldsymbol{\sigma}_{dev}$$
   In matrix form: $\boldsymbol{\nu}^m = 2\,\mathbf{G}\,\boldsymbol{\sigma}^m_{dev}$. **Careful:** the 6-vector $\boldsymbol{\nu}^m$ holds *tensor* components of $\boldsymbol{\nu}$; when reassembling the symmetric 3×3 flow direction, the off-diagonal entries appear twice ($\nu_{12} = \nu_{21}$, etc.) but are *not* doubled in value.
5. **Plastic incompressibility:** $\text{tr}\,\dot{\boldsymbol{\varepsilon}}^p = \lambda\, \boldsymbol{I} : 2\mathbb{G}:\boldsymbol{\sigma}_{dev} = 0$ since $\mathbb{G}:\boldsymbol{I} = \boldsymbol{0}$ — plastic flow is purely deviatoric (lecture eq. 6.102). State and later verify this.
6. **Karush–Kuhn–Tucker (loading/unloading) conditions:**
$$\lambda \ge 0, \qquad \Phi \le 0, \qquad \lambda\,\Phi = 0$$
7. **Dissipation:** $\mathcal{D} = \boldsymbol{\sigma} : \dot{\boldsymbol{\varepsilon}}^p = \lambda\,\boldsymbol{\sigma} : 2\mathbb{G} : \boldsymbol{\sigma} = 2\lambda\,[\Phi + 1] = 2\lambda \ge 0$ at yield ($\Phi = 0$) — thermodynamic consistency guaranteed. A one-line derivation worth showing.

### 1.5 Time discretization — Backward Euler return mapping (lecture §5.3, §6.6.2)

Discretize the flow rule with the implicit (Backward) Euler rule over $[t_n, t_{n+1}]$, $\mu := \Delta t\,\lambda_{n+1}$ (the incremental consistency parameter, often written $\Delta\lambda$):

$$\boldsymbol{\varepsilon}^p_{n+1} = \boldsymbol{\varepsilon}^p_n + \mu\,\boldsymbol{\nu}(\boldsymbol{\sigma}_{n+1})$$

**Elastic predictor (trial state).** Freeze plastic flow:

$$\boldsymbol{\sigma}^{tr} = \mathbb{E}^e : [\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}^p_n], \qquad \Phi^{tr} = \Phi(\boldsymbol{\sigma}^{tr})$$

- If $\Phi^{tr} \le 0$: step is elastic. Accept $\boldsymbol{\sigma}_{n+1} = \boldsymbol{\sigma}^{tr}$, $\boldsymbol{\varepsilon}^p_{n+1} = \boldsymbol{\varepsilon}^p_n$. Done.
- If $\Phi^{tr} > 0$: **plastic corrector** required.

**Plastic corrector.** Enforce simultaneously (lecture eq. 6.103):

$$\boldsymbol{R}_\sigma(\boldsymbol{\sigma}, \mu) = [\mathbb{E}^e]^{-1} : [\boldsymbol{\sigma} - \boldsymbol{\sigma}^{tr}] + \mu\,\boldsymbol{\nu}(\boldsymbol{\sigma}) = \boldsymbol{0}$$
$$R_\mu(\boldsymbol{\sigma}) = \Phi(\boldsymbol{\sigma}) = 0$$

(all quantities at $t_{n+1}$; the first equation is just the discretized flow rule rewritten in stresses).

**Key point of the whole assignment (Hint 1):** for von Mises perfect plasticity the corrector reduces to *radial return* with a closed-form $\mu$ (lecture eq. 5.38: $\mu = [\sigma_e^{tr} - \sigma_y]/3G$). For Hill's criterion the flow direction $\boldsymbol{\nu} = 2\mathbb{G}:\boldsymbol{\sigma}_{dev}$ is **not parallel** to $\boldsymbol{\sigma}^{tr}_{dev}$ (the return is no longer radial, because $\mathbb{G} \ne \tfrac{3}{2\sigma_y^2}\mathbb{I}^{sym}_{dev}$), so $\Phi_{n+1} = 0$ **cannot be solved explicitly for $\mu$**. You must iterate: **local Newton–Raphson**.

### 1.6 The local Newton–Raphson scheme — derive it explicitly

Two equivalent formulations. Pick **one**, derive it fully, and be able to defend the choice.

**Option A — coupled Newton in $(\boldsymbol{\sigma}, \mu)$ (the lecture's generic scheme, eq. 6.104–6.106).** Unknowns $X = (\boldsymbol{\sigma}^m, \mu)$, 7 scalar unknowns. Per iteration solve the linearized system

$$\begin{bmatrix} [\mathbb{E}^e]^{-1} + \mu\,\mathbf{N} & \boldsymbol{\nu} \\ \boldsymbol{\nu}^T & 0 \end{bmatrix} \begin{bmatrix} d\boldsymbol{\sigma} \\ d\mu \end{bmatrix} = -\begin{bmatrix} \boldsymbol{R}_\sigma \\ R_\mu \end{bmatrix}$$

where $\mathbf{N} = \partial\boldsymbol{\nu}/\partial\boldsymbol{\sigma} = 2\mathbb{G}$ is **constant** — a major simplification specific to quadratic yield functions. The lecture (eq. 6.106) shows the block-decomposition with $\mathbb{E}^e_a = [[\mathbb{E}^e]^{-1} + 2\mu\,\mathbb{G}]^{-1}$: compute $r_\sigma = \mathbb{E}^e_a : \boldsymbol{R}_\sigma$, $\boldsymbol{F} = \mathbb{E}^e_a : \boldsymbol{\nu}$, $h_a = \boldsymbol{\nu} : \boldsymbol{F}$, then $d\mu = [R_\mu - \boldsymbol{\nu} : r_\sigma]/h_a$ and $d\boldsymbol{\sigma} = -r_\sigma - \boldsymbol{F}\,d\mu$. In 6×6 matrix representation this is direct linear algebra.

**Option B — scalar Newton in $\mu$ alone (fits Hint 1's wording most directly).** Because $\mathbb{E}^e$ is isotropic and plastic flow is deviatoric, the discrete stress update can be condensed: from $\boldsymbol{R}_\sigma = 0$,

$$\boldsymbol{\sigma}^m_{dev}(\mu) = \left[\mathbf{I} + 2\mu\, \mathbf{D}\,\mathbf{G}\right]^{-1} \boldsymbol{\sigma}^{m,tr}_{dev}, \qquad \sigma_m = \sigma_m^{tr}$$

where $\mathbf{D}$ carries the deviatoric elastic moduli ($2G$ on normal components; the shear rows get the appropriate $2G$ factors — derive carefully with your component convention, this is where the shear factor-of-2 bookkeeping bites). Then solve the single scalar equation $f(\mu) := \Phi(\boldsymbol{\sigma}_{dev}(\mu)) = 0$ by Newton:

$$\mu^{(k+1)} = \mu^{(k)} - \frac{f(\mu^{(k)})}{f'(\mu^{(k)})}, \qquad f'(\mu) = \boldsymbol{\nu} : \frac{d\boldsymbol{\sigma}_{dev}}{d\mu}$$

with $\frac{d\boldsymbol{\sigma}_{dev}}{d\mu} = -[\mathbf{I} + 2\mu\,\mathbf{D}\,\mathbf{G}]^{-1}\, 2\,\mathbf{D}\,\mathbf{G}\; \boldsymbol{\sigma}_{dev}(\mu)$. Start from $\mu^{(0)} = 0$ (i.e. from the trial state, where $f(0) = \Phi^{tr} > 0$).

Properties to note for either option: $f(\mu)$ is monotonically decreasing and convex in $\mu$ (quadratic yield function, positive semi-definite $\mathbb{G}$), so Newton from $\mu = 0$ converges monotonically — mention this robustness argument in the presentation.

**Convergence criterion:** $|\Phi| < \text{tol}$ (e.g. $10^{-10}$, plus a residual-norm check on $\boldsymbol{R}_\sigma$ if Option A) and a max-iteration guard (~25) with a hard error if exceeded.

**After convergence:** update $\boldsymbol{\varepsilon}^p_{n+1} = \boldsymbol{\varepsilon}^p_n + \mu\,\boldsymbol{\nu}(\boldsymbol{\sigma}_{n+1})$ and store it as the new history.

### 1.7 Tangent operator

Per the assignment: **do not** derive/implement the algorithmic (ATS) tangent analytically. The tangent for the global driver iteration is obtained by **forward differences**:

$$\left[\mathbb{E}_a\right]_{ij} \approx \frac{\sigma_i(\boldsymbol{\varepsilon} + \delta\,\boldsymbol{e}_j,\ \text{history}_n) - \sigma_i(\boldsymbol{\varepsilon},\ \text{history}_n)}{\delta}$$

with perturbation $\delta \approx 10^{-7}$–$10^{-8}$ (relative to strain magnitude; too small → round-off, too large → truncation error). Crucial detail: every perturbed evaluation must start from the **same converged history** $\boldsymbol{\varepsilon}^p_n$ (do **not** let the perturbed calls overwrite the stored plastic strain). This is exactly why stress computation and tangent computation must be two separate, side-effect-free functions.

---

## Phase 2 — The algorithm, end to end (write this as a flowchart before coding)

Per global time step $t_n \to t_{n+1}$, per global Newton iteration of the driver:

1. Driver proposes a full strain tensor $\boldsymbol{\varepsilon}_{n+1}$ (prescribed components from the load curve, free components from the driver's Newton update).
2. **Material routine** (stateless w.r.t. its inputs: takes $\boldsymbol{\varepsilon}_{n+1}$, $\boldsymbol{\varepsilon}^p_n$; returns $\boldsymbol{\sigma}_{n+1}$, candidate $\boldsymbol{\varepsilon}^p_{n+1}$):
   1. trial stress; evaluate $\Phi^{tr}$,
   2. elastic if $\Phi^{tr} \le 0$; otherwise local Newton for $\mu$ (and $\boldsymbol{\sigma}$),
   3. return stress + updated internal variables.
3. **Numerical tangent**: 6 extra calls of the material routine with perturbed strain, same history.
4. Driver checks the residual (stress components that must vanish); iterates 1–3 until converged.
5. **Only after global convergence**: commit $\boldsymbol{\varepsilon}^p_{n+1}$ to history, log all output quantities, advance to the next step.

Quantities to log every step (you need them all for Phase 5/6): $t$, full $\boldsymbol{\varepsilon}$, $\boldsymbol{\varepsilon}^p$, $\boldsymbol{\sigma}$, $\Phi$, $\mu$, local-Newton iteration count, elastic/plastic flag.

---

## Phase 3 — Driver adaptation: equi-biaxial load case

The prescribed history (assignment eq. 13):

$$\bar{\varepsilon}_{n+1} := \varepsilon_{11}(t_{n+1}) = \varepsilon_{22}(t_{n+1})$$

- Modify `partition()` so components 11 **and** 22 are strain-driven; components 33, 12, 23, 13 are stress-free ($\sigma_{33} = \sigma_{12} = \sigma_{23} = \sigma_{13} = 0$) and solved for by the driver's global Newton.
- **Load case 1 — ramp (Fig. 1 left):** $\varepsilon_{11} = \varepsilon_{22}$ rises linearly from 0 to 0.05 over $t \in [0, 100]$ s.
- **Load case 2 — cyclic (Fig. 1 right):** triangular wave between +0.05 and −0.05 over $t \in [0, 400]$ s (ramp up to +0.05, down through zero to −0.05, back up — read the exact period off Fig. 1; it's a symmetric zig-zag).
- Time step: start with the tutorial's Δt; then check step-size sensitivity (Phase 5). Since Backward Euler is first-order, results should converge as Δt → 0.
- Note: rate-independent plasticity ⇒ "time" is a pseudo-time; only the load *path* matters, not the speed. Good interpretation point.

---

## Phase 4 — Implementation order (incremental, always testable)

Build in this sequence; each stage has an unambiguous pass criterion before moving on.

1. **Elasticity only.** Hooke's law + driver + numerical tangent, uniaxial and equi-biaxial. Pass: linear response; for equi-biaxial, analytically $\sigma_{11} = \sigma_{22} = \frac{E}{1-\nu}\,\bar\varepsilon$ and $\varepsilon_{33} = -\frac{2\nu}{1-\nu}\,\bar\varepsilon$. Check the numerical tangent equals $\mathbb{E}^e$ (assemble analytically once just for this test).
2. **Yield function module.** Compute F, G, H, L, M, N from the table; assemble $\mathbf{G}$; unit-test: (a) $\Phi(\text{uniaxial } \sigma_{11} = 580) = 0$ and likewise for all six yield stresses, (b) $\mathbf{G} \cdot [1,1,1,0,0,0]^T = \mathbf{0}$, (c) convexity conditions from §1.3 hold numerically.
3. **Von Mises perfect plasticity first** (radial return, closed-form $\mu$, lecture §5.4.2). Small, known-good stepping stone: verify uniaxial response yields at $\sigma_y$ and stays exactly on the yield surface.
4. **Hill return mapping** with local Newton. Pass criteria: $|\Phi_{n+1}| <$ tol at every plastic step; local Newton converges in ≲ 5–8 iterations (quadratic tail in the residual log — print it once for the presentation); $\text{tr}\,\boldsymbol{\varepsilon}^p = 0$ to machine precision.
5. **Von Mises recovery through the Hill machinery** (Task 3, second part): feed $F = G = H = 1/(2\sigma_y^2)$, $L = M = N = 3/(2\sigma_y^2)$ with a single $\sigma_y$ (e.g. 580) into your **Hill** implementation and compare stress/strain histories against your step-3 von Mises implementation. Pass: curves overlay to solver tolerance. This is the assignment's mandated verification.
6. **Both load cases** with the real anisotropic parameters; full logging.

---

## Phase 5 — Verification (feeds directly into the "verification" slide)

Checks, roughly in order of persuasive power:

1. **Von Mises recovery** (Phase 4.5) — the assignment's own criterion.
2. **Elastic limit:** for the equi-biaxial path, the first-yield point can be computed analytically. Pre-yield, $\sigma_{33} = 0$ and $\sigma_{11} = \sigma_{22} = \frac{E}{1-\nu}\bar\varepsilon =: s$; insert into Hill: $\Phi = 0$ at $s^* = 1/\sqrt{F\cdot 0 + G\,s^2/s^2\ldots}$ — concretely $\Phi = [G + H]\,s^2 - 1 = 0$ (the $F$-term vanishes since $\sigma_{11} = \sigma_{22}$), so $s^* = 1/\sqrt{G+H}$. Compare against the simulated onset of plasticity (kink in the σ–ε curve, first step with $\mu > 0$). Do the same check for von Mises parameters: $s^* = \sigma_y$.
3. **Consistency:** plot $\Phi(t)$ — must be $\le 0$ always, $= 0$ (to tolerance) whenever $\mu > 0$.
4. **Perfect plasticity signature:** during plastic flow the stress *point* may still move along the fixed yield surface (no hardening ⇒ the surface never grows, but the stress is not frozen — for this anisotropic biaxial case $\sigma_{11}$ and $\sigma_{22}$ redistribute). Show the stress path in the $\sigma_{11}$–$\sigma_{22}$ plane together with the yield ellipse (with $\sigma_{33} = 0$, shear = 0, the yield locus is an ellipse — plot it as an implicit contour of $\Phi = 0$). The path must ride exactly on the ellipse during plastic loading.
5. **Plastic incompressibility:** $\varepsilon^p_{11} + \varepsilon^p_{22} + \varepsilon^p_{33} = 0$ over the whole history.
6. **Dissipation:** $\boldsymbol{\sigma} : \Delta\boldsymbol{\varepsilon}^p \ge 0$ every step.
7. **Time-step convergence:** run the ramp with Δt, Δt/2, Δt/4; curves must converge (first-order BE).
8. **Tangent sanity:** the driver's global Newton should converge in few iterations; if it stalls, suspect the perturbation size or history leakage into the tangent evaluations.

---

## Phase 6 — Results & interpretation (Task 3)

Produce at least these plots (readable axes/fonts — an explicit grading criterion):

1. **Input:** $\varepsilon_{11}, \varepsilon_{22}$ vs. $t$ for both load cases (reproduce Fig. 1).
2. **Stress response:** $\sigma_{11}(t)$ and $\sigma_{22}(t)$ — the headline anisotropy result: identical strain input, **different stresses** ($\sigma_y^{11} = 580 \ne 510 = \sigma_y^{22}$). Under von Mises parameters the two curves must coincide.
3. **Plastic strains:** $\varepsilon^p_{11}(t)$, $\varepsilon^p_{22}(t)$ (and $\varepsilon^p_{33}(t)$ to show incompressibility). Interpret: unequal plastic flow despite equal total strain — the flow direction $2\mathbb{G}:\boldsymbol{\sigma}_{dev}$ is skewed by the anisotropy.
4. **Stress–strain:** $\sigma_{11}$ vs. $\varepsilon_{11}$ and $\sigma_{22}$ vs. $\varepsilon_{22}$, ramp and cyclic. Cyclic: elastic unloading with slope $\frac{E}{1-\nu}$ (biaxial modulus), re-yield in compression at the same surface (no Bauschinger effect — no kinematic hardening), closed hysteresis loops from the second cycle on.
5. **Yield locus plot:** Hill ellipse and von Mises ellipse in the $\sigma_{11}$–$\sigma_{22}$ plane (at $\sigma_{33}=0$) + simulated stress paths. Most instructive single figure for the presentation.
6. **Also show** $\varepsilon_{33}(t)$: elastic lateral contraction plus plastic contribution — discuss.

Interpretation points to prepare for the oral exam:

- Why is the return not radial for Hill? (flow direction vs. trial deviator)
- Why can't $\Delta\lambda$ be solved in closed form here, but can for von Mises? (For von Mises, $\boldsymbol{\nu} \parallel \boldsymbol{\sigma}_{dev} \parallel \boldsymbol{\sigma}^{tr}_{dev}$ ⇒ scalar equation linear in $\mu$ after normalization; for Hill the matrix $[\mathbf{I} + 2\mu\mathbf{D}\mathbf{G}]^{-1}$ makes $\Phi(\mu)$ a genuine nonlinear (rational) function.)
- Where does $\Phi = \sigma:\mathbb{G}:\sigma - 1$ sit in the lecture's taxonomy? (Quadratic-form anisotropy §6.5.1 → orthotropy via structural tensors $\boldsymbol{A}_i = \boldsymbol{a}_i \otimes \boldsymbol{a}_i$, 9 parameters → mean-stress independence, −3 constraints → 6 Hill parameters, eq. 6.57–6.68.)
- Units/normalization: Hill's $\Phi$ is dimensionless and quadratic in stress (yield surface at $\Phi = 0$ ⇔ value 1 of the quadratic form), whereas the lecture's von Mises $\Phi = \sigma_e - \sigma_y$ has stress units — hence the different $\lambda$ scaling noted in §1.3.
- Perfect plasticity limits: no hardening ⇒ under continued proportional straining stresses saturate on the yield surface; comparison hook to lecture §5.5 hardening if asked.

---

## Phase 7 — Deliverables checklist

- [ ] Derivation notes (Task 1): framework, parameter identification, return mapping, Newton scheme — clean enough to present.
- [ ] Code: well-structured (subroutines), commented, meaningful names — explicit grading criteria. Separate: material parameters / yield-function module / stress update / numerical tangent / driver / plotting.
- [ ] All plots from Phase 6, high resolution, readable fonts, labeled axes with units.
- [ ] Verification evidence from Phase 5 (at minimum: von Mises recovery overlay + $\Phi(t)$ + yield-locus ride).
- [ ] Presentation (~10 min): Title → Motivation (why anisotropy: rolled sheet metal, texture) → Model & relation to lecture (only core equations!) → Algorithm (flowchart, *not* code) → Results + verification → Conclusion. LaTeX/beamer recommended; cite Hill (1950), de Souza Neto et al. (2008), lecture notes.
- [ ] Code submitted by **Aug 28, 9:00** · slides by **Sep 2, 23:59** → tim.furlan@tu-dortmund.de.

---

## Known pitfalls (read again before debugging)

1. **Component ordering** 12/23/13 vs. Voigt 23/13/12 — a silent, plausible-looking wrongness.
2. **Shear factor 2**: never double shear stresses in $\boldsymbol{\sigma}^m$; the 2s live inside $\mathbf{G}$. But when mapping strain-like 6-vectors or building $\mathbf{D}$, engineering-vs-tensor shear conventions matter — derive once, unit-test with the pure-shear yield check (Phase 4.2a).
3. **History leakage in the numerical tangent**: perturbed stress calls must not update $\boldsymbol{\varepsilon}^p$.
4. **Commit history only after global convergence**, not inside driver iterations.
5. **Trial state from converged history** $\boldsymbol{\varepsilon}^p_n$, never from a mid-iteration candidate.
6. **Perturbation size** of the forward difference: too small kills the global Newton via noise.
7. **Local Newton start value** $\mu^{(0)} = 0$; guard against $\mu < 0$ (shouldn't happen if $\Phi^{tr} > 0$ and $f$ monotone — if it does, the bug is upstream).
8. Deviator projection: $\boldsymbol{\sigma}_{dev} = \boldsymbol{\sigma} - \tfrac{1}{3}\text{tr}(\boldsymbol{\sigma})\boldsymbol{I}$ — recompute after every stress update inside the local Newton.

## References

- [1] A. Menzel, *Theory of Materials*, Lecture Notes (kapAll-ho): §5.2–5.4 (perfect plasticity, von Mises prototype, radial return), §6.3 (BE integrator, Newton), §6.5 (anisotropic yield surfaces, orthotropy, Hill), §6.6 (prototype integrator eqs. 6.103–6.106, plate-with-hole example with **identical material parameters**).
- [2] R. Hill, *The Mathematical Theory of Plasticity*, Oxford University Press, 1950.
- [3] E.A. de Souza Neto, D. Perić, D.R.J. Owen, *Computational Methods for Plasticity*, Wiley, 2008 — ch. 6/7 for return mapping, ch. on anisotropic plasticity for Hill.
