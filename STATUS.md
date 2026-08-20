# Project Status

Snapshot of implementation progress against `ROADMAP.md` (the assignment plan) and the assignment PDF (`Pruefung02.pdf`). Updated 2026-08-17.

## Update log

- 2026-08-18 — `verify_von_mises_recovery` implemented in `plots_and_verifs.py` (was a diff-and-print sketch). Now runs four checks and returns a single pass/fail bool: (1) the two `G` tensors coincide — the algebraic content of eqs. (11)-(12); (2) stress / plastic-strain / total-strain histories overlay, reported as absolute *and* relative max differences; (3) the specialized run is isotropic (`sigma11 == sigma22` under equi-biaxial load); (4) every plastic step sits at `sigma_eq = sqrt(3/2 s:s) = sigma_y`, computed from the definition rather than through `G`, so it is not just the same code path again. Signature gained optional `plastic_hill` / `plastic_vm` / `sigma_y` / `rtol` / `show_plot`; `main.py`'s `task=1` passes all of them. Plot helper `_plot_von_mises_recovery` draws the mandated 3-panel overlay (stress, plastic strain, semilog pointwise difference). Verified end-to-end: ramp and cyclic both pass (checks 1-2 at exactly 0.0, 3-4 at ~1e-13, i.e. machine precision), and a negative control with the true anisotropic parameters fails all four as it should. Caveat recorded below.
- 2026-08-17 — `main.py` task switch collapsed from three cases to two: `task=1` now runs the Hill case plus the von Mises recovery check (formerly `task=2`'s job; the old standalone "just run Hill" `task=1` is gone), `task=2` runs the plots (formerly `task=3`). `time_points` is no longer a placeholder — `deltaT = 1.0` and `np.arange(...)` (0-100 for ramp, 0-400 for cyclic) are wired in for real. The `build_von_mises_params` keyword bug (below) is fixed. `plot_lateral_contraction` is no longer imported/called from `main.py` (still defined in `plots_and_verifs.py`, just unwired). Added `SEMESTER6_REFERENCE_CHECK.md`: audit of `Semester 6/Material theory/` confirming no plasticity/return-mapping/Hill/von-Mises material exists there to port — `Model.stress_update` has to be derived from the lecture notes, not adapted from tutorial code.
- 2026-08-12 — `Solver.Driver` implemented (`partition`, `residual`, `global_newton_step`, `step`, `run`), ported from the tutorial's `uniaxialStressDriver.uniaxialStress` (`Semester 6/Material theory/Ubung4_muster`) and generalized to the equi-biaxial prescribed set `{0, 1}` (= 11, 22). Verified structurally (imports, `partition()` splits correctly into free `[2,3,4,5]` / prescribed `[0,1]`) but not runnable end-to-end yet — it calls into `Model.stress_update`/`Model.numerical_tangent`, which are still stubs.
- 2026-08-12 — `verify_von_mises_recovery` and the six `plot_*` stubs moved out of `main.py` into a new `plots_and_verifs.py`; `main.py` now only imports them. `main()` restructured around a module-level `task` switch (`match task: case 1/2/3`) selecting between (1) running both load cases through Hill, (2) running the von Mises-recovery check, (3) running the ramp case and calling the plot functions.
- 2026-08-12 — `plots_and_verifs.py`: all six `plot_*` functions implemented (each unpacks `History.to_arrays()` and plots the relevant quantity vs. `t`, or stress vs. strain / the yield locus contour for `plot_stress_strain`/`plot_yield_locus`), each now has an axis title. `verify_von_mises_recovery` was briefly filled in (diff sigma/eps_p between Hill and von Mises histories, prints max diffs, overlays the sigma11 curves) but is back to a `pass`-style placeholder on disk — not yet reimplemented.
- 2026-08-12 — `main.py:31` bug fixed: `sy23` corrected from `600` to `530` (Table 1 of `Pruefung02.pdf`). The `main.py:43` `build_von_mises_params(sy11=580)` keyword-mismatch bug is still present (see Confirmed bugs below).
- 2026-08-12 — `utils.py` gained a `NewtonRaphsonSolver` class (generic damped Newton solver: `solve(func, jacobian, x0)`). Not yet wired into `Model.stress_update`/`trial_stress`, which are still stubs — presumably intended for the local return-mapping iteration.

## Architecture

Code is split into: `mat_params.py` (parameter containers), `utils.py` (tensor/Voigt-6 math), `Model.py` (material point: `Elastic`, `Plastic`, `Model`, `Strain`), `Solver.py` (parameter builders, `StepRecord`, `History`, `Driver` — the constitutive driver/global solver), `plots_and_verifs.py` (verification + plotting functions), `main.py` (material parameter values, load curves, task dispatch, entry point). `SEMESTER6_REFERENCE_CHECK.md` is a one-off audit doc (not code) checking whether the Semester 6 tutorial folder has anything plasticity-relevant left to port — conclusion: no.

## Done

**`mat_params.py`** — `ElasticParams`, `PlasticParams`, `MatParams`: plain parameter containers. Complete.

**`utils.py`** — tensor/Voigt-6 plumbing, fully implemented:
`tensor_to_voigt6`, `voigt6_to_tensor`, `double_contraction_2nd`, `double_contraction_4th_2nd`, `double_contraction_2nd_4th`, `double_contraction_voigt6`, `double_contraction_voigt6_4th_2nd`, `quadratic_form_voigt6`, `identity_2nd`, `identity_4th_sym`, `identity_4th_sym_dev`, `outer_2nd`, `trace`, `deviator`, `symmetrize`, `frobenius_norm`. Note: component order is (11, 22, 33, 12, 23, 13), not standard Voigt order — documented at the top of the file. Also has a standalone `NewtonRaphsonSolver` class (damped Newton, generic `func`/`jacobian`/`x0`) — implemented but not called from anywhere yet.

**`Model.py` — `Elastic`**: `stiffness_tensor`, `compliance_tensor`, `stress` all implemented (isotropic Hooke's law).

**`Model.py` — `Plastic`**: `calculate_G_tensor`, `yield_function` (`Phi = sigma_dev : G : sigma_dev - 1`), `flow_direction` (`nu = 2 * G : sigma_dev`) all implemented, matching this assignment's specific F/G/H, L/M/N convention from the PDF (eqs. 3-6) — not the generic Hill-1948 textbook layout; the two differ and it matters.

**`Model.py` — `Strain`**: value object (`eps_e`, `eps_p`, `eps`) plus `from_total_and_plastic` / `from_total_and_elastic` constructors. Complete.

**`Solver.py`** — `hill_params_from_yield_stresses`, `elastic_params_from_E_nu`, `build_von_mises_params` (F=G=H, L=M=N special case). Complete. `StepRecord` data container defined. `History.update_array` appends a record; `History.to_arrays` unpacks all records into per-quantity numpy arrays (stacks `eps`/`eps_p`/`sigma` as `(n_steps, 3, 3)`, scalars as `(n_steps,)`) — implemented and matches `StepRecord`'s field order.

**`Solver.py` — `Driver`**: strain/stress-mixed constitutive driver, implemented.
- `partition()` — splits the 6 Voigt components (`utils.VOIGT_INDEX_PAIRS` order) into `prescribed` (from `self.prescribed_components`, e.g. `{0, 1}` for equi-biaxial 11/22) and `free` (stress-driven, driven to zero).
- `residual()` — assembles the full strain tensor from `eps_bar` (prescribed) + a free-strain guess, calls `model.stress_update`, returns the stress residual restricted to the free components.
- `global_newton_step()` — gets the full 6x6 tangent from `model.numerical_tangent`, restricts it to the free x free block, solves one Newton update (mirrors the tutorial driver's `CKelBar`/`sigKelBar` partitioning).
- `step()` — iterates `global_newton_step` to convergence (residual checked *before* each update, same order as the tutorial code), returns a `StepRecord`.
- `run()` — marches through `time_points`, evaluates the load curve, calls `step`; only commits `eps_p` / warm-starts the next free-strain guess *after* convergence (Phase 2 rule: commit history only after global convergence).
- Fixes the interface contract `Model.stress_update(eps, eps_p_n) -> (sigma, eps_p_new, info)` and `Model.numerical_tangent(eps, eps_p_n) -> 6x6 ndarray`, with `info` carrying `phi`, `mu`, `n_local_iter`, `is_plastic` — documented in the `Driver` docstring for whoever implements `Model.py` next.

**`main.py` — `ramp_load`**: `0 -> 0.05` linear ramp over `t in [0, 100]`, matches Fig. 1 (left). Complete.

**`main.py` — `cyclic_load`**: triangular wave via closed form `ampl * (2/pi) * arcsin(sin(2*pi*f*t))`, `f = 0.01` Hz (period 100s, read off Fig. 1 right — 4 full cycles over 400s), amplitude 0.05, starts at 0 and rises first. Complete.

**`main.py` — `main()` task dispatch**: `task` (module-level int, currently `1`) selects behaviour via `match task: case 1/2` — `1` runs Hill + von Mises recovery check, `2` runs the plots. `time_points` is now a real array (`deltaT = 1.0`, `np.arange`: 0-100 for ramp, 0-400 for cyclic) — structurally complete and no longer blocked on time discretization. Still not runnable end-to-end: every path hits the `Model.stress_update` stub regardless.

**`plots_and_verifs.py` — plot functions**: `plot_load_history`, `plot_stress_response`, `plot_plastic_strains`, `plot_stress_strain`, `plot_yield_locus`, `plot_lateral_contraction` all implemented — each unpacks `History.to_arrays()`, plots the relevant quantity/quantities (vs. `t`, vs. strain, or as the `plastic.G_tensor` yield-locus contour with the stress path overlaid), and has axis labels, a legend, and a title. `plot_lateral_contraction` is implemented but currently unwired — `main.py`'s `task=2` case doesn't import or call it.

**`plots_and_verifs.py` — `verify_von_mises_recovery`**: implemented, with the private helpers `_difference` (absolute + relative max difference) and `_von_mises_equivalent_stress` (`sqrt(3/2 s:s)` over a whole history, computed from the definition, not via `G`) and the 3-panel overlay plot `_plot_von_mises_recovery`. Passes on ramp and cyclic; a negative control with the true anisotropic parameters fails it. See the To Do caveat about what the check does and does not prove.

## To Do

Roughly Phase 4 of `ROADMAP.md` — the return-mapping algorithm and the von Mises-recovery check are still stubs; plots are done:

- **`Model.trial_stress`, `Model.stress_update`, `Model.numerical_tangent`** — not implemented (`pass` stubs). This is the core of the assignment: elastic predictor, local Newton return mapping (von Mises closed-form and/or Hill iterative, likely via `utils.NewtonRaphsonSolver`), forward-difference tangent. `Solver.Driver` already assumes and calls into the interface described above — nothing can be simulated end-to-end until these exist.
- ~~**`verify_von_mises_recovery`**~~ — implemented and passing on both load cases (see the 2026-08-18 update-log entry). This was the *only* verification explicitly required by `Pruefung02.md` Task 3 ("check if under usage of equations (11) and (12), the von Mises flow criterion can be recovered").
  - **Caveat worth knowing before the oral exam:** checks 1 and 2 come out at *exactly* 0.0, not merely within tolerance. That is not suspicious, but it is weaker evidence than it looks. `hill_params_from_yield_stresses(sy, sy, sy, sy/sqrt(3), sy/sqrt(3), sy/sqrt(3))` and `build_von_mises_params(sy)` produce bit-identical floats, so both `Driver` runs execute the same code on the same `G` — the history comparison is a determinism check, not an independent one. ROADMAP Phase 4.3 assumed a *separate* closed-form radial-return von Mises implementation as the reference; that was never written. Checks 3 and 4 exist precisely to cover the gap (they test isotropy and the surface itself, independently of `G`), but building the standalone radial-return reference would make the ROADMAP's original argument properly. Cheap to add if time allows.
- **Plot functions** — implemented in `plots_and_verifs.py`, each unpacking `History.to_arrays()` and plotting with axis labels, legend, and title. `plot_lateral_contraction` is written but not yet called from `main.py`'s `task=2` case — worth wiring back in since Task 3's stress-strain/lateral-contraction plots are part of the deliverable. Not runnable end-to-end yet — blocked on `Model.stress_update`, same as everything else.
- **`time_points`** — done. `main.py` now builds a real array (`deltaT = 1.0`, `np.arange(0, 100+deltaT, deltaT)` for ramp / `np.arange(0, 400+deltaT, deltaT)` for cyclic). `ROADMAP.md` Phase 3's step-size-sensitivity check is still open but not blocking.
- **No unit tests** exist yet — `ROADMAP.md` §Phase 4.2 calls for checking `Phi` at each of the six uniaxial/shear yield stresses and `G · [1,1,1,0,0,0]^T = 0` before trusting the yield module further. Not required by `Pruefung02.md`, but cheap and worth doing once the yield module is exercised.

## Scope note (checked against `Pruefung02.md`, the actual exam task, 2026-08-12)

`ROADMAP.md` is a self-authored study plan and over-scoped the verification requirements. The exam document only mandates: (1) illustrative plots of strains/stresses over time, (2) von Mises recovery under eqs. (11)-(12). Removed `verify_elastic_limit`, `verify_consistency`, `verify_plastic_incompressibility`, `verify_dissipation`, `verify_timestep_convergence` from `main.py` as unrequired stubs — they were good practice but not deliverables, and `verify_consistency`/`verify_plastic_incompressibility` in particular can be revisited cheaply later for the oral exam if time allows.

## Confirmed bugs to fix

- ~~`main.py:31` sets `sy23 = 600`~~ — fixed, now `530` (matches Table 1 of `Pruefung02.pdf`, p.5).
- ~~`main.py:43` calls `build_von_mises_params(sy11=580)`~~ — fixed, now `build_von_mises_params(sigma_y=sy11)` (`main.py:67`), matching the `Solver.py:26` signature (`build_von_mises_params(sigma_y)`) and no longer hardcoding `580` separately from `sy11`.
- None currently confirmed open. Everything below "None" is speculative until `Model.stress_update` exists and the code can actually be run end-to-end — worth re-auditing once it does.
