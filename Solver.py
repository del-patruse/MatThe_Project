# Teo Tassa 03.08.2026
import numpy as np

from mat_params import ElasticParams, PlasticParams
from Model import Model
from utils import tensor_to_voigt6, voigt6_to_tensor


#-------------------------------------------
#defining material parameters
#-------------------------------------------
def hill_params_from_yield_stresses(sy11, sy22, sy33, sy12, sy23, sy13) -> PlasticParams:
    F = .5 * (1/sy11**2 + 1/sy22**2 - 1/sy33**2)
    G = .5 * (1/sy11**2 + 1/sy33**2 - 1/sy22**2)
    H = .5 * (1/sy22**2 + 1/sy33**2 - 1/sy11**2)
    L = 1/(2 * sy12**2)
    M = 1/(2 * sy23**2)
    N = 1/(2 * sy13**2)
    return PlasticParams(F=F, G=G, H=H, L=L, M=M, N=N)

def elastic_params_from_E_nu(E, nu) -> ElasticParams:
    return ElasticParams(E=E, nu=nu)


def build_von_mises_params(sigma_y) -> PlasticParams:
    F = G = H = 1/(2 * sigma_y**2)
    L = M = N = 3/(2 * sigma_y**2)
    return PlasticParams(F=F, G=G, H=H, L=L, M=M, N=N)



# ---------------------------------------------------------------------------
# Phase 2 — per-step log record
# ---------------------------------------------------------------------------
class StepRecord:
    def __init__(self, t, eps, eps_p, sigma, phi, mu, n_local_iter, is_plastic):
        self.t = t
        self.eps = eps
        self.eps_p = eps_p
        self.sigma = sigma
        self.phi = phi
        self.mu = mu
        self.n_local_iter = n_local_iter
        self.is_plastic = is_plastic


class History:
    def __init__(self):
        self.records = []

    def update_array(self, record: StepRecord):
        return self.records.append(record)

    def to_arrays(self):
        t_array = np.array([record.t for record in self.records])
        eps_array = np.array([record.eps for record in self.records])
        eps_p_array = np.array([record.eps_p for record in self.records])
        sigma_array = np.array([record.sigma for record in self.records])
        phi_array = np.array([record.phi for record in self.records])
        mu_array = np.array([record.mu for record in self.records])
        n_local_iter_array = np.array([record.n_local_iter for record in self.records])
        is_plastic_array = np.array([record.is_plastic for record in self.records])

        return t_array, eps_array, eps_p_array, sigma_array, phi_array, mu_array, n_local_iter_array, is_plastic_array


# ---------------------------------------------------------------------------
# Phase 0.2 / 2 / 3 — strain/stress-mixed constitutive driver
# ---------------------------------------------------------------------------
def _extract(value, indices):
    # Partition a 6-vector or 6x6 matrix (this project's Voigt-6 order,
    # see utils.VOIGT_INDEX_PAIRS) to the given component indices.
    if value.shape == (6,):
        return value[indices]
    elif value.shape == (6, 6):
        return value[np.ix_(indices, indices)]
    else:
        raise ValueError("Tensor shape error.")


class Driver:
    """
    Strain/stress-mixed constitutive driver, following the pattern from the exercise 

    Expected M

    """

    def __init__(self, model: Model, prescribed_components, tol=1e-8, max_iter=30):
        self.model = model
        self.prescribed_components = prescribed_components  # equi-biaxial: {0, 1} (=11, 22)
        self.tol = tol
        self.max_iter = max_iter

    def partition(self):
        # split the 6 components into strain-driven (prescribed) vs stress-driven (free)
        prescribed = np.array(sorted(self.prescribed_components))
        free = np.array([i for i in range(6) if i not in self.prescribed_components])
        return free, prescribed

    def residual(self, eps_free_guess, eps_bar, eps_p_n):
        # stress components that must vanish, given the current free-strain guess
        free, prescribed = self.partition()

        eps_voigt = np.zeros(6)
        eps_voigt[prescribed] = eps_bar
        eps_voigt[free] = eps_free_guess
        eps = voigt6_to_tensor(eps_voigt)

        sigma, eps_p_candidate, info = self.model.stress_update(eps, eps_p_n)
        R = _extract(tensor_to_voigt6(sigma), free)

        return R, eps, sigma, eps_p_candidate, info

    def global_newton_step(self, eps_free_guess, eps_bar, eps_p_n):
        # one global Newton update using Model.numerical_tangent
        free, _ = self.partition()

        R, eps, sigma, eps_p_candidate, info = self.residual(eps_free_guess, eps_bar, eps_p_n)
        C = self.model.numerical_tangent(eps, eps_p_n)
        C_free = _extract(C, free)

        eps_free_new = eps_free_guess - np.linalg.solve(C_free, R)

        return eps_free_new, R, eps, sigma, eps_p_candidate, info

    def step(self, t, eps_bar, eps_p_n, eps_free_guess=None, verbose=False) -> StepRecord:
        # iterate global_newton_step to convergence; return converged eps, sigma, eps_p
        free, _ = self.partition()
        if eps_free_guess is None:
            eps_free_guess = np.zeros(len(free))

        for i in range(self.max_iter):
            eps_free_new, R, eps, sigma, eps_p_candidate, info = self.global_newton_step(
                eps_free_guess, eps_bar, eps_p_n)

            if verbose:
                print('  t =', t, ' iter ', i, ', residual = ', "%.3e" % np.linalg.norm(R))

            if np.linalg.norm(R) <= self.tol:
                return StepRecord(t=t, eps=eps, eps_p=eps_p_candidate, sigma=sigma,
                                   phi=info['phi'], mu=info['mu'],
                                   n_local_iter=info['n_local_iter'],
                                   is_plastic=info['is_plastic'])

            eps_free_guess = eps_free_new

            if i == self.max_iter - 1:
                raise RuntimeError(f"Global Newton scheme did not converge at t={t}.")

    def run(self, load_curve, time_points) -> History:
        # march through time_points; commit eps_p to history only after global convergence
        history = History()
        eps_p_n = np.zeros((3, 3))
        eps_free_guess = None
        free, _ = self.partition()

        for t in time_points:
            eps_bar = load_curve(t)
            record = self.step(t, eps_bar, eps_p_n, eps_free_guess=eps_free_guess)
            history.update_array(record)

            # commit history and warm-start the next step's free-strain guess
            eps_p_n = record.eps_p
            eps_free_guess = _extract(tensor_to_voigt6(record.eps), free)

        return history
