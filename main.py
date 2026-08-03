# Teo Tassa 18.07.2026
import numpy as np
import matplotlib.pyplot as plt

from mat_params import MatParams, ElasticParams, PlasticParams
from Model import Model, Plastic, Strain



#task specific param vaLUES: 
E = 200000
nu = .3
sy11 = 580
sy22 = 510
sy33 = 460
sy12 = 590
sy23 = 600
sy13 = 550




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


#ToDo
def build_von_mises_params(sigma_y) -> PlasticParams:
    F = G = H = 1/(2 * sigma_y**2)
    L = M = N = 3/(2 * sigma_y**2)
    return PlasticParams(F=F, G=G, H=H, L=L, M=M, N=N)

hillMatParams =  MatParams(elastic_params=elastic_params_from_E_nu(E, nu), plastic_params=hill_params_from_yield_stresses(sy11, sy22, sy33, sy12, sy23, sy13))





#defining löoads
def ramp_load(t)-> float:
    # 0 -> 0.05 over t in [0, 100]
    m = 0.05 / 100
    return m * t
    pass 


def cyclic_load(t)-> float:
    
    pass  # not implemented yet


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

    def append(self, record: StepRecord):
        pass  # not implemented yet

    def to_arrays(self):
        # unpack self.records into per-quantity time-series arrays for plotting
        pass  # not implemented yet


# ---------------------------------------------------------------------------
# Phase 0.2 / 2 / 3 — strain/stress-mixed constitutive driver
# ---------------------------------------------------------------------------
class Driver:
    def __init__(self, model: Model, prescribed_components):
        self.model = model
        self.prescribed_components = prescribed_components  # equi-biaxial: {11, 22}

    def partition(self):
        # split the 6 components into strain-driven vs stress-driven (stress-free)
        pass  # not implemented yet

    def residual(self, eps_free_guess, eps_bar, eps_p_n):
        # stress components that must vanish, given the current free-strain guess
        pass  # not implemented yet

    def global_newton_step(self, eps_free_guess, eps_bar, eps_p_n):
        # one global Newton update using Model.numerical_tangent
        pass  # not implemented yet

    def step(self, eps_bar, eps_p_n):
        # iterate global_newton_step to convergence; return converged eps, sigma, eps_p
        pass  # not implemented yet

    def run(self, load_curve, time_points) -> History:
        # march through time_points; commit eps_p to history only after global convergence
        pass  # not implemented yet


# ---------------------------------------------------------------------------
# Phase 5 — verification
# ---------------------------------------------------------------------------
def verify_von_mises_recovery(history_hill: History, history_vm: History):
    pass  # not implemented yet


def verify_elastic_limit(history: History, params: MatParams):
    pass  # not implemented yet


def verify_consistency(history: History):
    # Phi <= 0 always, Phi ~= 0 whenever mu > 0
    pass  # not implemented yet


def verify_plastic_incompressibility(history: History):
    pass  # not implemented yet


def verify_dissipation(history: History):
    pass  # not implemented yet


def verify_timestep_convergence(driver: Driver, load_curve, dt_values):
    pass  # not implemented yet


# ---------------------------------------------------------------------------
# Phase 6 — plots
# ---------------------------------------------------------------------------
def plot_load_history(history: History):
    pass  # not implemented yet


def plot_stress_response(history: History):
    pass  # not implemented yet


def plot_plastic_strains(history: History):
    pass  # not implemented yet


def plot_stress_strain(history: History):
    pass  # not implemented yet


def plot_yield_locus(history: History, plastic: Plastic):
    pass  # not implemented yet


def plot_lateral_strain(history: History):
    pass  # not implemented yet


def main():
    params = hillMatParams
    model = Model(params)
    driver = Driver(model, prescribed_components={0, 1})

    # ramp_history = driver.run(ramp_load, time_points=...)
    # cyclic_history = driver.run(cyclic_load, time_points=...)

    # vm_params = build_von_mises_params(sigma_y=580)
    # vm_history = Driver(Model(vm_params), prescribed_components={0, 1}).run(ramp_load, time_points=...)
    # verify_von_mises_recovery(ramp_history, vm_history)

    # plot_load_history(ramp_history)
    # plot_stress_response(ramp_history)
    # ...

    plt.show()


if __name__ == "__main__":
    main()
