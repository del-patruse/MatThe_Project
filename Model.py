# Teo Tassa 17.07.2026
import numpy as np
from mat_params import MatParams, ElasticParams, PlasticParams
from utils import voigt6_to_tensor, double_contraction_voigt6, double_contraction_voigt6_4th_2nd, quadratic_form_voigt6, tensor_to_voigt6, identity_2nd, identity_4th_sym, identity_4th_sym_dev, outer_2nd, double_contraction_2nd, double_contraction_4th_2nd, double_contraction_2nd_4th, trace, deviator, symmetrize, frobenius_norm, NewtonRaphsonSolver


'''Isotropic Hooke elasticity (E, nu)'''
class Elastic:
    def __init__(self, params: ElasticParams):
        self.params = params

    def stiffness_tensor(self) -> np.ndarray:
   
        E = self.params.E
        nu = self.params.nu

        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
        mu = E / (2 * (1 + nu))

        E_tens = np.array([[lam + 2*mu, lam,        lam,         0,    0,     0],
                            [lam,        lam + 2*mu, lam,        0,    0,     0],
                            [lam,        lam,        lam + 2*mu, 0,    0,     0],
                            [0,          0,          0,          2*mu, 0,     0],
                            [0,          0,          0,          0,    2*mu,  0],
                            [0,          0,          0,          0,    0,     2*mu]])
        return E_tens

    def compliance_tensor(self) -> np.ndarray:
        return np.linalg.inv(self.stiffness_tensor())

    def stress(self, eps_e) -> np.ndarray:

        return voigt6_to_tensor(double_contraction_voigt6_4th_2nd(self.stiffness_tensor(), tensor_to_voigt6(eps_e)))


'''Hill's anisotropic yield criterion (F, G, H, L, M, N)'''
class Plastic:
    def __init__(self, params: PlasticParams):
        self.params = params
        self.G_tensor = self.calculate_G_tensor()

    def calculate_G_tensor(self):
        # Calculate the G tensor based on the parameters
        F = self.params.F
        G = self.params.G
        H = self.params.H
        L = self.params.L
        M = self.params.M
        N = self.params.N

        # Fourth-order Hill structural tensor in Voigt notation (component order
        # 11, 22, 33, 12, 23, 13), following THIS assignment's convention (ROADMAP.md,
        # assignment eqs. 3-6): Phi = F(s11-s22)^2 + G(s11-s33)^2 + H(s22-s33)^2
        #                             + 2L s12^2 + 2M s23^2 + 2N s13^2
        G_tens = np.array( [[F+G, -F,   -G,   0,    0,    0],
                            [-F,  F+H,  -H,   0,    0,    0],
                            [-G,  -H,   G+H,  0,    0,    0],
                            [0,   0,    0,    2*L,  0,    0],
                            [0,   0,    0,    0,    2*M,  0],
                            [0,   0,    0,    0,    0,    2*N]])
        return G_tens

    def yield_function(self, sigma_dev):
        # Phi = sigma_dev : G : sigma_dev - 1
        v = tensor_to_voigt6(sigma_dev)
        Phi = quadratic_form_voigt6(v, self.G_tensor) - 1.0
        return Phi

    def flow_direction(self, sigma_dev):
        # nu = dPhi/dsigma = 2 * G : sigma_dev
        v = tensor_to_voigt6(sigma_dev)
        nu = voigt6_to_tensor(2.0 * double_contraction_voigt6_4th_2nd(self.G_tensor, v))
        return nu


'''Material model: composes elasticity + Hill plasticity via return mapping'''
class Model:
    def __init__(self, params: MatParams):
        self.params = params
        self.elastic = Elastic(params.elastic_params)
        self.plastic = Plastic(params.plastic_params)

    def trial_stress(self, eps, eps_p_n):
        # sigma_tr = E^e : (eps - eps_p_n), frozen plastic flow
        return self.elastic.stress(eps - eps_p_n)

    def stress_update(self, eps, eps_p_n):
        # Material routine (side-effect-free): elastic predictor, local Newton
        # return mapping if plastic, returns (sigma_n1, eps_p_n1)
        sigma_trial = self.trial_stress(eps, eps_p_n)
        phi_trial = self.plastic.yield_function(deviator(sigma_trial))

        if phi_trial <= 0.0:
            info = {'phi': phi_trial, 'mu': 0.0, 'n_local_iter': 0, 'is_plastic': False}
            return sigma_trial, np.array(eps_p_n, copy=True), info

        stiffness = self.elastic.stiffness_tensor()
        G_tensor = self.plastic.G_tensor
        sigma_trial_voigt = tensor_to_voigt6(sigma_trial)
        identity = np.eye(6)
        iteration_count = 0

        def corrected_stress(mu):
            system = identity + 2.0 * mu * stiffness.dot(G_tensor)
            return np.linalg.solve(system, sigma_trial_voigt)

        def consistency(mu_array):
            nonlocal iteration_count
            iteration_count += 1
            sigma_voigt = corrected_stress(mu_array[0])
            return np.array([quadratic_form_voigt6(sigma_voigt, G_tensor) - 1.0])

        def consistency_jacobian(mu_array):
            mu = mu_array[0]
            system = identity + 2.0 * mu * stiffness.dot(G_tensor)
            sigma_voigt = np.linalg.solve(system, sigma_trial_voigt)
            derivative = -np.linalg.solve(system, 2.0 * stiffness.dot(G_tensor).dot(sigma_voigt))
            return np.array([[2.0 * sigma_voigt.dot(G_tensor).dot(derivative)]])

        solver = NewtonRaphsonSolver(tol=1e-10, max_iter=30)
        mu = float(solver.solve(consistency, consistency_jacobian, np.array([0.0]))[0])
        sigma = voigt6_to_tensor(corrected_stress(mu))
        eps_p = eps_p_n + mu * self.plastic.flow_direction(deviator(sigma))
        phi = self.plastic.yield_function(deviator(sigma))
        info = {'phi': phi, 'mu': mu, 'n_local_iter': iteration_count, 'is_plastic': True}
        return sigma, eps_p, info

    def numerical_tangent(self, eps, eps_p_n):
        # Forward-difference tangent via repeated stress_update calls,
        # all from the same eps_p_n
        delta = 1e-8
        sigma, _, _ = self.stress_update(eps, eps_p_n)
        sigma_voigt = tensor_to_voigt6(sigma)
        tangent = np.zeros((6, 6))

        for i in range(6):
            perturbation = np.zeros(6)
            perturbation[i] = delta
            sigma_perturbed, _, _ = self.stress_update(
                eps + voigt6_to_tensor(perturbation), eps_p_n)
            tangent[:, i] = (tensor_to_voigt6(sigma_perturbed) - sigma_voigt) / delta

        return tangent


'''
Class Strain containing the elastic and plastic strain components
'''
class Strain:
    def __init__(self, eps_e, eps_p):
        self.eps_e = eps_e
        self.eps_p = eps_p
        self.eps = eps_e + eps_p

    @staticmethod
    def from_total_and_plastic(eps, eps_p):
        return Strain(eps_e=eps - eps_p, eps_p=eps_p)

    @staticmethod
    def from_total_and_elastic(eps, eps_e):
        return Strain(eps_e=eps_e, eps_p=eps - eps_e)
