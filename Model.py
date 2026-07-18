# Teo Tassa 17.07.2026
import numpy as np
from mat_params import MatParams, ElasticParams, PlasticParams


'''Isotropic Hooke elasticity (E, nu)'''
class Elastic:
    def __init__(self, params: ElasticParams):
        self.params = params

    def stiffness_tensor(self):
        # E^e as 6x6 matrix in the project's component order (11,22,33,12,23,13)
        pass  # not implemented yet

    def compliance_tensor(self):
        # [E^e]^-1
        pass  # not implemented yet

    def stress(self, eps_e):
        # sigma = E^e : eps_e
        pass  # not implemented yet


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

        #fourth order structural tensor in voigt notation:
        G_tens = np.array( [[F+G, -F,   -G,   0,    0,    0],
                            [-F,  F+H,  -H,   0,    0,    0],
                            [-G,  -H,   G+H,  0,    0,    0],
                            [0,   0,    0,    2*L,  0,    0],
                            [0,   0,    0,    0,    2*M,  0],
                            [0,   0,    0,    0,    0,    2*N]])#placeholder for now
        return G_tens

    def yield_function(self, sigma_dev):
        # Phi = sigma_dev : G : sigma_dev - 1
        return 0 #placeholder

    def flow_direction(self, sigma_dev):
        # nu = dPhi/dsigma = 2 * G : sigma_dev
        pass  # not implemented yet


'''Material model: composes elasticity + Hill plasticity via return mapping'''
class Model:
    def __init__(self, params: MatParams):
        self.params = params
        self.elastic = Elastic(params.elastic_params)
        self.plastic = Plastic(params.plastic_params)

    def trial_stress(self, eps, eps_p_n):
        # sigma_tr = E^e : (eps - eps_p_n), frozen plastic flow
        pass  # not implemented yet

    def stress_update(self, eps, eps_p_n):
        # Material routine (side-effect-free): elastic predictor, local Newton
        # return mapping if plastic, returns (sigma_n1, eps_p_n1)
        pass  # not implemented yet

    def numerical_tangent(self, eps, eps_p_n):
        # Forward-difference tangent via repeated stress_update calls,
        # all from the same eps_p_n
        pass  # not implemented yet


'''
Class Strain containing the elastic and plastic strain components
'''
class Strain:
    def __init__(self, eps_e, eps_p):
        self.eps_e = eps_e
        self.eps_p = eps_p
        self.eps = eps_e + eps_p
