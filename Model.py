Teo Tassa 17.07.2026
import numpy as np

class Model:
    def __init__(self, params: MatParams):
        self.params = params
        self.G_tensor = self.calculate_G_tensor()

    def flow_criterion(self, sigma_dev):

        return 0 #placeholder
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
                            [0,   0     0,    0,    0,    2*N]])#placeholder for now
        return G_tens



class Strain
     def __init__(self, eps_e, eps_p):
        self.eps_e = eps_e
        self.eps_p = eps_p 

