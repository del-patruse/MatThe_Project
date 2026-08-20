#Teo Tassqa 16.07.26

import numpy as np

#this would be a struct in c# but dirty python doesn't like that
class PlasticParams:
    def __init__(self, F, G, H, L, M, N ): #
        self.F = F
        self.G = G
        self.H = H
        self.L = L
        self.M = M
        self.N = N

class ElasticParams:
    def __init__(self, E, nu):
        self.E = E
        self.nu = nu


class MatParams:
    def __init__(self, elastic_params: ElasticParams, plastic_params: PlasticParams):
        self.elastic_params = elastic_params
        self.plastic_params = plastic_params


