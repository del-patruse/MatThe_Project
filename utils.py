import numpy as np

# Component order used throughout this project for the 6-vector
# representation of a symmetric 3x3 tensor: (11, 22, 33, 12, 23, 13).
# NOTE: this is NOT standard Voigt order (which puts 23 before 13/12).
VOIGT_INDEX_PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))

    
def identity_2nd(): 
    return np.eye(3)


def identity_4th_sym() -> np.ndarray:
    # I^sym_ijkl = 1/2 (delta_ik delta_jl + delta_il delta_jk)
    I = np.eye(3)
    return 0.5 * (np.einsum('ik,jl->ijkl', I, I) + np.einsum('il,jk->ijkl', I, I))


def identity_4th_sym_dev() -> np.ndarray:
    # Deviatoric projector: I^sym_dev = I^sym - 1/3 (I x I)
    return identity_4th_sym() - (1.0 / 3.0) * outer_2nd(identity_2nd(), identity_2nd())


def outer_2nd(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    # Dyadic product of two 2nd-order tensors -> 4th-order tensor (A x B)_ijkl = A_ij B_kl
    return np.einsum('ij,kl->ijkl', A, B)


def double_contraction_2nd(A: np.ndarray, B: np.ndarray) -> float:
    # A : B = A_ij B_ij
    return float(np.tensordot(A, B, axes=2))


def double_contraction_4th_2nd(C: np.ndarray, A: np.ndarray) -> np.ndarray:
    # (C : A)_ij = C_ijkl A_kl
    return np.tensordot(C, A, axes=([2, 3], [0, 1]))


def double_contraction_2nd_4th(A: np.ndarray, C: np.ndarray) -> np.ndarray:
    # (A : C)_kl = A_ij C_ijkl
    return np.tensordot(A, C, axes=([0, 1], [0, 1]))


def trace(T: np.ndarray) -> float:
    return float(np.trace(T))


def deviator(T: np.ndarray) -> np.ndarray:
    return T - (trace(T) / 3.0) * identity_2nd()


def symmetrize(T: np.ndarray) -> np.ndarray:
    return 0.5 * (T + T.T)


def frobenius_norm(T: np.ndarray) -> float:
    return float(np.sqrt(double_contraction_2nd(T, T)))


def tensor_to_voigt6(T: np.ndarray) -> np.ndarray:
    # Symmetric 3x3 tensor -> 6-vector, pure tensor components (no shear factor of 2).
    return np.array([T[i, j] for i, j in VOIGT_INDEX_PAIRS])


def voigt6_to_tensor(v: np.ndarray) -> np.ndarray:
    # Inverse of tensor_to_voigt6: 6-vector of pure tensor components -> symmetric 3x3 tensor.
    T = np.zeros((3, 3))
    for (i, j), value in zip(VOIGT_INDEX_PAIRS, v):
        T[i, j] = value
        T[j, i] = value
    return T


def double_contraction_voigt6(vA: np.ndarray, vB: np.ndarray) -> float:
    # Tensor double contraction A:B from two 6-vectors of pure tensor components
    # (shear entries counted twice, since off-diagonal tensor entries appear twice in A_ij B_ij).
    normal = np.dot(vA[:3], vB[:3])
    shear = np.dot(vA[3:], vB[3:])
    return float(normal + 2.0 * shear)


def quadratic_form(v: np.ndarray, M: np.ndarray) -> float:
    # v^T M v
    return float(v @ M @ v)
