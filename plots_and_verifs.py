# Teo Tassa 11.08.2026
import numpy as np
import matplotlib.pyplot as plt

from mat_params import MatParams
from Model import Plastic
from Solver import History


# ---------------------------------------------------------------------------
# plot Hill and von Mises results together
# ---------------------------------------------------------------------------
def verify_von_mises_recovery(history_hill: History, history_vm: History):
    t, _, eps_p_hill, sigma_hill, _, _, _, _ = history_hill.to_arrays()
    _, _, eps_p_vm, sigma_vm, _, _, _, _ = history_vm.to_arrays()
    plot_von_mises_recovery(t, sigma_hill, sigma_vm, eps_p_hill, eps_p_vm)


def plot_von_mises_recovery(t, sigma_hill, sigma_vm, eps_p_hill, eps_p_vm):
    # overlay figure for Hill and von Mises
    fig, (ax_sigma, ax_eps_p) = plt.subplots(2, 1)

    ax_sigma.plot(t, sigma_hill[:, 0, 0], '-', label="sigma11")
    ax_sigma.plot(t, sigma_hill[:, 1, 1], '-', label="sigma22")
    ax_sigma.plot(t, sigma_vm[:, 0, 0], '--', label="sigma11 (von Mises)")
    ax_sigma.plot(t, sigma_vm[:, 1, 1], '--', label="sigma22 (von Mises)")
    ax_sigma.set_ylabel("stress")
    ax_sigma.legend()
    ax_sigma.set_title("von Mises recovery: stress overlay")

    ax_eps_p.plot(t, eps_p_hill[:, 0, 0], '-', label="eps_p11 (Hill)")
    ax_eps_p.plot(t, eps_p_hill[:, 2, 2], '-', label="eps_p33 (Hill)")
    ax_eps_p.plot(t, eps_p_vm[:, 0, 0], '--', label="eps_p11 (von Mises)")
    ax_eps_p.plot(t, eps_p_vm[:, 2, 2], '--', label="eps_p33 (von Mises)")
    ax_eps_p.set_ylabel("plastic strain")
    ax_eps_p.legend()
    ax_eps_p.set_title("von Mises recovery: plastic strain overlay")

    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def plot_load_history(history: History):
    t, eps, _, _, _, _, _, _ = history.to_arrays()

    plt.figure()
    plt.plot(t, eps[:, 0, 0], label="eps11")
    plt.plot(t, eps[:, 1, 1], label="eps22")
    plt.xlabel("t")
    plt.ylabel("prescribed strain")
    plt.legend()
    plt.title("Load history over time")
    plt.show()


def plot_stress_response(history: History):
    t, _, _, sigma, _, _, _, _ = history.to_arrays()

    plt.figure()
    plt.plot(t, sigma[:, 0, 0], label="sigma11")
    plt.plot(t, sigma[:, 1, 1], label="sigma22")
    plt.xlabel("t")
    plt.ylabel("stress")
    plt.legend()
    plt.title("Stress response over time")
    plt.show()


def plot_plastic_strains(history: History):
    t, _, eps_p, _, _, _, _, _ = history.to_arrays()

    plt.figure()
    plt.plot(t, eps_p[:, 0, 0], label="eps_p11")
    plt.plot(t, eps_p[:, 1, 1], label="eps_p22")
    plt.plot(t, eps_p[:, 2, 2], label="eps_p33")
    plt.xlabel("t")
    plt.ylabel("plastic strain")
    plt.legend()
    plt.title("Plastic strains over time")
    plt.show()


def plot_stress_strain(history: History):
    t, eps, _, sigma, _, _, _, _ = history.to_arrays()

    plt.figure()
    plt.plot(eps[:, 0, 0], sigma[:, 0, 0], label="sigma11 vs eps11")
    plt.plot(eps[:, 1, 1], sigma[:, 1, 1], label="sigma22 vs eps22")
    plt.xlabel("strain")
    plt.ylabel("stress")
    plt.legend()
    plt.title("Stress-strain response")
    plt.show()




def plot_lateral_contraction(history: History):
    t, eps, _, _, _, _, _, _ = history.to_arrays()

    plt.figure()
    plt.plot(t, eps[:, 2, 2], label="eps33")
    plt.xlabel("t")
    plt.ylabel("eps33")
    plt.legend()
    plt.title("Lateral contraction over time")
    plt.show()
