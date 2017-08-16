import numpy as np
import scipy.special
from microvessel_chaste.utility import *


def get_P_sprout(pc, c_p=None):

    c_50 = 0.65  # nM
    if c_p is None:
        c_p = pc.get_parameter("PelletConcentration").value.Convert(1.0e-6*mole_per_metre_cubed)
    print c_p

    h = pc.get_parameter("PelletHeight").value.Convert(1.0*1.e-6*metres)
    epsilon = pc.get_parameter("LimbalOffset").value.Convert(1.0*1.e-6*metres)
    P_max = pc.get_parameter("SproutingProbability").value.Convert((1.0/3600.0)*per_second)
    print P_max
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    return P


def get_tip_density_transient(x, t, pc, c_p=None):

    P = get_P_sprout(pc, c_p)
    v = pc.get_parameter("TipVelocity").value.Convert((1.e-6/3600.0)*metre_per_second)
    # n_max = 1.0/40.0
    n_max = 1.0/(20.0*30.0*30.0)
    L_0 = 350.0
    n = (n_max/(1.0-v/(P*L_0)))*(np.exp(-v*(t-x/v)/L_0)-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x > v * t
    n[out_of_bound_indices] = 0.0
    return n


def get_tip_density_low_velocity(x, t, pc, c_p=None):

    P = get_P_sprout(pc, c_p)
    print P, t
    v = pc.get_parameter("TipVelocity").value.Convert((1.e-6/3600.0)*metre_per_second)
    n_max = 1.0/(20.0*100.0)
    n = n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x > v * t
    n[out_of_bound_indices] = 0.0
    return n


def get_tip_density_high_velocity(x, t, pc, c_p=None):

    x = x - 100.0
    P = get_P_sprout(pc, c_p)
    v = pc.get_parameter("TipVelocity").value.Convert((1.e-6/3600.0)*metre_per_second)
    n_max = 1.0/(20.0*100.0)
    n = n_max*np.exp(-P*(t-x/v))
    print "comps", x, v*t
    out_of_bound_indices = x > (v * t)
    n[out_of_bound_indices] = 0.0
    return n

def get_line_density_high_velocity(x, t, pc, c_p=None):

    x = x - 100.0
    P = get_P_sprout(pc, c_p)
    v = pc.get_parameter("TipVelocity").value.Convert((1.e-6/3600.0)*metre_per_second)
    n_max = 1.0/(20.0*100.0)
    print "vel", v
    rho = (v/P)*n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x > (v * t)
    rho[0] += 0.01
    rho[0] /= 2.0
    rho[out_of_bound_indices] = 0.0
    return rho


def get_pde_solution(x, times, diffusivity, permeability, c_0):
    solutions = []
    for eachTime in times:
        h = permeability/diffusivity
        L = np.sqrt(eachTime*diffusivity)
        erfc_term_1 = scipy.special.erfc(x/(2.0*L))
        erfc_term_2 = scipy.special.erfc(x/(2.0*L) + h*L)
        exp_term = np.exp(h*x+h*h*diffusivity*eachTime)
        solutions.append(c_0*(erfc_term_1-exp_term*erfc_term_2))
    return solutions


# c_0 = 1.0
# D = 1.0
# permeability = 1.0
# times = np.linspace(0.01, 100, 10)
# x = np.linspace(0.01, 100, 101)
# 
# solutions = get_pde_solution(x, times, D, permeability, c_0)
# 
# fig, ax = plt.subplots()
# ax.set_ylim([0,1])
# ax.set_xlim([0,50])
# 
# for eachSoln in solutions:
#     ax.plot(x, eachSoln, color='blue')
# 
# plt.show()