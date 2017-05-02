"""
Analytical solutions for vessel line density
"""

import numpy as np

def get_density_profile(x, t, c_50, c_p, h, epsilon, n_max, v, P_max, beta):
    
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    n = np.exp(-(beta/v)*x)*n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x >v*t  
    n[out_of_bound_indices] = 0.0
    return n

def get_density_profile_with_sink(x, t, c_50, c_p, h, epsilon, n_max, v, P_max, beta):
    
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    n = n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x >v*t  
    n[out_of_bound_indices] = 0.0
    return n
