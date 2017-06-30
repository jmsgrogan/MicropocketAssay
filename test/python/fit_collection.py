import numpy as np
import matplotlib.pyplot as plt

def simple_evaluate_fit(t, x, v, rho_max, k, a, alpha):
    
    x_0 = v*t
    z = 1.0/(1.0 + np.exp(k*(x-x_0))) - alpha/(1.0+np.exp(k*(x-x_0+0.3)))
    z = z*rho_max*(t/(a+t))    
    return z

x_vals = np.linspace(0.0, 1.0, 100)
t_vals = np.linspace(0.0, 1.0, 6)

rho_max = 1.0
k = 13.0
alpha = 0.0
v = 1.17
ap = 0.296

for eachTime in t_vals:
    z = simple_evaluate_fit(eachTime, x_vals, v, rho_max, k, ap, alpha)
    plt.plot(x_vals, z)
    
plt.show()
    

