import numpy as np
import matplotlib.pyplot as plt

def simple_evaluate_fit(t, x, v, rho_max, k, a, alpha, beta):
    
    x_0 = v*t
    z = 1.0/(1.0 + np.exp(k*(x-x_0))) - (0.4)/(1.0+np.exp(k*(x-x_0+0.7)))
    z = z*rho_max*(t/(a+t))    
    return z

x_vals = np.linspace(0.0, 2.0, 100)
t_vals = np.linspace(0.0, 12.0, 6)

rho_max = 1.0
k = 10.0
beta = 1.0
alpha = 0.0
v = 0.1
ap = 24.0

for eachTime in t_vals:
    z = simple_evaluate_fit(eachTime, x_vals, v, rho_max, k, ap, alpha, beta)
    plt.plot(x_vals, z)
    
plt.show()
    

