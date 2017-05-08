import matplotlib.pyplot as plt
import numpy as np

def get_density_profile(x, t, c_50, c_p, h, epsilon, n_max, v, P_max):
    
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    n = n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x >v*t  
    n[out_of_bound_indices] = 0.0
    return n

def vary_P(workdir):
    
    x = np.linspace(0.0, 1000.0, 1000)
    times = np.linspace(0.0, 48.0, 4)
    
    c_50 = 0.3
    c_p = 0.3
    h = 1000.0
    epsilon = 200.0
    n_max = 1.0
    v = 20.0
    P_max_vals = [0.1, 0.2, 0.5, 1.0]
    colors = ["black", "red", "green", "blue"]
    
    _, ax = plt.subplots()
    for idx, P_max in enumerate(P_max_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_density_profile(x, eachTime, c_50, c_p, h, epsilon, n_max, v, P_max)])   
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])
            
def vary_c_50(workdir):
    
    x = np.linspace(0.0, 1000.0, 1000)
    times = np.linspace(0.0, 48.0, 4)
    
    c_50 = 0.3
    c_p = 0.3
    h = 1000.0
    epsilon = 200.0
    n_max = 1.0
    v = 20.0
    P_max = 0.5
    c_50_vals = [0.01, 0.1, 0.5, 1.0]
    colors = ["black", "red", "green", "blue"]
    
    _, ax = plt.subplots()
    for idx, c_50 in enumerate(c_50_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_density_profile(x, eachTime, c_50, c_p, h, epsilon, n_max, v, P_max)])   
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])
            
def vary_h(workdir):
    
    x = np.linspace(-1000.0, 1000.0, 1000)
    times = np.linspace(0.0, 48.0, 4)
    
    c_50 = 0.3
    c_p = 0.3
    h = 1000.0
    epsilon = 200.0
    n_max = 1.0
    v = 20.0
    P_max = 0.5
    h_vals = [500.0, 750.0, 1000.0, 1250.0]
    colors = ["black", "red", "green", "blue"]
    
    _, ax = plt.subplots()
    for idx, h in enumerate(h_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_density_profile(x, eachTime, c_50, c_p, h, epsilon, n_max, v, P_max)])   
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])

workdir = "/home/grogan/test/"
vary_P(workdir)
vary_c_50(workdir)
vary_h(workdir)
    
plt.show()