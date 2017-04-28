import matplotlib.pyplot as plt
import numpy as np

def get_density_profile(x, t, c_50, c_p, h, epsilon, n_max, v, P_max, beta):
    
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    n = np.exp(-(beta/v)*x)*n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x >v*t  
    n[out_of_bound_indices] = 0.0
    return n

def do_plot(workdir):
    
    x = np.linspace(0.0, 1000.0, 1000)
    times = np.linspace(0.0, 48.0, 4)
    
    c_50 = 0.3
    c_p = 0.3
    h = 1000.0
    epsilon = 200.0
    n_max = (1.0/40.0)
    v = 20.0
    P_max = 0.5
    beta = 0.04
    colors = ["black", "red", "green", "blue"]
    
    _, ax = plt.subplots()

    results = [] 
    for eachTime in times:
        results.append([eachTime, get_density_profile(x, eachTime, c_50, c_p, h, epsilon, n_max, v, P_max, beta)])   
    for eachResult in results:
        ax.plot(x, eachResult[1], color='black')
       
    results = [] 
    for eachTime in times:
        results.append([eachTime, get_density_profile(x, eachTime, c_50, c_p, h, epsilon, n_max, v, P_max, 0.0)])   
    for eachResult in results:
        ax.plot(x, eachResult[1], color=colors[1])
            
workdir = "/home/grogan/test/"
do_plot(workdir)
    
plt.show()