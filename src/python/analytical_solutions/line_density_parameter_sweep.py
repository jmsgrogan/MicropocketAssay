import matplotlib.pyplot as plt
import numpy as np


def get_density_profile(x, t, c_50, c_p, h, epsilon, n_max, v, P_max):

    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    #n = n_max*(1.0-np.exp(-P*(t - x/v)))
    n_0 = 0.025
    v = 20.0
    n = n_max*(v/P)*(1-np.exp((P/v)*(x-v*t)))
    out_of_bound_indices = x > v*t
    n[out_of_bound_indices] = 0.0
    return n


def get_tip_density_profile(x, t, c_50, c_p, h, epsilon, n_max, v, P_max):

    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    h_ref = 500.0
    print P, v, v/(P*h_ref)
    #n = (n_max/(1.0-v/(P*h_ref)))*(np.exp(-v*(t-x/v)/h_ref)-np.exp(-P*(t-x/v)))
    n = n_max*np.exp(-P*(t-x/v))
    out_of_bound_indices = x >v*t
    
    #n[out_of_bound_indices] = 0.0
    return n

def vary_P(workdir):

    x = np.linspace(0.0, 1000.0, 1000)
    times = np.linspace(0.0, 50.0, 4)

    c_50 = 0.3
    c_p = 0.3
    h = 1000.0
    epsilon = 200.0
    n_max = 1.0/40.0
    v = 20.0
    P_max_vals = [0.01, 0.02, 0.05, 1.0]
    colors = ["black", "red", "green", "blue"]

    _, ax = plt.subplots()
    for idx, P_max in enumerate(P_max_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_density_profile(x, eachTime, c_50,
                                                          c_p, h, epsilon,
                                                          n_max, v, P_max)])
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])


def vary_c_p(workdir):

    x = np.linspace(0.0, 1000.0, 1000)
    times = np.linspace(0.0, 50.0, 5)

    c_50 = 0.5
    c_p = 0.3
    h = 1000.0
    epsilon = 100.0
    n_max = 1.0/40.0
    v = 20.0
    P_max = 0.1
    c_p_vals = [0.1]
    colors = ["black", "red", "green", "blue", "black"]

    _, ax = plt.subplots()
    for idx, c_p in enumerate(c_p_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_density_profile(x, eachTime, c_50,
                                                          c_p, h, epsilon,
                                                          n_max, v, P_max)])
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])


def vary_c_p_tip(workdir):

    x = np.linspace(0.0, 1000.0, 1000)
    times = np.linspace(0.0, 50.0, 8)

    c_50 = 0.5
    c_p = 0.3
    h = 1000.0
    epsilon = 100.0
    n_max = 1.0/40.0
    v = 20.0
    P_max = 2.0
    c_p_vals = [0.1]
    colors = ["black", "red", "green", "blue", "black"]

    _, ax = plt.subplots()
    for idx, c_p in enumerate(c_p_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_tip_density_profile(x, eachTime, c_50,
                                                          c_p, h, epsilon,
                                                          n_max, v, P_max)])
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])
            
def vary_c_p_tip_x0(workdir):

    x = 0.0
    times = np.linspace(0.0, 100.0, 80)

    c_50 = 0.5
    c_p = 0.3
    h = 1000.0
    epsilon = 100.0
    n_max = 1.0/40.0
    v = 20.0
    P_max = 2.0
    c_p_vals = [0.1, 0.4, 0.8, 1.2, 1.6]
    colors = ["black", "red", "green", "blue", "purple"]

    _, ax = plt.subplots()
    for idx, c_p in enumerate(c_p_vals):
        results = []
        for eachTime in times:
            results.append(get_tip_density_profile(x, eachTime, c_50,
                                                          c_p, h, epsilon,
                                                          n_max, v, P_max))
        ax.plot(times, np.array(results), color=colors[idx])


def vary_h(workdir):

    x = np.linspace(-1000.0, 1000.0, 1000)
    times = np.linspace(0.0, 48.0, 4)

    c_50 = 0.3
    c_p = 0.3
    h = 1000.0
    epsilon = 200.0
    n_max = 1.0
    v = 20.0
    P_max = 2.0
    h_vals = [500.0, 750.0, 1000.0, 1250.0]
    colors = ["black", "red", "green", "blue"]

    _, ax = plt.subplots()
    for idx, h in enumerate(h_vals):
        results = []
        for eachTime in times:
            results.append([eachTime, get_density_profile(x, eachTime, c_50,
                                                          c_p, h, epsilon,
                                                          n_max, v, P_max)])
        for eachResult in results:
            ax.plot(x, eachResult[1], color=colors[idx])


def plot_max_density_variation(workdir):

    times = np.linspace(0.0, 50.0, 10) # hours

    c_50 = 0.5
    c_p = 0.1
    h = 5000.0
    epsilon = 200.0
    n_max = 1.0/40.0
    v = 20.0
    P_max = 0.5*(1.0/(1.0 + (0.65/2.0)*11.0))

    _, ax = plt.subplots()
    vals = 0.025*(1.0 - np.exp(-P_max*times))
    ax.plot(times, vals)


workdir = "/home/grogan/test/"
#plot_max_density_variation(workdir)
#vary_P(workdir)
#vary_c_p(workdir)
#vary_c_p_tip(workdir)
vary_c_p_tip_x0(workdir)
#vary_h(workdir)

plt.show()