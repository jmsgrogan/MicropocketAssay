import matplotlib.pyplot as plt
import numpy as np
import csv

import analytical_solutions.network_density
import postprocessing.network_density

def get_density_profile(x, t, c_50, c_p, h, epsilon, n_max, v, P_max):
    
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    n = n_max*(1.0-np.exp(-P*(t-x/v)))
    out_of_bound_indices = x >v*t  
    n[out_of_bound_indices] = 0.0
    return n

workdir = "/home/grogan/test/"

case_dir = workdir+"TestOneDimensionalDomain/TestVesselOnly/sampled_line_density.txt"
results_1d, locations_1d = postprocessing.network_density.process_csv()
with open(case_dir, 'rb') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    
    for idx, eachrow in enumerate(csv_reader):
        
        if idx==0:
            locations_1d = eachrow[1:-1]
            locations_1d = [float(x) for x in locations_1d]
        else:
            time = float(eachrow[0])
            samples = eachrow[1:-1]
            samples = [float(x) for x in samples]
            results_1d.append([time, samples])
            
case_dir = workdir+"TestTwoDimensionalDomain/TestVesselOnly/sampled_line_density.txt"
results_2d = []
locations_2d = []
with open(case_dir, 'rb') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    
    for idx, eachrow in enumerate(csv_reader):
        
        if idx==0:
            locations_2d = eachrow[1:-1]
            locations_2d = [float(x) for x in locations_2d]
        else:
            time = float(eachrow[0])
            samples = eachrow[1:-1]
            samples = [float(x) for x in samples]
            results_2d.append([time, samples])    
            
fig, ax = plt.subplots()
#ax.set_ylim([0, 0.04])
ax.set_xlim([0, 1000])

x = np.linspace(0.0, 1000.0, 1000)

c_50 = 0.5
c_p = 0.3
h = 1000.0
epsilon = 200.0
n_max = (1.0/40.0)
v = 20.0
P_max = 0.5

sampling_freq = 20
for idx, eachResult in enumerate(results_1d[::sampling_freq]):
    analytical = analytical_solutions.network_density.get_density_profile(x, eachResult[0], c_50, c_p, h, epsilon, n_max, v, P_max)
    ax.plot(x, analytical, color='red')
    ax.plot(locations_1d, np.array(eachResult[1]), color='black')
    ax.plot(locations_2d, np.array(results_2d[idx*sampling_freq][1]), color='green')
plt.show()