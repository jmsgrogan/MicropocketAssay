import matplotlib.pyplot as plt
import numpy as np
import csv

workdir = "/home/grogan/test/"


concentration_samples = []
locations = []

case_dir = workdir+"TestTwoDimensionalDomainCircular/TestTransportOnly/sample_values.txt"
#case_dir = workdir+"TestTwoDimensionalDomain/TestTransportOnly/sample_values.txt"
#case_dir = workdir+"TestOneDimensionalDomain/TestTransportOnly/sample_values.txt"
with open(case_dir, 'rb') as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    current_day = -1
    for idx, eachrow in enumerate(csv_reader):
        if idx>0:
            file_day = int(round(float(eachrow[0])/24.0))
            if file_day>current_day:
                print file_day
                concentration_samples.append([float(i) for i in eachrow[1:-1]])
                current_day = file_day
        else:
            locations = np.array([float(i) for i in eachrow[1:-1]])
            
fig, ax = plt.subplots()
ax.set_ylim([0,0.8])
ax.set_xlim([-1, 1])

locations = (2.0*locations)/np.max(locations) - 1.0

ax.plot(locations, concentration_samples[1], color='black')
ax.plot(locations, concentration_samples[2], color='black')
ax.plot(locations, concentration_samples[3], color='black')
ax.plot(locations, concentration_samples[4], color='black')
ax.plot(locations, concentration_samples[5], color='black')
#ax.plot(locations, concentration_samples[6], color='blue')
plt.show()