import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

def get_distance_circle(L, d0=40.0, R0=1300):
    return (d0/R0)*np.ones(L.size)

def get_distance_hemisphere(L, d0=40.0, R0=1300):
    return (d0/R0)*np.sin(L/R0)

def get_density_circle(L, d0=40.0, R0=1300):
    return 1.0/((1.0-L/R0)*(1.0-L/R0))

def get_density_hemisphere(L, d0=40.0, R0=1300):
    return np.tan(L/R0)*(1.0/np.cos(L/R0))/R0

R0 = 1300
locations = np.linspace(0.0, 0.5*R0, 100)
            
fig, ax = plt.subplots()
#ax.set_ylim([0,0.05])
ax.set_xlim([0, 0.5])
ax.set_xlabel(r'$\frac{L}{R_0}$')
ax.set_ylabel(r'$\frac{-d \rho}{dL}\left(\frac{1}{\rho_0}\right)$')

ax.plot(locations/R0, np.zeros(locations.size), color='red', label="Front On")
ax.plot(locations/R0, get_density_circle(locations), color='black', label="Top Down")
ax.plot(locations/R0, get_density_hemisphere(locations), color='green', label = "Hemisphere")
plt.legend()
plt.tight_layout()
plt.show()