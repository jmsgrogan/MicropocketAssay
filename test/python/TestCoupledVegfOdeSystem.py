import numpy as np
import scipy.special
import matplotlib.pyplot as plt

def get_analytical_solution(x, times, diffusivity, permeability, c_0):
    
    solutions = []
    for eachTime in times:
        h = permeability/diffusivity
        L = np.sqrt(eachTime*diffusivity)
        erfc_term_1 = scipy.special.erfc(x/(2.0*L))
        erfc_term_2 = scipy.special.erfc(x/(2.0*L) + h*L)
        exp_term = np.exp(h*x+h*h*diffusivity*eachTime)
        solutions.append(c_0*(erfc_term_1-exp_term*erfc_term_2))
    return solutions

c_0 = 1.0
D = 1.0
permeability = 1.0
times = np.linspace(0.01, 100, 10)
x = np.linspace(0.01,100,101)

solutions = get_analytical_solution(x, times, D, permeability, c_0)

fig, ax = plt.subplots()
ax.set_ylim([0,1])
ax.set_xlim([0,50])

for eachSoln in solutions:
    ax.plot(x, eachSoln, color='blue')

plt.show()


