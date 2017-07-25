import numpy as np
from scipy.signal import convolve

a = np.array([0.0, 0.0, 1.0, 1.0, 0.0])

out = convolve(a, np.array([1,2,1]), 'same')/4.0
       
print out