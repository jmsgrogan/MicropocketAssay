import numpy as np
from symfit import parameters, variables, Fit

x = np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
t = np.array([0.0, 0.0, 1.0, 1.0, 2.0, 2.0])
z = np.array([0.0, 0.0, 1.0, 1.0, 4.0, 4.0])

Z, X, T = variables('Z, X, T')
a, b = parameters('a, b')

model = {Z: a * X + b * T**2}
fit = Fit(model, X=x, T=t, Z=z)
fit_result = fit.execute()

print(fit_result)



