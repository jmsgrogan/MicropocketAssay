import numpy as np
from symfit import Parameter, variables, Fit, exp


def fit(locations, values):

    # Fit the line density
    cumulative_x = []
    cumulative_t = []
    cumulative_z = []
    for eachResult in values:
        cumulative_x.extend(list(locations))
        cumulative_t.extend(list(np.ones(len(eachResult[1]))*eachResult[0]))
        cumulative_z.extend(eachResult[1])

        Z, X, T = variables('Z, X, T')
        v = Parameter(value=1.0)
        n = Parameter(value=1/40.0)
        #p = Parameter(value = 0.1, min = 0.0, max = 100.0)
#       model = {Z: n*(1.0 - exp(-p*(T-X/v)))}
        model = {Z: n*(1.0 - v*T*(X-200.0)/1000.0)}

        fit = Fit(model, X=np.array(cumulative_x), T=np.array(cumulative_t), Z=np.array(cumulative_z))
        fit_result = fit.execute()
        n_fit = fit_result.value(n)
        v_fit = fit_result.value(v)

        print(fit_result)

    return n_fit, v_fit
