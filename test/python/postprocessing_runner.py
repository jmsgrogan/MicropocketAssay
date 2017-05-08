import matplotlib.pyplot as plt
import numpy as np
from symfit import Parameter, variables, Fit, exp

import chaste.core
import analytical_solutions.network_density
import postprocessing.network_density

def run(work_dir, domain_types, num_repeats):
    
    # Get the density profiles
    for eachDomainType in domain_types:
        sample_number = 0
        results_dir = work_dir + "/" + eachDomainType.replace(" ", "")+"/Run" + str(sample_number)+"/"

        line_density_locations, line_densities = postprocessing.network_density.process_csv(results_dir+"Sampled_Line_density.txt")
        tip_density_locations, tip_densities = postprocessing.network_density.process_csv(results_dir+"Sampled_Tip_density.txt")
        
        # Plot the line density 
        fig, ax = plt.subplots()
        #ax.set_ylim([0, 0.04])
        ax.set_xlim([0, 1000])
        sampling_freq = 20
        for idx, eachResult in enumerate(line_densities[::sampling_freq]):
            ax.plot(line_density_locations, np.array(eachResult[1]), color='black')

#         # Fit the line density
        cumulative_x = []
        cumulative_t = []
        cumulative_z = []
        for eachResult in line_densities:
            cumulative_x.extend(list(line_density_locations))
            cumulative_t.extend(list(np.ones(len(eachResult[1]))*eachResult[0]))
            cumulative_z.extend(eachResult[1])
#             break
#             
#         # Need to remove 0 valued regions
        Z, X, T = variables('Z, X, T')
        v = Parameter(value = 1.0)
        n = Parameter(value = 1/40.0)
        #p = Parameter(value = 0.1, min = 0.0, max = 100.0)
#       model = {Z: n*(1.0 - exp(-p*(T-X/v)))}
        model = {Z: n*(1.0 - v*T*(X-200.0)/1000.0)}

        fit = Fit(model, X=np.array(cumulative_x), T=np.array(cumulative_t), Z=np.array(cumulative_z))
        fit_result = fit.execute()
        n_fit = fit_result.value(n)
        v_fit = fit_result.value(v)
        
        print(fit_result)

#         # Plot the tip density 
#         fig, ax = plt.subplots()
#         #ax.set_ylim([0, 0.04])
#         ax.set_xlim([0, 1000])
#         sampling_freq = 20
#         for idx, eachResult in enumerate(tip_densities[::sampling_freq]):
#             ax.plot(tip_density_locations, np.array(eachResult[1]), color='black')            
        plt.show()

if __name__ == '__main__':

    file_handler = chaste.core.OutputFileHandler("Python/Cornea/TestSimulation/", False)
    work_dir = file_handler.GetOutputDirectoryFullPath()

    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    domain_types = ["Planar 2D"]
    num_repeats = 1

    run(work_dir, domain_types, num_repeats)
