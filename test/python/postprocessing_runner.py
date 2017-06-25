import matplotlib.pyplot as plt
import numpy as np

import chaste.core
import cornea.postprocessing.sampled_quantity
import cornea.postprocessing.do_fitting

class OutputParameter():
    
    def __init__(self, name, title, limits=None, line_color="black", sampling_frequency=10):
        
        self.name = name
        self.title = title
        self.limits = limits
        self.line_color = line_color
        self.sampling_frequency = sampling_frequency

def make_figure(work_dir, locations, values, eachParam, x_title=r"Position - $\mu m$"):
    
    fig, ax = plt.subplots()
    
    if eachParam.limits is not None:
        if eachParam.limits[0] is not None:
            ax.set_xlim(eachParam.limits[0])
        if eachParam.limits[1] is not None:
            ax.set_ylim(eachParam.limits[1])
    
    ax.axvline(100.0, color='C0', linestyle='--', lw=1)
    ax.axvline(1100.0, color='C3', linestyle='--', lw=1)
    ax.set_xlabel(x_title) 
    ax.set_ylabel(eachParam.title) 
    
    max_result = 0.0
    colormap = plt.cm.viridis #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 1,len(values[::eachParam.sampling_frequency]))]
    
    for idx, eachResult in enumerate(values[::eachParam.sampling_frequency]):
        ax.plot(locations, np.array(eachResult[1]), color=colors[idx], lw=1)
        local_max = np.max(np.array(eachResult[1]))
        if local_max > max_result:
            max_result = local_max
            
    ylim = ax.get_ylim()
    ax.set_ylim([0, ylim[1]])
    
    ax.annotate('Limbus', xy=(150, max_result*0.99), color="C0")
    ax.annotate('Pellet', xy=(950, max_result*0.99), color="C3")
        
    fig.savefig(work_dir,bbox_inches='tight',dpi=300) 
    
def run(work_dir, domain_types, output_params, num_repeats):
    
    # Get the density profiles
    for eachDomainType in domain_types:
        sample_number = 0
        
        for eachParam in output_params:
        
            results_dir = work_dir + "/DomainType_" + eachDomainType.replace(" ", "")+"/Run_" + str(sample_number)+"/"
            results_dir += "Sampled_" + eachParam.name
            locations, values  = cornea.postprocessing.sampled_quantity.process_csv(results_dir + ".txt")

            #cornea.postprocessing.do_fitting.fit(locations, values)
            make_figure(results_dir+".png", locations, values, eachParam)

if __name__ == '__main__':

    #file_handler = chaste.core.OutputFileHandler("Python/Cornea/ParamSweep_WithConsumption/ParamName_sproutingprobability/ParamValue_0", False)
    file_handler = chaste.core.OutputFileHandler("Python/Cornea/TestSimulationFixedGradient", False)
    work_dir = file_handler.GetOutputDirectoryFullPath()

    domain_types = ["Planar_2D", "Planar_3D", "Circle_2D", "Circle_3D", "Hemisphere"]
    
    #domain_types = ["Planar_2D"]
    
    output_params = [OutputParameter(name = "Line_density", title = r"Line Density - $\mu m$ per $\mu m^3$",
                                     limits = [[0, 1200], None], line_color = "C0", sampling_frequency=2),
                     OutputParameter(name = "Tip_density", title = r"Tip Density - $\mu m^{-3}$",
                                     limits = [[0, 1200], None], line_color = "C0", sampling_frequency=2),
                     OutputParameter(name = "Branch_density", title = r"Branch Density - $\mu m^{-3}$",
                                     limits = [[0, 1200], None], line_color = "C0", sampling_frequency=2),
                     OutputParameter(name = "PDE", title = r"Concentration - nanomolar",
                                     limits = [[0, 1200], None], line_color = "C0", sampling_frequency=2),]
    
#     output_params = [OutputParameter(name = "PDE", title = r"Concentration - nanomolar",
#                                      limits = [[0, 1200], None], line_color = "C0", sampling_frequency=10),]
    
    num_repeats = 1
    run(work_dir, domain_types, output_params, num_repeats)
