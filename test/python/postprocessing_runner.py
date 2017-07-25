from scipy.signal import convolve
import pickle
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import chaste.core
from microvessel_chaste.utility import *
from cornea.postprocessing.sampled_quantity import process_csv
from cornea.parameters.parameter_collection import SimulationParameterCollection


class OutputParameter():

    def __init__(self, name, title, limits=None, line_color="black",
                 sampling_frequency=10):
        self.name = name
        self.title = title
        self.limits = limits
        self.line_color = line_color
        self.sampling_frequency = sampling_frequency


def merge_images_x(output_path, paths):

    images = map(Image.open, paths)
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]
    new_im.save(output_path)


def merge_images_y(output_path, paths):

    images = map(Image.open, paths)
    widths, heights = zip(*(i.size for i in images))
    max_width = max(widths)
    total_height = sum(heights)

    new_im = Image.new('RGB', (max_width, total_height))
    y_offset = 0
    for im in images:
        new_im.paste(im, (0, y_offset))
        y_offset += im.size[1]
    new_im.save(output_path)


def plot_single_profile(work_dir, time_series, eachParam,
                        x_title=r"Position - $\mu m$"):

    left_line_props = {"xloc": 100.0, "color": 'C0'}
    right_line_props = {"xloc": 1100.0, "color": 'C3'}

    fig, ax = plt.subplots()
    if eachParam.limits is not None:
        if eachParam.limits[0] is not None:
            ax.set_xlim(eachParam.limits[0])
        if eachParam.limits[1] is not None:
            ax.set_ylim(eachParam.limits[1])
    ax.axvline(left_line_props['xloc'], color=left_line_props['color'],
               linestyle='--', lw=1)
    ax.axvline(right_line_props['xloc'], color=right_line_props['color'],
               linestyle='--', lw=1)

    ax.set_xlabel(x_title)
    ax.set_ylabel(eachParam.title)
    max_result = 0.0
    colormap = plt.cm.viridis
    colorscale = np.linspace(0, 1, len(time_series))
    colors = [colormap(i) for i in colorscale]

    for idx, eachTime in enumerate(time_series):
        results = eachTime["results"]
        locations = eachTime["locations"]
        smooth_result = convolve(np.array(results), np.array([1,2,1]), 'same')/4.0
        smooth_result[0] = results[0]
        smooth_result[-1] = results[-1]
        ax.plot(locations, smooth_result, color=colors[idx], lw=1)
        local_max = np.max(smooth_result)
        if local_max > max_result:
            max_result = local_max

    ylim = ax.get_ylim()
    ax.set_ylim([0, ylim[1]])
    ax.annotate('Limbus', xy=(150, max_result*0.99), color="C0")
    ax.annotate('Pellet', xy=(950, max_result*0.99), color="C3")
    fig.savefig(work_dir+"/" + eachParam.name + ".png", bbox_inches='tight', dpi=90)


def plot_domain_summary(work_dir, output_params, results, domain_types, num_repeats):

    output_vars = [{"name": "location_min",
                    "axis": "Front Position (um)"},
                   {"name": "location_mid",
                    "axis": "Mid Position (um)"},
                   {"name": "location_max",
                    "axis": "Max Position (um)"},
                   {"name": "density_max",
                    "axis": "Max Density (um^-2)"}, ]

    outer_plot_list = []
    summary_results = {}
    
    for eachParam in output_params:
        summary_results[eachParam.name] = {}
        plot_list = []
        for eachVar in output_vars:
            summary_results[eachParam.name][eachVar["name"]] = {}
            fig, ax = plt.subplots()
            ax.set_xlabel("Time (hr)")
            ax.set_ylabel(eachVar["axis"])
            if "location" in eachVar["name"]:
                ax.set_ylim([0, 1100])
            for eachDomainType in domain_types:
                summary_results[eachParam.name][eachVar["name"]][eachDomainType] = {"midpoints": [],
                                                                       "slopes": []}
                for idx in range(num_repeats):
                    time_series = results[eachDomainType][eachParam]["output"][idx]
                    times = []
                    local_results = []
                    for eachTime in time_series:
                        times.append(eachTime["time"])
                        local_results.append(eachTime[eachVar["name"]])
    
                    # Get midpoint and slope
                    dt = times[1] - times[0]
                    smooth_result = convolve(np.array(local_results), np.array([1,2,1]), 'same')/4.0
                    smooth_result[0] = local_results[0]
                    smooth_result[-1] = local_results[-1]
                    mid_index = int(round(len(smooth_result)/2.0))
                    md_val = smooth_result[mid_index]
                    lower_index = int(round(len(smooth_result)/4.0))
                    upper_index = int(round(3.0*len(smooth_result)/4.0))
                    lower_val = smooth_result[lower_index]
                    lower_time = float(lower_index)*dt
                    upper_val = smooth_result[upper_index]
                    upper_time = float(upper_index)*dt
                    slope = (upper_val-lower_val)/(upper_time-lower_time)
                    summary_results[eachParam.name][eachVar["name"]][eachDomainType]["slopes"].append(slope)
                    summary_results[eachParam.name][eachVar["name"]][eachDomainType]["midpoints"].append(md_val)
                    if idx == 0:
                        ax.plot(np.array(times), smooth_result, lw=1, label=eachDomainType)
            ax.legend(loc='upper center')
            plot_path = work_dir + eachVar["name"] + "_" + eachParam.name + ".png"
            plot_list.append(plot_path)
            fig.savefig(plot_path, bbox_inches='tight', dpi=90)
        merge_path = work_dir + eachParam.name + "_merge.png"
        outer_plot_list.append(merge_path)
        merge_images_x(merge_path, plot_list)
    merge_path = work_dir + "merge.png"
    merge_images_y(merge_path, outer_plot_list)
    return summary_results


def domain_is_3d(domain_type):
    return ("3D" in domain_type or "Hemi" in domain_type)


def load_results(work_dir, domain_types, output_params, num_repeats):
    results = {}
    for eachDomainType in domain_types:
        results[eachDomainType] = {}
        domain_string = "/DomainType_" + eachDomainType.replace(" ", "")
        for eachParam in output_params:
            results[eachDomainType][eachParam] = {"params": [],
                                                  "output": []}
        for idx in range(num_repeats):
            run_string = "/Run_" + str(idx)
            pc_dir = work_dir + domain_string + run_string
            pc = SimulationParameterCollection()
            pc.load(pc_dir + "/input_parameters.p")

            for eachParam in output_params:
                results[eachDomainType][eachParam]["params"].append(pc)

                results_dir = pc_dir + "/Sampled_" + eachParam.name
                locations, values = process_csv(results_dir + ".txt")

                thickness = pc.get_parameter("CorneaThickness").value
                thickness = thickness.Convert(1.e-6*metres)
                sampled_values = values[::eachParam.sampling_frequency]
                result_left_offset = 0
                local_time_series = []
                for idx, eachResult in enumerate(sampled_values):
                    local_results = {"location_min": 0.0,
                                     "location_mid": 0.0,
                                     "location_max": 0.0,
                                     "density_max": 0.0,
                                     "locations": [],
                                     "results": [],
                                     "time": eachResult[0]}  # hr
                    result = np.array(eachResult[1])
                    if "density" in eachParam.name and domain_is_3d(eachDomainType):
                        result *= thickness
                    offset_result = result[result_left_offset:]
                    offset_locations = locations[result_left_offset:]
                    local_results["results"] = offset_result
                    local_results["locations"] = offset_locations
                    local_results["density_max"] = np.max(offset_result)
                    max_arg = np.argmax(offset_result)
                    local_results["location_max"] = offset_locations[max_arg]
                    local_results["location_mid"] = offset_locations[(np.abs(offset_result-offset_result[max_arg]/2.0)).argmin()]
                    local_results["location_min"] = offset_locations[(np.abs(offset_result-0.01*offset_result[max_arg])).argmin()]
                    local_time_series.append(local_results)
                results[eachDomainType][eachParam]["output"].append(local_time_series)
    return results


def plot_top_level_domain_summary(work_dir, summary_results, domain_types,
                                  output_params, results, num_repeats):

    names  = ["Line", "Tip", "Branch"]
    for eachName in names:
        fig, ax = plt.subplots()
        ax.set_ylabel(eachName + " Velocity (um/hr)")
        ind = np.arange(len(domain_types))
        width = 0.3
        ax.set_xticks(ind)
        ax.set_xticklabels([x for x in domain_types])
        vals = ["min", "mid", "max"]
        colors = ["r", "g", "b"]
        for idx in range(len(vals)):
            av_vels = []
            std_vels = []
            for eachDomain in domain_types:
                velocities = summary_results[eachName + "_density"]["location_"+vals[idx]][eachDomain]["slopes"]
                av_vels.append(np.mean(np.array(velocities)))
                std_vels.append(np.std(np.array(velocities)))
            ax.bar(ind +float(idx)*width, np.array(av_vels), width, color=colors[idx], yerr=np.array(std_vels), alpha=0.3)
        ax.set_ylim([0, 25])
        fig.savefig(work_dir+"/top_level_" + eachName +"_velocities.png", bbox_inches='tight', dpi=90)

    for eachName in names:
        fig, ax = plt.subplots()
        ax.set_ylabel(eachName + " Location (um)")
        ind = np.arange(len(domain_types))
        width = 0.3
        ax.set_xticks(ind)
        ax.set_xticklabels([x for x in domain_types])
        vals = ["min", "mid", "max"]
        colors = ["r", "g", "b"]
        for idx in range(len(vals)):
            av_vels = []
            std_vels = []
            for eachDomain in domain_types:
                velocities = summary_results[eachName + "_density"]["location_"+vals[idx]][eachDomain]["midpoints"]
                av_vels.append(np.mean(np.array(velocities)))
                std_vels.append(np.std(np.array(velocities)))
            ax.bar(ind +float(idx)*width, np.array(av_vels), width, color=colors[idx], yerr=np.array(std_vels), alpha=0.3)
        ax.set_ylim([0, 1000])
        fig.savefig(work_dir+"/top_level_" + eachName +"_locations.png", bbox_inches='tight', dpi=90)
        
    for eachName in names:
        fig, ax = plt.subplots()
        ax.set_ylabel(eachName + " Density Rate of Increase (um^-2hr^-1)")
        ind = np.arange(len(domain_types))
        width = 0.7
        ax.set_xticks(ind)
        ax.set_xticklabels([x for x in domain_types])
        av_vels = []
        std_vels = []
        for eachDomain in domain_types:
            velocities = summary_results[eachName + "_density"]["density_max"][eachDomain]["slopes"]
            av_vels.append(np.mean(np.array(velocities)))
            std_vels.append(np.std(np.array(velocities)))
        ax.bar(ind, np.array(av_vels), width, color="r", yerr=np.array(std_vels), alpha=0.3)
        #ax.set_ylim([0, 1000])
        fig.savefig(work_dir+"/top_level_" + eachName +"_density_rate.png", bbox_inches='tight', dpi=90)
        
    for eachName in names:
        fig, ax = plt.subplots()
        ax.set_ylabel(eachName + " Density(um^-2)")
        ind = np.arange(len(domain_types))
        width = 0.7
        ax.set_xticks(ind)
        ax.set_xticklabels([x for x in domain_types])
        av_vels = []
        std_vels = []
        for eachDomain in domain_types:
            velocities = summary_results[eachName + "_density"]["density_max"][eachDomain]["midpoints"]
            av_vels.append(np.mean(np.array(velocities)))
            std_vels.append(np.std(np.array(velocities)))
        ax.bar(ind, np.array(av_vels), width, color="r", yerr=np.array(std_vels), alpha=0.3)
        #ax.set_ylim([0, 1000])
        fig.savefig(work_dir+"/top_level_" + eachName +"_density.png", bbox_inches='tight', dpi=90)


def run(work_dir, domain_types, output_params,
        num_repeats, study_name):

    results = load_results(work_dir + "/" + study_name + "/",
                           domain_types, output_params, num_repeats)

    # Plot single profiles
    for eachDomainType in domain_types:
        print "Single Profiles for: " + eachDomainType
        domain_string = "/DomainType_" + eachDomainType.replace(" ", "")
        outer_figure_paths = []
        for idx in range(num_repeats):
            run_string = "/Run_" + str(idx) + "/"
            figure_paths = []
            for eachParam in output_params:
                local_work_dir = work_dir + "/" + study_name + domain_string + run_string
                figure_paths.append(local_work_dir + "/" + eachParam.name + ".png")
                time_series = results[eachDomainType][eachParam]["output"][idx]
                plot_single_profile(local_work_dir, time_series, eachParam)
            merge_path = work_dir + "/" + study_name + domain_string + run_string + "merge.png"
            outer_figure_paths.append(merge_path)
            merge_images_x(merge_path, figure_paths)
        merge_path = work_dir + "/" + study_name + domain_string + "/merge.png" 
        merge_images_y(merge_path, outer_figure_paths)

    summary_results = plot_domain_summary(work_dir + "/" + study_name + "/", 
                                          output_params, results, domain_types, num_repeats)

    plot_top_level_domain_summary(work_dir + "/" + study_name + "/", summary_results, 
                                  domain_types, output_params, results, num_repeats)


if __name__ == '__main__':

    limits = [[0, 1200], None]
    sampling_frequency = 2

    output_params = [OutputParameter(name="Line_density",
                                     title=r"Line Density - $\mu m$ per $\mu m^2$",
                                     limits=limits,
                                     line_color="C0",
                                     sampling_frequency=sampling_frequency),
                     OutputParameter(name="Tip_density",
                                     title=r"Tip Density - $\mu m^{-2}$",
                                     limits=limits,
                                     line_color="C0",
                                     sampling_frequency=sampling_frequency),
                     OutputParameter(name="Branch_density",
                                     title=r"Branch Density - $\mu m^{-1}$",
                                     limits=limits,
                                     line_color="C0",
                                     sampling_frequency=sampling_frequency),
                     OutputParameter(name="PDE",
                                     title=r"Concentration - nanomolar",
                                     limits=limits,
                                     line_color="C0",
                                     sampling_frequency=sampling_frequency), ]

#    work_dir = "Python/Cornea/Study_fixed_gradient_e682e545-d635-427c-ac4e-a38419daf014"
    work_dir = "Python/Cornea/Study_fixed_gradient_c3e35591-035d-4eae-a83b-e40a9422af3f"
    #work_dir = "Python/Cornea/Study_pde_e477273c-fa5e-4eec-bbfe-bc2776e42654"
    work_dir = "Python/Cornea/Study_fg_vary_cpa238f435-d103-4cb3-b583-d38a94720733"
    work_dir = "Python/Cornea/Study_fg_vary_cp_low_chem42d39da3-ef5f-4431-90b7-eae88adabdc9"
    work_dir = "Python/Cornea/Study_fg_vary_cp_sprout_no_ana8012721a-f37b-42ad-98cd-dfe5bac37faf"
    file_handler = chaste.core.OutputFileHandler(work_dir, False)
    work_dir = file_handler.GetOutputDirectoryFullPath()

    f = open(work_dir + "study_data.p", 'r')
    study_data = pickle.load(f)

    domain_types = study_data["domain_types"]
    num_repeats = len(study_data["random_seeds"])
    #num_repeats = 3
    study_names = study_data["study_names"]
    count = 0
    for eachStudy in study_names:
        print "Processing Study: " + eachStudy
        run(work_dir, domain_types, output_params,
            num_repeats=num_repeats, study_name=eachStudy)
        count = count + 1
#         if count == 1:
#             break
