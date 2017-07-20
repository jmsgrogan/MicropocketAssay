import pickle
import matplotlib.pyplot as plt
import numpy as np

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
        ax.plot(locations, results, color=colors[idx], lw=1)
        local_max = np.max(results)
        if local_max > max_result:
            max_result = local_max

    ylim = ax.get_ylim()
    ax.set_ylim([0, ylim[1]])
    ax.annotate('Limbus', xy=(150, max_result*0.99), color="C0")
    ax.annotate('Pellet', xy=(950, max_result*0.99), color="C3")
    fig.savefig(work_dir+"/" + eachParam.name + ".png",
                bbox_inches='tight', dpi=90)


def plot_domain_summary(work_dir, output_params, results, domain_types):

    output_vars = [{"name": "location_min",
                    "axis": "Front Position (um)"},
                   {"name": "location_mid",
                    "axis": "Mid Position (um)"},
                   {"name": "location_max",
                    "axis": "Max Position (um)"},
                   {"name": "density_max",
                    "axis": "Max Density (um^-2)"}, ]

    for eachParam in output_params:
        for eachVar in output_vars:
            fig, ax = plt.subplots()
            ax.set_xlabel("Time (hr)")
            ax.set_ylabel(eachVar["axis"])
            for eachDomainType in domain_types:
                print results[eachDomainType][eachParam]["output"]
                time_series = results[eachDomainType][eachParam]["output"][0]
                times = []
                results = []
                for eachTime in time_series:
                    times.append(eachTime["time"])
                    results.append(eachTime[eachVar["name"]])
                ax.plot(np.array(times), np.array(results), lw=1, label=eachDomainType)
            ax.legend(loc='upper center')
            fig.savefig(work_dir + eachVar["name"] + "_" + eachParam.name + ".png",
                        bbox_inches='tight', dpi=300)


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
                result_left_offset = 2
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
                    local_results["results"] = result[result_left_offset:]
                    local_results["locations"] = locations[result_left_offset:]
                    local_results["density_max"] = np.max(result)
                    max_arg = np.argmax(result)
                    local_results["location_max"] = max_arg
                    local_results["location_mid"] = (np.abs(result-result[max_arg]/2.0)).argmin()
                    local_results["location_min"] = (np.abs(result-0.01*result[max_arg])).argmin()
                    local_time_series.append(local_results)
                results[eachDomainType][eachParam]["output"].append(local_time_series)
    return results


def run(work_dir, domain_types, output_params,
        num_repeats, study_name):

    results = load_results(work_dir + "/" + study_name + "/",
                           domain_types, output_params, num_repeats)

    # Plot single profiles
    for eachDomainType in domain_types:
        print "Single Profiles for: " + eachDomainType
        domain_string = "/DomainType_" + eachDomainType.replace(" ", "")
        for idx in range(num_repeats):
            run_string = "/Run_" + str(idx) + "/"
            for eachParam in output_params:
                local_work_dir = work_dir + "/" + study_name + domain_string + run_string
                time_series = results[eachDomainType][eachParam]["output"][idx]
                plot_single_profile(local_work_dir, time_series, eachParam)

    #plot_domain_summary(work_dir + "/" + study_name + "/",
    #                    output_params, results, domain_types)


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
                                     title=r"Branch Density - $\mu m^{-2}$",
                                     limits=limits,
                                     line_color="C0",
                                     sampling_frequency=sampling_frequency),
                     OutputParameter(name="PDE",
                                     title=r"Concentration - nanomolar",
                                     limits=limits,
                                     line_color="C0",
                                     sampling_frequency=sampling_frequency), ]

    work_dir = "Python/Cornea/Study_fixed_gradient_d3d211d0-0887-43f6-b380-78fb09700f13"
    file_handler = chaste.core.OutputFileHandler(work_dir, False)
    work_dir = file_handler.GetOutputDirectoryFullPath()

    f = open(work_dir + "study_data.p", 'r')
    study_data = pickle.load(f)

    domain_types = study_data["domain_types"]
    num_repeats = len(study_data["random_seeds"])
    study_names = study_data["study_names"]
    for eachStudy in study_names:
        print "Processing Study: " + eachStudy
        run(work_dir, domain_types, output_params,
            num_repeats=num_repeats, study_name=eachStudy)
        break
