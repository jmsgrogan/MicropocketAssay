import os
import matplotlib.pyplot as plt
import numpy as np

import chaste
from microvessel_chaste.utility import *
from cornea.postprocessing.plotting_tools import *


def get_analytical_solution(x, t, pc, c_p=None):

    c_50 = 0.5
    if c_p is None:
        c_p = pc.get_parameter("PelletConcentration").value.Convert(1.0e-7*mole_per_metre_cubed)
    print c_p
    h = pc.get_parameter("PelletHeight").value.Convert(1.e-6*metres)
    epsilon = pc.get_parameter("LimbalOffset").value.Convert(1.e-6*metres)
    n_max = 1.0/40.0
    v = pc.get_parameter("TipVelocity").value.Convert((1.e-6/3600.0)*metre_per_second)
    P_max = pc.get_parameter("SproutingProbability").value.Convert((1.0/3600.0)*per_second)
    x = x - epsilon

    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    
    n_max = 1.0/(20.0*30.0*30.0)
    L_0=350.0
    #n = n_max*np.exp(-P*(t-x/v))
    #n = n_max*(1.0-np.exp(-P*(t-x/v)))
    n = (n_max/(1.0-v/(P*L_0)))*(np.exp(-v*(t-x/v)/L_0)-np.exp(-P*(t-x/v)))
    
    out_of_bound_indices = x > v * t
    n[out_of_bound_indices] = 0.0
    return n


def do_density_line_plots(rc, x_title=r"Position - $\mu m$"):

    left_line_props = {"xloc": 100.0, "color": 'C0'}
    right_line_props = {"xloc": 1100.0, "color": 'C3'}
    colormap = plt.cm.viridis
    dpi = 90

    for eachStudy in rc.study_data["study_names"]:
        print "Density Line Plots for: " + eachStudy
        for eachDomain in rc.study_data["domain_types"]:
            if "Planar_2D" not in eachDomain:
                continue
            print "Density Line Plots for: " + eachDomain
            num_random = len(rc.study_data["random_seeds"])
            outer_figure_paths = []
            for idx in range(num_random):
                local_work_dir = rc.get_path(eachStudy, eachDomain, str(idx))
                figure_paths = []
                for eachParam in rc.parameters:
                    param = eachParam.name
                    if len(rc.results[eachStudy][eachDomain][param]["output"])<=idx:
                        continue
                    fig_path = local_work_dir + "/" + param + ".png"
                    figure_paths.append(fig_path)
                    time_series = rc.results[eachStudy][eachDomain][param]["output"][idx]
                    pc = rc.results[eachStudy][eachDomain][param]["params"][idx]
                    fig, ax = plt.subplots()
                    if eachParam.limits is not None:
                        if eachParam.limits[0] is not None:
                            ax.set_xlim(eachParam.limits[0])
                        if eachParam.limits[1] is not None:
                            ax.set_ylim(eachParam.limits[1])
                    height = pc.get_parameter("PelletHeight").value
                    height = height.Convert(1.e-6*metres)
                    offset = pc.get_parameter("LimbalOffset").value
                    offset = offset.Convert(1.e-6*metres)
                    xline = height + offset
                    if "Hemisphere" in eachDomain:
                        radius = pc.get_parameter("CorneaRadius").value
                        radius = radius.Convert(1.e-6*metres)
                        xline = radius*np.arcsin(xline/radius)
                    ax.axvline(left_line_props['xloc'], 
                               color=left_line_props['color'],
                               linestyle='--', lw=1)
                    ax.axvline(xline, 
                               color=right_line_props['color'],
                               linestyle='--', lw=1)
                    ax.set_xlabel(x_title)
                    ax.set_ylabel(rc.density_plot_keys[param])
                    max_result = 0.0
                    colorscale = np.linspace(0, 1, len(time_series))
                    colors = [colormap(i) for i in colorscale]
                    for jdx, eachTime in enumerate(time_series):
                        results = eachTime["results"]
                        locations = eachTime["locations"]
                        smooth_result = smooth_results(np.array(results))
                        if "Tip" in eachParam.name:
                            ana_results = get_analytical_solution(locations, 
                                                                  eachTime["time"],
                                                                  pc)
                            ax.plot(locations, ana_results, color='black', lw=1)
                        ax.plot(locations, smooth_result, color=colors[jdx], lw=1)
                        local_max = np.max(smooth_result)
                        if local_max > max_result:
                            max_result = local_max
                    ylim = ax.get_ylim()
                    ax.set_ylim([0, ylim[1]])
                    ax.annotate('Limbus', xy=(150, max_result*0.99), color="C0")
                    ax.annotate('Pellet', xy=(950, max_result*0.99), color="C3")
                    fig.savefig(fig_path, bbox_inches='tight', dpi=dpi)
                merge_path = local_work_dir + "/merge.png"
                outer_figure_paths.append(merge_path)
                merge_images_x(merge_path, figure_paths)
            merge_path = rc.get_path(eachStudy, eachDomain) + "/merge.png"
            merge_images_y(merge_path, outer_figure_paths)


def do_density_domain_comparisons(rc):

    dpi = 90
    ymax = 1100

    for eachStudy in rc.study_data["study_names"]:
        print "Density Domain Plots for: " + eachStudy
        outer_plot_list = []
        work_dir = rc.get_path(eachStudy)
        for eachParam in rc.parameters:
            plot_list = []
            output_vars = rc.location_plot_keys.keys()
            param = eachParam.name
            for eachVar in output_vars:
                fig, ax = plt.subplots()
                ax.set_xlabel("Time (hr)")
                ax.set_ylabel(rc.location_plot_keys[eachVar])
                if "location" in eachVar:
                    ax.set_ylim([0, ymax])
                for eachDomain in rc.study_data["domain_types"]:
                    num_random = len(rc.study_data["random_seeds"])
                    res = rc.results[eachStudy][eachDomain][param]["output"]
                    for idx in range(num_random):
                        if len(res)<=idx:
                            continue
                        time_series = res[idx]
                        times = []
                        local_results = []
                        for eachTime in time_series:
                            times.append(eachTime["time"])
                            local_results.append(eachTime[eachVar])
                        smooth_result = smooth_results(local_results)
                        if idx == 0:
                            ax.plot(np.array(times), smooth_result, lw=1,
                                    label=eachDomain)
                ax.legend(loc='upper center')
                plot_path = work_dir + eachVar + "_" + param + ".png"
                plot_list.append(plot_path)
                fig.savefig(plot_path, bbox_inches='tight', dpi=dpi)
            merge_path = work_dir + param + "_merge.png"
            outer_plot_list.append(merge_path)
            merge_images_x(merge_path, plot_list)
        merge_path = work_dir + "merge.png"
        merge_images_y(merge_path, outer_plot_list)


def do_density_position_summary(rc):

    names = ["Line", "Tip", "Branch"]
    vals = ["min", "mid", "max"]
    colors = ["r", "g", "b"]
    width = 0.3
    domain_types = rc.study_data["domain_types"]
    dpi = 90
    alpha = 0.3

    for eachStudy in rc.study_data["study_names"]:
        print "Density Position Summary for: ", eachStudy
        for eachName in names:
            fig, ax = plt.subplots()
            ax.set_ylabel(eachName + " Velocity (um/hr)")
            ind = np.arange(len(domain_types))
            ax.set_xticks(ind)
            ax.set_xticklabels([x for x in domain_types])
            for idx in range(len(vals)):
                av_vels = []
                std_vels = []
                var_label = "location_" + vals[idx]
                param_label = eachName + "_density"
                for eachDomain in domain_types:
                    summaries = rc.results[eachStudy][eachDomain][param_label]["summaries"]
                    velocities = summaries[var_label]["slopes"]
                    av_vels.append(np.mean(np.array(velocities)))
                    std_vels.append(np.std(np.array(velocities)))
                ax.bar(ind + float(idx)*width, np.array(av_vels), width,
                       color=colors[idx], yerr=np.array(std_vels), alpha=alpha)
            ax.set_ylim([0, 25])
            work_dir = rc.get_path(eachStudy)
            fig.savefig(work_dir+"/top_level_" + eachName + "_velocities.png",
                        bbox_inches='tight', dpi=dpi)

        for eachName in names:
            fig, ax = plt.subplots()
            ax.set_ylabel(eachName + " Location (um)")
            ind = np.arange(len(domain_types))
            ax.set_xticks(ind)
            ax.set_xticklabels([x for x in domain_types])
            for idx in range(len(vals)):
                av_vels = []
                std_vels = []
                var_label = "location_"+vals[idx]
                param_label = eachName + "_density"
                for eachDomain in domain_types:
                    summaries = rc.results[eachStudy][eachDomain][param_label]["summaries"]
                    velocities = summaries[var_label]["midpoints"]
                    av_vels.append(np.mean(np.array(velocities)))
                    std_vels.append(np.std(np.array(velocities)))
                ax.bar(ind + float(idx)*width, np.array(av_vels), width,
                       color=colors[idx], yerr=np.array(std_vels), alpha=alpha)
            ax.set_ylim([0, 1000])
            work_dir = rc.get_path(eachStudy)
            fig.savefig(work_dir+"/top_level_" + eachName + "_locations.png",
                        bbox_inches='tight', dpi=dpi)

        for eachName in names:
            fig, ax = plt.subplots()
            ax.set_ylabel(eachName + " Density Rate of Increase (um^-2hr^-1)")
            ind = np.arange(len(domain_types))
            ax.set_xticks(ind)
            ax.set_xticklabels([x for x in domain_types])
            av_vels = []
            std_vels = []
            var_label = "density_max"
            param_label = eachName + "_density"
            for eachDomain in domain_types:
                summaries = rc.results[eachStudy][eachDomain][param_label]["summaries"]
                velocities = summaries[var_label]["slopes"]
                av_vels.append(np.mean(np.array(velocities)))
                std_vels.append(np.std(np.array(velocities)))
            ax.bar(ind, np.array(av_vels), width, color="r",
                   yerr=np.array(std_vels), alpha=alpha)
            work_dir = rc.get_path(eachStudy)
            fig.savefig(work_dir+"/top_level_" + eachName + "_density_rate.png",
                        bbox_inches='tight', dpi=dpi)

        for eachName in names:
            fig, ax = plt.subplots()
            ax.set_ylabel(eachName + " Density(um^-2)")
            ind = np.arange(len(domain_types))
            ax.set_xticks(ind)
            ax.set_xticklabels([x for x in domain_types])
            av_vels = []
            std_vels = []
            var_label = "density_max"
            param_label = eachName + "_density"
            for eachDomain in domain_types:
                summaries = rc.results[eachStudy][eachDomain][param_label]["summaries"]
                velocities = summaries[var_label]["midpoints"]
                av_vels.append(np.mean(np.array(velocities)))
                std_vels.append(np.std(np.array(velocities)))
            ax.bar(ind, np.array(av_vels), width, color="r",
                   yerr=np.array(std_vels), alpha=alpha)
            work_dir = rc.get_path(eachStudy)
            fig.savefig(work_dir+"/top_level_" + eachName + "_density.png",
                        bbox_inches='tight', dpi=dpi)


def do_line_density_overview(rc):

    dpi = 90
    fig, ax = plt.subplots()
    ax.set_ylabel("Line Density (um per um^-2)")
    ax.set_xlabel("Time (hr)")
    c_p = [0.1, 0.5, 1.0, 1.5, 2.0]

    for idx, eachStudy in enumerate(rc.study_data["study_names"]):
        if len(rc.results[eachStudy]["Planar_2D"]["Line_density"]["output"]) ==0:
            continue
        time_series = rc.results[eachStudy]["Planar_2D"]["Line_density"]["output"][0]
        pc = rc.results[eachStudy]["Planar_2D"]["Line_density"]["params"][0]
        times = []
        local_results = []
        analytical = []
        for eachTime in time_series:
            times.append(eachTime["time"])
            local_results.append(eachTime["density_max"])
            full_ana = get_analytical_solution(np.array([0]), eachTime["time"],
                                               pc, c_p=c_p[idx])
            analytical.append(full_ana[0])
        smooth_result = smooth_results(local_results)
        ax.plot(np.array(times), smooth_result, lw=1, label=eachStudy)
        #ax.plot(np.array(times), analytical, lw=1, label=eachStudy, color="black")
    work_dir = rc.get_path()
    fig.savefig(work_dir+"/study_level_line_density.png",
                bbox_inches='tight', dpi=dpi)
    
def get_analytical_xo_tip(t, pc):

    c_50 = 0.5
    c_p = pc.get_parameter("PelletConcentration").value.Convert(1.0e-7*mole_per_metre_cubed)
    h = pc.get_parameter("PelletHeight").value.Convert(1.e-6*metres)
    epsilon = pc.get_parameter("LimbalOffset").value.Convert(1.e-6*metres)

    v = pc.get_parameter("TipVelocity").value.Convert((1.e-6/3600.0)*metre_per_second)
    P_max = pc.get_parameter("SproutingProbability").value.Convert((1.0/3600.0)*per_second)
    alpha = (h+epsilon)/(0+epsilon)
    P = P_max*(1.0/(1.0+(c_50/c_p)*alpha))
    
    s = 20.0
    n_max = 1.0/(s*40.0*40.0)
    L_0=500.0
    
    #n = n_max*(1.0-np.exp(-P*t))
    n = (n_max/(1.0-v/(P*L_0)))*(np.exp(-v*t/L_0)-np.exp(-P*t))

    return n

def do_tip_density_x0(rc):

    dpi = 90
    fig, ax = plt.subplots()
    ax.set_ylabel("Tip Density (per um^-3)")
    ax.set_xlabel("Time (hr)")
    c_p = [0.1, 0.5, 1.0, 1.5, 2.0]

    for idx, eachStudy in enumerate(rc.study_data["study_names"]):
        if len(rc.results[eachStudy]["Planar_2D"]["Tip_density"]["output"]) == 0:
            continue
        time_series = rc.results[eachStudy]["Planar_2D"]["Tip_density"]["output"][0]
        pc = rc.results[eachStudy]["Planar_2D"]["Tip_density"]["params"][0]
        times = []
        local_results = []
        analytical = []
        for eachTime in time_series:
            times.append(eachTime["time"])
            local_results.append(eachTime["results"][0])
            full_ana = get_analytical_xo_tip(eachTime["time"], pc)
            analytical.append(full_ana)
        smooth_result = smooth_results(local_results)
        smooth_result = smooth_results(smooth_result)
        ax.plot(np.array(times), smooth_result, lw=1, label=eachStudy)
        ax.plot(np.array(times), analytical, lw=1, label=eachStudy, color="black")
    work_dir = rc.get_path()
    fig.savefig(work_dir+"/xo_tip_density.png",
                bbox_inches='tight', dpi=dpi)


def do_front_position_overview(rc):

    dpi = 90
    fig, ax = plt.subplots()
    ax.set_ylabel("Front Position (um)")
    ax.set_xlabel("Time (hr)")
    c_p = [0.1, 0.5, 1.0, 1.5, 2.0]

#     for idx, eachStudy in enumerate(rc.study_data["study_names"]):
#         time_series = rc.results[eachStudy]["Planar_2D"]["Line_density"]["output"][0]
#         pc = rc.results[eachStudy]["Planar_2D"]["Line_density"]["params"][0]
#         times = []
#         local_results = []
#         analytical = []
#         for eachTime in time_series:
#             times.append(eachTime["time"])
#             local_results.append(eachTime["location_min"])
#         smooth_result = smooth_results(local_results)
#         ax.plot(np.array(times), smooth_result, lw=1, label=eachStudy)

    for idx, eachStudy in enumerate(rc.study_data["study_names"]):
        if len(rc.results[eachStudy]["Planar_2D"]["Tip_density"]["output"]) ==0:
            continue
        time_series = rc.results[eachStudy]["Planar_2D"]["Tip_density"]["output"][0]
        pc = rc.results[eachStudy]["Planar_2D"]["Tip_density"]["params"][0]
        times = []
        local_results = []
        analytical = []
        for eachTime in time_series:
            times.append(eachTime["time"])
            local_results.append(eachTime["location_min"])
        smooth_result = smooth_results(local_results)
        ax.plot(np.array(times), smooth_result, lw=1, label=eachStudy)
    work_dir = rc.get_path()
    fig.savefig(work_dir+"/study_level_front_position.png",
                bbox_inches='tight', dpi=dpi)


if __name__ == '__main__':

    work_dir = "Python/Cornea/Study_fg_vary_cpe4a6fe97-bd46-48a3-bc10-59a5ce1eda83"
#     work_dir = "Python/Cornea/Study_fg_vary_cp_randomd730ffa6-e920-4642-ad72-3d16baf302c7"
#     work_dir = "Python/Cornea/Study_pde_vary_h04ac8e66-da53-41ba-b161-2d70fc876a35"
    #work_dir = None

    if work_dir is None:
        work_dir = get_most_recently_modified_dir("Python/Cornea/")
    print "Working from: ", work_dir

    rc = ResultsCollection(work_dir)
    rc.load_results()

    do_density_line_plots(rc)
    #do_tip_density_x0(rc)
    #do_density_domain_comparisons(rc)
    #do_density_position_summary(rc)

    #do_line_density_overview(rc)
    #do_front_position_overview(rc)
