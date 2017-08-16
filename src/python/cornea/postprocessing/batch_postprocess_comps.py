import os
import matplotlib.pyplot as plt
import numpy as np
import chaste
from microvessel_chaste.utility import *
import cornea.analytical_solutions
from cornea.postprocessing.plotting_tools import *
from cornea.postprocessing import plot_collection


class PostProcessingTaskManager(object):

    def __init__(self, relative_work_dir):

        self.tasks = []
        self.study_data = None
        self.parameters = [OutputParameter(name="Line_density"),
                           OutputParameter(name="Tip_density"),
                           OutputParameter(name="Branch_density"),
                           OutputParameter(name="PDE"), ]
        self.break_indices = []
        file_handler = chaste.core.OutputFileHandler(relative_work_dir, False)
        self.work_dir = file_handler.GetOutputDirectoryFullPath()
        self.load_study_data()

    def load_study_data(self):

        f = open(self.work_dir + "study_data.p", 'r')
        self.study_data = pickle.load(f)

    def setup_density_line_plots(self):

        for eachStudy in self.study_data["study_names"]:
            for eachDomain in self.study_data["domain_types"]:
                num_random = len(self.study_data["random_seeds"])
                for idx in range(num_random):
                    local_work_dir = get_path(self.work_dir, eachStudy, eachDomain, str(idx))
                    for eachParam in self.parameters:
                        fig_dir = local_work_dir + "/density_plots/"
                        if not os.path.exists(fig_dir):
                            os.makedirs(fig_dir)
                        fig_path = fig_dir + "/" + eachParam.name + ".png"
                        task = plot_collection.DensityLinePlot(self.work_dir,
                                                               eachStudy,
                                                               eachDomain,
                                                               idx,
                                                               eachParam,
                                                               fig_path)
#                         if "Tip" in eachParam.name:
#                             task.analytical_solution = 'get_tip_density_high_velocity'
#                         if "Line" in eachParam.name:
#                             task.analytical_solution = 'get_line_density_high_velocity'
                        self.tasks.append(task)

    def setup_line_density_plot_merge(self):

        for eachStudy in self.study_data["study_names"]:
            for eachDomain in self.study_data["domain_types"]:
                num_random = len(self.study_data["random_seeds"])
                outer_figure_paths = []
                for idx in range(num_random):
                    local_work_dir = get_path(self.work_dir,
                                              eachStudy,
                                              eachDomain,
                                              str(idx))
                    figure_paths = []
                    for eachParam in self.parameters:
                        fig_path = local_work_dir + "/density_plots/" + eachParam.name + ".png"
                        figure_paths.append(fig_path)
                    merge_path = local_work_dir + "/density_plots/merge.png"
                    outer_figure_paths.append(merge_path)
                    task = plot_collection.MergePlots(self.work_dir,
                                                      merge_path,
                                                      figure_paths,
                                                      merge_axis=0)
                    self.tasks.append(task)
                merge_path = get_path(self.work_dir, eachStudy, eachDomain) + "/merge.png"
                task = plot_collection.MergePlots(self.work_dir,
                                                  merge_path,
                                                  outer_figure_paths,
                                                  merge_axis=1)
                self.tasks.append(task)


    def setup_density_domain_comparisons(self):

        for eachStudy in self.study_data["study_names"]:
            work_dir = rc.get_path(eachStudy)
            for eachParam in rc.parameters:
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
                    fig.savefig(plot_path, bbox_inches='tight', dpi=dpi)

    def setup_density_domain_comparisons_merge():

        for eachStudy in rc.study_data["study_names"]:
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

    work_dir = "Python/Cornea/Study_fg_vary_cpc0d02423-8f2d-42ab-880b-5607f1659282"

    if work_dir is None:
        work_dir = get_most_recently_modified_dir("Python/Cornea/")
    print "Working from: ", work_dir

    tm = PostProcessingTaskManager(work_dir)
    tm.setup_density_line_plots()
