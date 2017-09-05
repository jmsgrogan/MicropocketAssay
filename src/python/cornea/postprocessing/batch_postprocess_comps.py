import pickle
import matplotlib
import matplotlib.pyplot as plt
import chaste
from microvessel_chaste.utility import *
import cornea.analytical_solutions
from cornea.postprocessing.plotting_tools import *
from cornea.postprocessing import plot_collection

# Matplotlib global settings
matplotlib.rcParams.update({'font.size': 18})
plt.locator_params(nticks=4)


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

    def setup_box_plots(self):

        for eachStudy in self.study_data["study_names"]:
            local_work_dir = get_path(self.work_dir, eachStudy)
            fig_dir = local_work_dir + "/box_plots/"
            if not os.path.exists(fig_dir):
                os.makedirs(fig_dir)
            task = plot_collection.BoxPlot(self.work_dir,
                                           eachStudy,
                                           self.study_data["domain_types"],
                                           len(self.study_data["random_seeds"]),
                                           self.parameters,
                                           fig_dir)
            self.tasks.append(task)

