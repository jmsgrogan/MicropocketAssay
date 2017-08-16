import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from microvessel_chaste.utility import *

from cornea.parameters.parameter_collection import SimulationParameterCollection
import cornea.analytical_solutions.solution_collection
from cornea.postprocessing import plotting_tools


_density_plot_keys = {"Line_density": r"Line Density - $\mu m$ per $\mu m^2$",
                      "Tip_density": r"Tip Density - $\mu m^{-2}$",
                      "Branch_density": r"Branch Density - $\mu m^{-1}$",
                      "PDE": r"Concentration - nanomolar"}

_location_plot_keys = {"location_min": "Front Position (um)",
                       "location_mid": "Mid Position (um)",
                       "location_max": "Max Position (um)",
                       "density_max": "Max Density (um^-2)"}


class PostProcessingTask(object):

    def __init__(self, work_dir, colormap=plt.cm.viridis, resolution=90):

        self.work_dir = work_dir
        self.colormap = colormap
        self.resolution = resolution

    def generate(self):

        self.load_data()

        # Make plot

        self.write()

    def load_data(self):

        pass

    def write(self):

        pass


class MergePlots(PostProcessingTask):

    def __init__(self, work_dir, filename, merge_files, merge_axis=0):

        super(MergePlots, self).__init__(work_dir)

        self.filename = filename
        self.merge_files = merge_files
        self.merge_axis = merge_axis

    def load_data(self):

        checked_files = []
        for eachFile in self.merge_files:
            if os.path.isfile(eachFile):
                checked_files.append(eachFile)
        self.merge_files = checked_files

    def write(self):

        images = map(Image.open, self.merge_files)
        print "Merging plots to: ", self.filename
        if len(images) > 0:
            widths, heights = zip(*(i.size for i in images))
            if self.merge_axis == 0:
                total_width = sum(widths)
                max_height = max(heights)
                new_im = Image.new('RGB', (total_width, max_height))
                x_offset = 0
                for im in images:
                    new_im.paste(im, (x_offset, 0))
                    x_offset += im.size[0]
                new_im.save(self.filename)
            else:
                max_width = max(widths)
                total_height = sum(heights)
                new_im = Image.new('RGB', (max_width, total_height))
                y_offset = 0
                for im in images:
                    new_im.paste(im, (0, y_offset))
                    y_offset += im.size[1]
                new_im.save(self.filename)


class DensityLinePlot(PostProcessingTask):

    def __init__(self, work_dir, study, domain, run_number, param, filename):

        super(DensityLinePlot, self).__init__(work_dir)

        self.left_line_loc = 0.0
        self.left_line_color = 'C0'
        self.right_line_loc = 1000.0
        self.right_line_color = 'C3'
        self.result_left_offset = 0
        self.study = study
        self.domain = domain
        self.param = param
        self.run_number = run_number
        self.pc = None
        self.results = None
        self.fig = None
        self.x_title = r"Position - $\mu m$"
        self.analytical_solution = None
        self.filename = filename

    def correct_for_thickness(self, results):

        if "density" in self.param.name and self.domain_is_3d(self.domain):
            thickness = self.pc.get_parameter("CorneaThickness").value
            thickness = thickness.Convert(1.0e-6*metres)
            results *= thickness
        return results

    def domain_is_3d(self, domain_type):
        return ("3D" in domain_type or "Hemi" in domain_type)

    def load_data(self):

        simulation_dir = plotting_tools.get_path(self.work_dir,
                                                 self.study,
                                                 self.domain,
                                                 str(self.run_number))

        self.pc = SimulationParameterCollection()
        if not os.path.isfile(simulation_dir + "input_parameters.p"):
            return
        self.pc.load(simulation_dir + "input_parameters.p")

        results_dir = plotting_tools.get_path(self.work_dir,
                                              self.study,
                                              self.domain,
                                              str(self.run_number),
                                              self.param.name)
        locations, values = plotting_tools.process_csv(results_dir + ".txt")
        sampled_values = values[::self.param.sampling_frequency]
        self.results = []
        for eachResult in sampled_values:
            current_time = eachResult[0]
            density_values = np.array(eachResult[1])
            density_values = self.correct_for_thickness(density_values)
            offset_result = density_values[self.result_left_offset:]
            offset_locations = locations[self.result_left_offset:]
            self.results.append([current_time, offset_result, offset_locations])

    def get_line_properties(self):

        height = self.pc.get_parameter("PelletHeight").value
        height = height.Convert(1.0e-6*metres)
        offset = self.pc.get_parameter("LimbalOffset").value
        offset = offset.Convert(1.0e-6*metres)
        self.left_line_loc = offset
        self.right_line_loc = height + offset
        if "Hemisphere" in self.domain:
            radius = self.pc.get_parameter("CorneaRadius").value
            radius = radius.Convert(1.0e-6*metres)
            self.right_line_loc = radius*np.arcsin(self.right_line_loc/radius)

    def generate(self):

        self.load_data()
        if self.results is None:
            return

        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        if self.param.limits is not None:
            if self.param.limits[0] is not None:
                ax.set_xlim(self.param.limits[0])
            if self.param.limits[1] is not None:
                ax.set_ylim(self.param.limits[1])
        self.get_line_properties()

        self.fig.ax.axvline(self.left_line_loc,
                            color=self.left_line_color,
                            linestyle='--', lw=1)
        self.fig.ax.axvline(self.right_line_loc,
                            color=self.right_line_color,
                            linestyle='--', lw=1)

        ax.set_xlabel(self.x_title)
        ax.set_ylabel(_density_plot_keys[self.param.name])
        max_result = 0.0
        colorscale = np.linspace(0, 1, len(self.results))
        colors = [self.colormap(i) for i in colorscale]

        for jdx, eachTimeStep in enumerate(self.results):
            time, results, locations = eachTimeStep
            smooth_result = plotting_tools.smooth_results(np.array(results))
            if self.analytical_solution is not None:
                analytical_method = getattr(cornea.analytical_solutions.solution_collection,
                                            self.analytical_solution)
                ana_results = analytical_method(locations, time, self.pc)
                self.fig.ax.plot(locations, ana_results, color='black', lw=1)
            self.fig.ax.plot(locations, smooth_result, color=colors[jdx], lw=1)
            local_max = np.max(smooth_result)
            if local_max > max_result:
                max_result = local_max
        ylim = ax.get_ylim()
        #self.fig.ax.set_ylim([0, ylim[1]])
        if "Tip" in self.param.name:
            self.fig.ax.set_ylim([0, 3.e-5])
        if "Line" in self.param.name:
            self.fig.ax.set_ylim([0, 0.018])
        self.fig.ax.annotate('Limbus', xy=(150, max_result*0.99), color="C0")
        self.fig.ax.annotate('Pellet', xy=(950, max_result*0.99), color="C3")
        self.write()

    def write(self):

        print "Saving density line plot to: ", self.filename
        self.fig.savefig(self.filename, bbox_inches='tight', dpi=self.resolution)


class DomainDensityComparison(PostProcessingTask):

    def __init__(self, work_dir, study, domains, num_random, param,
                 filename, location_key):

        super(DomainDensityComparison, self).__init__(work_dir)

        self.right_line_loc = 1000.0
        self.result_left_offset = 0
        self.study = study
        self.domains = domains
        self.param = param
        self.location_key = location_key
        self.num_random = num_random
        self.pc = None
        self.results = None
        self.fig = None
        self.x_title = r"Position - $\mu m$"
        self.analytical_solution = None
        self.filename = filename

    def get_line_properties(self):

        height = self.pc.get_parameter("PelletHeight").value
        height = height.Convert(1.0e-6*metres)
        offset = self.pc.get_parameter("LimbalOffset").value
        offset = offset.Convert(1.0e-6*metres)
        self.right_line_loc = height + offset
        if "Hemisphere" in self.domain:
            radius = self.pc.get_parameter("CorneaRadius").value
            radius = radius.Convert(1.0e-6*metres)
            self.right_line_loc = radius*np.arcsin(self.right_line_loc/radius)

    def load_data(self):

        simulation_dir = plotting_tools.get_path(self.study,
                                                 self.domain,
                                                 str(self.run_number))
        self.pc = SimulationParameterCollection()
        if not os.path.isfile(simulation_dir + "input_parameters.p"):
            return
        self.pc.load(simulation_dir + "input_parameters.p")

        results_dir = plotting_tools.get_path(self.study, self.domain,
                                              str(self.run_number),
                                              self.param.name)
        locations, values = plotting_tools.process_csv(results_dir + ".txt")
        sampled_values = values[::self.param.sampling_frequency]
        self.results = []
        for eachResult in sampled_values:
            current_time = eachResult[0]
            density_values = np.array(eachResult[1])
            density_values = self.correct_for_thickness(density_values)
            offset_result = density_values[self.result_left_offset:]
            offset_locations = locations[self.result_left_offset:]
            self.results.append([current_time, offset_result, offset_locations])

    def generate(self):

        self.load_data()
        if self.results is None:
            return

        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        self.fig.ax.set_xlabel("Time (hr)")
        self.fig.ax.set_ylabel(_location_plot_keys[self.location_key])
        if "location" in self.location_key:
            ax.set_ylim([0, self.right_line_loc])
        for eachDomain in self.domains:

            res = rc.results[eachStudy][eachDomain][param]["output"]
            for idx in range(self.num_random):
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
                    ax.plot(np.array(times), smooth_result, lw=1, label=eachDomain)
        ax.legend(loc='upper center')
        self.write()

    def write(self):

        self.fig.savefig(self.filename, bbox_inches='tight', dpi=self.resolution)