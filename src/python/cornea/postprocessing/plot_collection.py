import os
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt
import vtk
from PIL import Image
from microvessel_chaste.utility import *

from cornea.parameters.parameter_collection import SimulationParameterCollection
import cornea.analytical_solutions.solution_collection
from cornea.postprocessing import plotting_tools
from cornea.postprocessing import sampling_grid

# Global matplotlib settings
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['axes.linewidth'] = 2.0  # set the value globally
plt.locator_params(nticks=4)

_density_plot_keys = {"Line_density": r"Line Density - $\mu m$ per $\mu m^3$",
                      "Tip_density": r"Tip Density - $\mu m^{-3}$",
                      "PDE": r"Concentration - nanomolar"}

_location_plot_keys = {"location_min": "Front Position (um)",
                       "location_mid": "Mid Position (um)",
                       "location_max": "Max Position (um)",
                       "density_max": "Max Density (um^-3)"}


class PostProcessingTask(object):

    def __init__(self, work_dir, colormap=plt.cm.viridis, resolution=180):

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
                            color=self.left_line_color, lw=3.0)
        self.fig.ax.axvline(self.right_line_loc,
                            color=self.right_line_color, lw=3.0)

        #ax.set_xlabel(self.x_title)
        #ax.set_ylabel(_density_plot_keys[self.param.name])
        max_result = 0.0
        colorscale = np.linspace(0, 1, len(self.results))
        colors = [self.colormap(i) for i in colorscale]

        result_factor = 1.0
        if "Tip" in self.param.name:
            result_factor = 1.e-8
        elif "Line" in self.param.name:
            result_factor = 1.e-5

        for jdx, eachTimeStep in enumerate(self.results):
            time, results, locations = eachTimeStep
            smooth_result = plotting_tools.smooth_results(np.array(results))
            if "Tip" in self.param.name:
                print "Max density: ", self.domain, " t: ", time, " d: ", 1.e6*np.max(smooth_result)
            
            if self.analytical_solution is not None:
                analytical_method = getattr(cornea.analytical_solutions.solution_collection,
                                            self.analytical_solution)
                ana_results = analytical_method(locations, time, self.pc)
                self.fig.ax.plot(locations, ana_results/result_factor, color='black', lw=3.0)
            self.fig.ax.plot(locations, smooth_result/result_factor, color=colors[jdx], lw=3.0)
            local_max = np.max(smooth_result/result_factor)
            if local_max > max_result:
                max_result = local_max
        #ylim = ax.get_ylim()
        #self.fig.ax.set_ylim([0, ylim[1]])
        if "Tip" in self.param.name:
            y_max = 600
        elif "Line" in self.param.name:
            y_max = 120
        else:
            y_max = ax.get_ylim()[1]

        self.fig.ax.set_xticks([0, 250, 500, 750, 1000])
        self.fig.ax.set_ylim([0, y_max])
        self.write()

    def write(self):

        print "Saving density line plot to: ", self.filename
        self.fig.savefig(self.filename, bbox_inches='tight', dpi=self.resolution)


class BoxPlot(PostProcessingTask):

    def __init__(self, work_dir, study, domains, num_random, params, filename):

        super(BoxPlot, self).__init__(work_dir)

        self.study = study
        self.colormap = plt.cm.Accent
        self.domains = domains
        self.params = params
        self.num_random = num_random
        self.result_left_offset = 0
        self.results = None
        self.fig = None
        self.x_title = r""
        self.filename = filename
        self.work_dir = work_dir

    def load_data(self):

        self.results = {}

        params = ["Tip", "Line"]
        for eachDomain in self.domains:
            self.results[eachDomain] = []
            for jdx, param_name in enumerate(params):
                self.results[eachDomain].append({"density_max": [],
                                                 "location_max": [],
                                                 "location_mid": [],
                                                 "location_min": []})
                for idx in range(self.num_random):
                    simulation_dir = plotting_tools.get_path(self.work_dir,
                                                             self.study,
                                                             eachDomain,
                                                             str(idx))
                    pc = SimulationParameterCollection()
                    if not os.path.isfile(simulation_dir + "input_parameters.p"):
                        continue
                    pc.load(simulation_dir + "input_parameters.p")

                    locations, values = plotting_tools.process_csv(simulation_dir + "Sampled_" + param_name + "_density.txt")
                    sampled_values = values[::self.params[jdx].sampling_frequency]
                    last_time_result = sampled_values[-1]
                    time = last_time_result[0]
                    print "density sample time: ", time
                    profile = np.array(last_time_result[1])

                    offset_result = profile[self.result_left_offset:]
                    smooth_result = plotting_tools.smooth_results(np.array(offset_result))
                    offset_locations = locations[self.result_left_offset:]

                    self.results[eachDomain][jdx]["density_max"].append(np.max(smooth_result))
                    max_arg = np.argmax(smooth_result)
                    self.results[eachDomain][jdx]["location_max"].append(offset_locations[max_arg])
                    mid_val = smooth_result[max_arg]/2.0

                    diff = 1.e12
                    result_size = len(smooth_result)
                    mid_index = result_size-1
                    for kdx in range(result_size):
                        check_index = result_size-1-kdx
                        if abs(smooth_result[check_index] - mid_val) < diff:
                            diff = abs(smooth_result[check_index] - mid_val)
                            if check_index >= max_arg:
                                mid_index = check_index
                    self.results[eachDomain][jdx]["location_mid"].append(offset_locations[mid_index]) 

                    min_val = 0.01*smooth_result[max_arg]
                    diff = 1.e12
                    min_index = result_size-1
                    for kdx in range(result_size):
                        check_index = result_size-1-kdx
                        if abs(smooth_result[check_index] - min_val) <= diff:
                            diff = abs(smooth_result[check_index] - min_val)
                            if check_index >= mid_index:
                                min_index = check_index
                    self.results[eachDomain][jdx]["location_min"].append(offset_locations[min_index]) 

    def generate(self):

        self.load_data()
        if self.results is None:
            return

        colorscale = np.linspace(0, 1, len(self.domains))
        colors = [self.colormap(i) for i in colorscale]
        width = 1.0
        alpha = 0.1
        
        domain_abbreviations = {"Planar_2D": "P2D",
                       "Planar_2D_Finite": "P2DF",
                       "Circle_2D": "C2D",
                       "Planar_3D":  "P3D",
                       "Planar_3D_Finite": "P3DF",
                       "Circle_3D": "C3D",
                       "Hemisphere": "H"}

        # Positions
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_xlabel("Location (um)")
        ind = np.arange(len(self.domains))
        ax.set_yticks(ind)
        ax.set_yticklabels([domain_abbreviations[x] for x in self.domains])
        for idx, eachDomain in enumerate(self.domains):
            print "Loc ", eachDomain
            results = self.results[eachDomain][0]["location_max"]
            print "max ", np.mean(results)
            ax.barh(ind[idx], np.mean(results), width, color='navy', xerr=np.std(results), alpha=alpha, edgecolor='black')
            results = self.results[eachDomain][0]["location_mid"]
            print "mid ", np.mean(results)
            ax.barh(ind[idx], np.mean(results), width, color='navy', xerr=np.std(results), alpha=alpha, edgecolor='black')            
            results = self.results[eachDomain][0]["location_min"]
            print "min ", np.mean(results)
            ax.barh(ind[idx], np.mean(results), width, color='navy', xerr=np.std(results), alpha=alpha, edgecolor='black')
            ax.invert_yaxis()             
        ax.set_xlim([0, 1200])
        self.fig.savefig(self.filename+"locations.png", bbox_inches='tight', dpi=self.resolution)
        
        # Density
        width = 6.0
        alpha = 0.5
        
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_xlabel("Tip Density (1e6 um^-3)")
        ax.set_xlim([0, 1.5])  # Fig 5
        #ax.set_xlim([0, 20])  # Fig 4
        ind = 6.0*np.arange(len(self.domains))
        tick_marks = list(ind)
        tick_marks.extend(list(ind+48.0))
        
        #ax.set_xticks(tick_marks)
        tick_labels = [domain_abbreviations[x] for x in self.domains]
        tick_labels.extend([domain_abbreviations[x] for x in self.domains])
        ax.get_yaxis().set_visible(False)
        #ax.set_xticklabels(tick_labels)
        for idx, eachDomain in enumerate(self.domains):
            results = np.array(self.results[eachDomain][0]["density_max"])
            print "Density ", eachDomain, " mean: ",  1.e6*np.mean(results), " err: ",  1.e6*np.std(results)
            ax.barh(ind[idx] + 0.0*width, 1.e6*np.mean(results), width, color='navy', xerr=1.e6*np.std(results), alpha=alpha, edgecolor='black')
            ax.invert_yaxis()   
        self.fig.savefig(self.filename+"tip_density.png", bbox_inches='tight', dpi=self.resolution)        

        width = 6.0
        alpha = 0.5

        # Density
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_xlabel("Line Density (1e3 um^-2)")
        ax.set_xlim([0, 3.0])   # Fig 4
        ax.set_xlim([0, 0.3])   # Fig 5
        ind = 6.0*np.arange(len(self.domains))
        tick_marks = list(ind)
        tick_marks.extend(list(ind+48.0))
        
        #ax.set_xticks(tick_marks)
        tick_labels = [domain_abbreviations[x] for x in self.domains]
        tick_labels.extend([domain_abbreviations[x] for x in self.domains])
        ax.get_yaxis().set_visible(False)
        #ax.set_xticklabels(tick_labels)
        for idx, eachDomain in enumerate(self.domains):
            results = np.array(self.results[eachDomain][1]["density_max"])
            ax.barh(ind[idx] + 0.0*width, 1.e3*np.mean(results), width, color='lightsteelblue', xerr=1.e3*np.std(results), alpha=alpha, edgecolor='black')
            ax.invert_yaxis()   
        self.fig.savefig(self.filename+"line_density.png", bbox_inches='tight', dpi=self.resolution)

    def write(self):
        pass


class PdePlot(PostProcessingTask):

    def __init__(self, work_dir, study, domains, num_random, params, filename):

        super(PdePlot, self).__init__(work_dir)

        self.study = study
        self.colormap = plt.cm.Accent
        self.domains = domains
        self.params = params
        self.num_random = num_random
        self.result_left_offset = 0
        self.results = None
        self.fig = None
        self.x_title = r""
        self.filename = filename
        self.work_dir = work_dir

    def load_data(self):

        self.results = {}

        params = ["Tip", "Line"]
        for eachDomain in self.domains:
            self.results[eachDomain] = {"line_frac": [],
                                        "angle": []}
            for idx in range(self.num_random):
                simulation_dir = plotting_tools.get_path(self.work_dir,
                                                         self.study,
                                                         eachDomain,
                                                         str(idx))
                pc = SimulationParameterCollection()
                if not os.path.isfile(simulation_dir + "input_parameters.p"):
                    continue
                pc.load(simulation_dir + "input_parameters.p")
                last_file = max(glob.glob(simulation_dir + "/sampled_density*.vtu"))
                print "last", last_file
                line_fraction, angle = sampling_grid.GetDensityMetrics(last_file, eachDomain, pc)
                sampling_grid.DoLineSampling(last_file, eachDomain, pc)
                self.results[eachDomain]["line_frac"].append(line_fraction)
                self.results[eachDomain]["angle"].append(angle)

    def generate(self):

        self.load_data()
        if self.results is None:
            return

        colorscale = np.linspace(0, 1, len(self.domains))
        colors = [self.colormap(i) for i in colorscale]
        width = 1.0
        alpha = 0.6

        domain_abbreviations = {"Planar_2D": "P2D",
                       "Planar_2D_Finite": "P2DF",
                       "Circle_2D": "C2D",
                       "Planar_3D":  "P3D",
                       "Planar_3D_Finite": "P3DF",
                       "Circle_3D": "C3D",
                       "Hemisphere": "H"}

        # Density
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_ylabel("Vascularized Fraction")
        ind = np.arange(len(self.domains))
        ax.set_xticks(ind)
        ax.set_ylim([0, 0.9])
        ax.set_xticklabels([domain_abbreviations[x] for x in self.domains])
        for idx, eachDomain in enumerate(self.domains):
            results = np.array(self.results[eachDomain]["line_frac"])
            print "vasc frac ", eachDomain, " : ", np.mean(results)
            ax.bar(ind[idx] + 0.0*width, np.mean(results), width, color=colors[idx], yerr=np.std(results), alpha=alpha, edgecolor='black')                    
        self.fig.savefig(self.filename+"/vascularized_fraction.png", bbox_inches='tight', dpi=self.resolution)

        # Angle
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_ylabel("Opening Angle")
        ind = np.arange(len(self.domains))
        ax.set_xticks(ind + width)
        ax.set_xticklabels([domain_abbreviations[x] for x in self.domains])
        for idx, eachDomain in enumerate(self.domains):
            results = np.array(self.results[eachDomain]["angle"])
            ax.bar(ind[idx] + 0.0*width, np.mean(results), width, color=colors[idx], yerr=np.std(results), alpha=alpha, edgecolor='black')                    
        self.fig.savefig(self.filename+"/opening_angle.png", bbox_inches='tight', dpi=self.resolution)

    def write(self):
        pass


class MaxTipDensityPlot(PostProcessingTask):

    def __init__(self, work_dir, study, domains, num_random, params, filename):

        super(MaxTipDensityPlot, self).__init__(work_dir)

        self.study = study
        self.colormap = plt.cm.Accent
        self.domains = domains
        self.params = params
        self.num_random = num_random
        self.result_left_offset = 0
        self.results = None
        self.fig = None
        self.x_title = r""
        self.filename = filename
        self.work_dir = work_dir

    def load_data(self):

        self.results = {}

        for eachDomain in self.domains:
            self.results[eachDomain] = {"times": [],
                                        "max_vals": []}
            simulation_dir = plotting_tools.get_path(self.work_dir,
                                                     self.study,
                                                     eachDomain,
                                                     str(0))
            pc = SimulationParameterCollection()
            if not os.path.isfile(simulation_dir + "input_parameters.p"):
                continue
            pc.load(simulation_dir + "input_parameters.p")
            locations, values = plotting_tools.process_csv(simulation_dir + "Sampled_Tip_density.txt")
            sampled_values = values[::1]

            for eachTimeStep in sampled_values:
                self.results[eachDomain]["times"].append(eachTimeStep[0])
                profile = np.array(eachTimeStep[1])
                offset_result = profile[self.result_left_offset:]
                smooth_result = plotting_tools.smooth_results(np.array(offset_result))
                self.results[eachDomain]["max_vals"].append(np.max(smooth_result))

    def forceAspect(self, ax,aspect=1):
        #aspect is width/height
        scale_str = ax.get_yaxis().get_scale()
        xmin,xmax = ax.get_xlim()
        ymin,ymax = ax.get_ylim()
        if scale_str=='linear':
            asp = abs((xmax-xmin)/(ymax-ymin))/aspect
        elif scale_str=='log':
            asp = abs((scipy.log(xmax)-scipy.log(xmin))/(scipy.log(ymax)-scipy.log(ymin)))/aspect
        ax.set_aspect(asp)
    
    def generate(self):

        self.load_data()
        if self.results is None:
            return

        colorscale = np.linspace(0, 1, len(self.domains))
        colors = [self.colormap(i) for i in colorscale]

        # Density
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Max Tip Density x 10-6 um^3")
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 90])
        self.forceAspect(ax, 2)
        for idx, eachDomain in enumerate(self.domains):
            times = np.array(self.results[eachDomain]["times"])
            vals = 1e6*np.array(self.results[eachDomain]["max_vals"])
            ax.plot(times, vals, color=colors[idx], lw=3, label=eachDomain.replace("_",""))    
        handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles, labels) 
        #ax.legend(frameon=False)               
        self.fig.savefig(self.filename+"/max_tip_density.png", bbox_inches='tight', dpi=self.resolution)
    def write(self):
        pass
    

class MaxConcPlot(PostProcessingTask):

    def __init__(self, work_dir, study, domains, num_random, params, filename):

        super(MaxConcPlot, self).__init__(work_dir)

        self.study = study
        self.colormap = plt.cm.Accent
        self.domains = domains
        self.params = params
        self.num_random = num_random
        self.result_left_offset = 0
        self.results = None
        self.fig = None
        self.x_title = r""
        self.filename = filename
        self.work_dir = work_dir

    def load_data(self):

        self.results = {}

        for eachDomain in self.domains:
            self.results[eachDomain] = {"times": [],
                                        "max_vals": []}
            simulation_dir = plotting_tools.get_path(self.work_dir,
                                                     self.study,
                                                     eachDomain,
                                                     str(0))
            pc = SimulationParameterCollection()
            if not os.path.isfile(simulation_dir + "input_parameters.p"):
                continue
            pc.load(simulation_dir + "input_parameters.p")
            extension = "vtu"
            if ("Planar" in eachDomain) and (not "Finite" in eachDomain):
                extension = "vti"

            vegf_files = glob.glob(simulation_dir + "/Vegf_Solution*."+extension)
            vegf_files.sort(key=lambda f: int(filter(str.isdigit, os.path.basename(str(f)))))
            for eachFile in vegf_files:
                time = filter(lambda x: x.isdigit(), os.path.basename(eachFile))
                print os.path.basename(eachFile), time
                self.results[eachDomain]["times"].append(time)

                if extension=="vtu":
                    reader = vtk.vtkXMLUnstructuredGridReader()
                    reader.SetFileName(eachFile)
                    reader.Update()
                else:
                    reader = vtk.vtkXMLImageDataReader()
                    reader.SetFileName(eachFile)
                    reader.Update()
                data_set = reader.GetOutput()
                solution = data_set.GetPointData().GetArray("vegf")
                max_val = 0.0
                for idx in range(solution.GetNumberOfTuples()):
                    point_soln = solution.GetTuple1(idx)
                    if point_soln > max_val:
                        max_val = point_soln
                self.results[eachDomain]["max_vals"].append(max_val)

    def forceAspect(self, ax,aspect=1):
        #aspect is width/height
        scale_str = ax.get_yaxis().get_scale()
        xmin,xmax = ax.get_xlim()
        ymin,ymax = ax.get_ylim()
        if scale_str=='linear':
            asp = abs((xmax-xmin)/(ymax-ymin))/aspect
        elif scale_str=='log':
            asp = abs((scipy.log(xmax)-scipy.log(xmin))/(scipy.log(ymax)-scipy.log(ymin)))/aspect
        ax.set_aspect(asp)
    
    def generate(self):

        self.load_data()
        if self.results is None:
            return

        colorscale = np.linspace(0, 1, len(self.domains))
        colors = [self.colormap(i) for i in colorscale]

        # Density
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Concentration (nM)")
        ax.set_ylim([0, 20])
        ax.set_xlim([0, 90])
        self.forceAspect(ax, 2)
        for idx, eachDomain in enumerate(self.domains):
            times = np.array(self.results[eachDomain]["times"])
            vals = 1e-6*np.array(self.results[eachDomain]["max_vals"])
            ax.plot(times, vals, color=colors[idx], lw=3, label=eachDomain.replace("_",""))    
        handles, labels = ax.get_legend_handles_labels()
#         ax.legend(handles, labels) 
#         ax.legend(frameon=False)               
        self.fig.savefig(self.filename+"/concentrations.png", bbox_inches='tight', dpi=self.resolution)
    def write(self):
        pass
    
    
class FrontPosPlot(PostProcessingTask):

    def __init__(self, work_dir, study, domains, num_random, params, filename):

        super(FrontPosPlot, self).__init__(work_dir)

        self.study = study
        self.colormap = plt.cm.Accent
        self.domains = domains
        self.params = params
        self.num_random = num_random
        self.result_left_offset = 0
        self.results = None
        self.fig = None
        self.x_title = r""
        self.filename = filename
        self.work_dir = work_dir

    def load_data(self):

        self.results = {}

        for eachDomain in self.domains:
            self.results[eachDomain] = {"times": [],
                                        "max_vals": []}
            simulation_dir = plotting_tools.get_path(self.work_dir,
                                                     self.study,
                                                     eachDomain,
                                                     str(0))
            pc = SimulationParameterCollection()
            if not os.path.isfile(simulation_dir + "input_parameters.p"):
                continue
            pc.load(simulation_dir + "input_parameters.p")
            locations, values = plotting_tools.process_csv(simulation_dir + "Sampled_Tip_density.txt")
            sampled_values = values[::1]

            for eachTimeStep in sampled_values:
                self.results[eachDomain]["times"].append(eachTimeStep[0])
                profile = np.array(eachTimeStep[1])
                offset_result = profile[self.result_left_offset:]
                smooth_result = plotting_tools.smooth_results(np.array(offset_result))
                offset_locations = locations[self.result_left_offset:]

                max_arg = np.argmax(smooth_result)
                mid_val = smooth_result[max_arg]/2.0

                diff = 1.e12
                result_size = len(smooth_result)
                mid_index = result_size-1
                for kdx in range(result_size):
                    check_index = result_size-1-kdx
                    if abs(smooth_result[check_index] - mid_val) < diff:
                        diff = abs(smooth_result[check_index] - mid_val)
                        if check_index >= max_arg:
                            mid_index = check_index

                min_val = 0.01*smooth_result[max_arg]
                diff = 1.e12
                min_index = result_size-1
                for kdx in range(result_size):
                    check_index = result_size-1-kdx
                    if abs(smooth_result[check_index] - min_val) <= diff:
                        diff = abs(smooth_result[check_index] - min_val)
                        if check_index >= mid_index:
                            min_index = check_index
                self.results[eachDomain]["max_vals"].append(offset_locations[min_index])                 

    def forceAspect(self, ax,aspect=1):
        #aspect is width/height
        scale_str = ax.get_yaxis().get_scale()
        xmin,xmax = ax.get_xlim()
        ymin,ymax = ax.get_ylim()
        if scale_str=='linear':
            asp = abs((xmax-xmin)/(ymax-ymin))/aspect
        elif scale_str=='log':
            asp = abs((scipy.log(xmax)-scipy.log(xmin))/(scipy.log(ymax)-scipy.log(ymin)))/aspect
        ax.set_aspect(asp)
    
    def generate(self):

        self.load_data()
        if self.results is None:
            return

        colorscale = np.linspace(0, 1, len(self.domains))
        colors = [self.colormap(i) for i in colorscale]

        # Density
        self.fig, ax = plt.subplots()
        self.fig.ax = ax
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Location (um)")
#         ax.set_ylim([0, 1])
#         ax.set_xlim([0, 90])
        #self.forceAspect(ax, 2)
        for idx, eachDomain in enumerate(self.domains):
            times = np.array(self.results[eachDomain]["times"])
            vals = np.array(self.results[eachDomain]["max_vals"])
            ax.plot(times, vals, color=colors[idx], lw=3, label=eachDomain.replace("_",""))    
        handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles, labels) 
        #ax.legend(frameon=False)               
        self.fig.savefig(self.filename+"/front_location.png", bbox_inches='tight', dpi=self.resolution)
    def write(self):
        pass
        