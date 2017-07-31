import os
import glob
import pickle
import copy
import numpy as np
from scipy.signal import convolve
from PIL import Image

import chaste.core
from microvessel_chaste.utility import *
from cornea.postprocessing.sampled_quantity import process_csv
from cornea.parameters.parameter_collection import SimulationParameterCollection


def smooth_results(results):

    smooth_result = convolve(np.array(results), np.array([1, 2, 1]), 'same')/4.0
    smooth_result[0] = results[0]
    smooth_result[-1] = results[-1]
    return smooth_result


def get_most_recently_modified_dir(search_dir):

    file_handler = chaste.core.OutputFileHandler(search_dir, False)
    work_dir = file_handler.GetOutputDirectoryFullPath()
    recent_dir = max(glob.glob(os.path.join(work_dir, '*/')), key=os.path.getmtime)
    return os.path.relpath(recent_dir, file_handler.GetChasteTestOutputDirectory())


class OutputParameter():

    def __init__(self, name, limits=None, line_color="C0", sampling_frequency=1):
        self.name = name
        self.limits = limits
        if limits is None:
            self.limits = [[0, 1200], None]
        self.line_color = line_color
        self.sampling_frequency = sampling_frequency


class ResultsCollection(object):

    def __init__(self, relative_work_dir):

        file_handler = chaste.core.OutputFileHandler(relative_work_dir, False)
        self.work_dir = file_handler.GetOutputDirectoryFullPath()
        self.study_data = None
        self.results = {}
        self.base_result_template = {"location_min": 0.0,
                                     "location_mid": 0.0,
                                     "location_max": 0.0,
                                     "density_max": 0.0,
                                     "locations": [],
                                     "results": [],
                                     "time": 0.0}
        self.location_plot_keys = {"location_min": "Front Position (um)",
                                   "location_mid": "Mid Position (um)",
                                   "location_max": "Max Position (um)",
                                   "density_max": "Max Density (um^-2)"}
        self.density_plot_keys = {"Line_density": r"Line Density - $\mu m$ per $\mu m^2$",
                                  "Tip_density": r"Line Density - $\mu m$ per $\mu m^2$",
                                  "Branch_density": r"Branch Density - $\mu m^{-1}$",
                                  "PDE": r"Concentration - nanomolar"}
        self.parameters = [OutputParameter(name="Line_density"),
                           OutputParameter(name="Tip_density"),
                           OutputParameter(name="Branch_density"),
                           OutputParameter(name="PDE"), ]
        self.result_left_offset = 0
        self.load_study_data()

    def load_study_data(self):

        f = open(self.work_dir + "study_data.p", 'r')
        self.study_data = pickle.load(f)

    def setup_results(self):

        for eachStudy in self.study_data["study_names"]:
            self.results[eachStudy] = {}
            for eachDomain in self.study_data["domain_types"]:
                self.results[eachStudy][eachDomain] = {}
                for eachParam in self.parameters:
                    self.results[eachStudy][eachDomain][eachParam.name] = {"params": [],
                                                                           "output": [],
                                                                           "summaries": {}}

    def domain_is_3d(self, domain_type):
        return ("3D" in domain_type or "Hemi" in domain_type)

    def get_path(self, study=None, domain=None, run=None, params=None):

        string = self.work_dir
        if study is not None:
            string += "/" + study + "/"
        if domain is not None:
            string += "/DomainType_" + domain.replace(" ", "") + "/"
        if run is not None:
            string += "/Run_" + run + "/"
        if params is not None:
            string += "/Sampled_" + params
        return string

    def get_nearest_index(self, value, result, locations):
        return locations[(np.abs(result-value)).argmin()]

    def load_results(self):
        self.setup_results()

        for eachStudy in self.study_data["study_names"]:
            for eachDomain in self.study_data["domain_types"]:
                num_random = len(self.study_data["random_seeds"])
                for idx in range(num_random):
                    pc_dir = self.get_path(eachStudy, eachDomain, str(idx))
                    pc = SimulationParameterCollection()
                    pc.load(pc_dir + "input_parameters.p")
                    for eachParam in self.parameters:
                        self.results[eachStudy][eachDomain][eachParam.name]["params"].append(pc)
                        results_dir = self.get_path(eachStudy, eachDomain, str(idx), eachParam.name)
                        locations, values = process_csv(results_dir + ".txt")
                        thickness = pc.get_parameter("CorneaThickness").value
                        thickness = thickness.Convert(1.e-6*metres)
                        sampled_values = values[::eachParam.sampling_frequency]
                        local_time_series = []
                        for eachResult in sampled_values:
                            local_results = copy.deepcopy(self.base_result_template)
                            local_results["time"] = eachResult[0]
                            result = np.array(eachResult[1])
                            if "density" in eachParam.name and self.domain_is_3d(eachDomain):
                                result *= thickness
                            offset_result = result[self.result_left_offset:]
                            offset_locations = locations[self.result_left_offset:]
                            local_results["results"] = offset_result
                            local_results["locations"] = offset_locations
                            local_results["density_max"] = np.max(offset_result)
                            max_arg = np.argmax(offset_result)
                            local_results["location_max"] = offset_locations[max_arg]
                            mid_val = offset_result[max_arg]/2.0
                            mid_index = (np.abs(offset_result-mid_val)).argmin()
                            if mid_index < max_arg:
                                mid_index = max_arg
                            local_results["location_mid"] = offset_locations[mid_index]
                            min_val = 0.01*offset_result[max_arg]
                            min_index = (np.abs(offset_result-min_val)).argmin()
                            if min_index < mid_index:
                                min_index = mid_index
                            local_results["location_min"] = offset_locations[min_index]
                            local_time_series.append(local_results)
                        self.results[eachStudy][eachDomain][eachParam.name]["output"].append(local_time_series)

        # Calculate summaries
        for eachStudy in self.study_data["study_names"]:
            for eachParam in self.parameters:
                param = eachParam.name
                output_vars = self.location_plot_keys.keys()
                for eachVar in output_vars:
                    for eachDomain in self.study_data["domain_types"]:
                        self.results[eachStudy][eachDomain][param]["summaries"][eachVar] = {"midpoints": [],
                                                                                            "slopes": []}
                        num_random = len(self.study_data["random_seeds"])
                        for idx in range(num_random):
                            time_series = self.results[eachStudy][eachDomain][param]["output"][idx]
                            times = []
                            local_results = []
                            for eachTime in time_series:
                                times.append(eachTime["time"])
                                local_results.append(eachTime[eachVar])
                            # Get midpoints and slopes
                            dt = times[1] - times[0]
                            smooth_result = smooth_results(local_results)
                            mid_index = int(round(len(smooth_result)/2.0))
                            md_val = smooth_result[mid_index]
                            lower_index = int(round(len(smooth_result)/4.0))
                            upper_index = int(round(3.0*len(smooth_result)/4.0))
                            lower_val = smooth_result[lower_index]
                            lower_time = float(lower_index)*dt
                            upper_val = smooth_result[upper_index]
                            upper_time = float(upper_index)*dt
                            slope = (upper_val-lower_val)/(upper_time-lower_time)
                            self.results[eachStudy][eachDomain][param]["summaries"][eachVar]["slopes"].append(slope)
                            self.results[eachStudy][eachDomain][param]["summaries"][eachVar]["midpoints"].append(md_val)   


def merge_images_x(output_path, input_paths):

    images = map(Image.open, input_paths)
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]
    new_im.save(output_path)


def merge_images_y(output_path, input_paths):

    images = map(Image.open, input_paths)
    widths, heights = zip(*(i.size for i in images))
    max_width = max(widths)
    total_height = sum(heights)

    new_im = Image.new('RGB', (max_width, total_height))
    y_offset = 0
    for im in images:
        new_im.paste(im, (0, y_offset))
        y_offset += im.size[1]
    new_im.save(output_path)
