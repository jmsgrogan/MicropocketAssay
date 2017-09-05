import os
import csv
import glob
import numpy as np
from scipy.signal import convolve

import chaste.core
from microvessel_chaste.utility import *


def smooth_results(results):

    smooth_result = convolve(np.array(results), np.array([1, 2, 1]), 'same')/4.0
    smooth_result[0] = results[0]
    smooth_result[-1] = results[-1]
    return smooth_result


def process_csv(file_name):

    results = []
    locations = []
    with open(file_name, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        for idx, eachrow in enumerate(csv_reader):
            if idx == 0:
                locations = eachrow[1:-1]
                locations = [np.abs(float(x)) for x in locations]
            else:
                time = float(eachrow[0])
                samples = eachrow[1:-1]
                samples = [float(x) for x in samples]
                results.append([time, samples])
    return np.array(locations), results


def get_most_recently_modified_dir(search_dir):

    file_handler = chaste.core.OutputFileHandler(search_dir, False)
    work_dir = file_handler.GetOutputDirectoryFullPath()
    recent_dir = max(glob.glob(os.path.join(work_dir, '*/')), key=os.path.getmtime)
    return os.path.relpath(recent_dir, file_handler.GetChasteTestOutputDirectory())


def get_path(work_dir, study=None, domain=None, run=None, params=None):

    this_string = work_dir
    if study is not None:
        this_string += "/" + study + "/"
    if domain is not None:
        this_string += "/DomainType_" + domain.replace(" ", "") + "/"
    if run is not None:
        this_string += "/Run_" + run + "/"
    if params is not None:
        this_string += "/Sampled_" + params
    return this_string


class OutputParameter():

    def __init__(self, name, limits=None, line_color="C0", sampling_frequency=3):
        self.name = name
        self.limits = limits
        if limits is None:
            self.limits = [[0, 1200], None]
        self.line_color = line_color
        self.sampling_frequency = sampling_frequency

 