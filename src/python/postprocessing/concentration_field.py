"""
Module for postprocessing csv files output from simulations
"""

import csv
import numpy as np

def process_csv(file_name):
    
    concentrations = []
    locations = []
    with open(file_name, 'rb') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        current_day = -1
        for idx, eachrow in enumerate(csv_reader):
            if idx>0:
                file_day = int(round(float(eachrow[0])/24.0))
                if file_day>current_day:
                    print file_day
                    concentrations.append([float(i) for i in eachrow[1:-1]])
                    current_day = file_day
            else:
                locations = np.array([float(i) for i in eachrow[1:-1]])
    return locations, np.array(concentrations)