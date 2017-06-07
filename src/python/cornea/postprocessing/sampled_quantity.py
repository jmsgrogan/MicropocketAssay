"""
Module for postprocessing csv files output from simulations
"""

import csv
import numpy as np

def process_csv(file_name):
    
    results = []
    locations = []
    with open(file_name, 'rb') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        for idx, eachrow in enumerate(csv_reader):
            if idx==0:
                locations = eachrow[1:-1]
                locations = [float(x) for x in locations]
            else:
                time = float(eachrow[0])
                samples = eachrow[1:-1]
                samples = [float(x) for x in samples]
                results.append([time, samples])
                
    return np.array(locations), results