#!/usr/bin/env python

import sys
import os
import random
import time
from mpi4py import MPI
from microvessel_chaste.utility import *

import parameters.parameter_collection
from parameters.parameter_collection import SimulationParameterCollection, Parameter  

def launch(work_dir):
    
    random.seed(1234)
    
    collection = SimulationParameterCollection()
    collection.add_parameter(Parameter("pellet height", 1.0e-3*metre()))
    collection.add_parameter(Parameter("cornea radius", 1.3e-3*metre()))
    collection.add_parameter(Parameter("cornea thickness", 80.0e-6*metre()))    
    collection.add_parameter(Parameter("pellet thickness", 40.0e-6*metre()))
    collection.add_parameter(Parameter("grid spacing", 40.0e-6*metre()))  
    collection.add_parameter(Parameter("node spacing", 40.0e-6*metre()))  
    collection.add_parameter(Parameter("limbal offset", 200.0e-6*metre()))  
    collection.add_parameter(Parameter("density grid spacing", 40.0e-6*metre())) 
    collection.add_parameter(Parameter("sample spacing x", 20.0e-6*metre()))  
    collection.add_parameter(Parameter("sample spacing y", 20.0e-6*metre()))  
    collection.add_parameter(Parameter("sample spacing z", 20.0e-6*metre()))  
    collection.add_parameter(Parameter("use pellet", False))  
    collection.add_parameter(Parameter("pellet radius", 300.0e-6*metre()))  
    collection.add_parameter(Parameter("use finite pellet width", False)) 
    collection.add_parameter(Parameter("sprouting probability", 0.5 /(3600.0*second()), 
                                       min_val = 0.1, max_val = 10.0) )   
    
    # Set up the study
    study = parameters.parameter_collection.Study(work_dir, collection)
    study.range = [["sprouting probability", 5]]
    study.random_realisations = 3
    task_list = study.get_task_list()

    # Make one MPI subprocess per processor
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    rank = comm.rank
    size = comm.size
    
    new_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
    new_rank = new_comm.Get_rank()
    
    # Work through the tasks
    num_tasks = len(task_list)
    task_indices = []
    if num_tasks < size:
        if rank < num_tasks-1:
            task_indices.append(rank)
    else:
        tasks_per_proc = int(num_tasks/size)
        if rank<size - 1:
            task_indices = range(rank*tasks_per_proc, (rank+1)*tasks_per_proc)
        else:
            task_indices = range(rank*tasks_per_proc, num_tasks)            
            
    print task_indices
    for eachTaskIndex in task_indices:
        eachTask = task_list[eachTaskIndex]
        
        # Make the workdir
        if not os.path.exists(eachTask[0]):
            os.makedirs(eachTask[0])
            
        # dump the parameter collection
        eachTask[1].save(eachTask[0] + "/input_parameters.p")
        time.sleep(1)
        
        # Launch the task
        print "Launching task " + str(eachTaskIndex) + " on rank:", rank
        new_comm.Spawn("python", args=["simulation_runner.py", "-i "+ eachTask[0]], maxprocs=1)
        print "Completed task " + str(eachTaskIndex) + " on rank:", rank
        
    # Shutdown
    #comm.Disconnect()

if __name__ == "__main__":
    work_dir = "/home/grogan/test/"
    launch(work_dir)