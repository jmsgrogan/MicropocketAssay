#!/usr/bin/env python

import sys
import os
import random
from random import shuffle
import time
from mpi4py import MPI

import chaste.core
import cornea.parameters.parameter_collection
import cornea.parameters.default_parameters
import cornea.simulations

def get_task_indices(rank, size, task_list):
    
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
    return task_indices

def launch(work_dir):
    
    random.seed(1234)
     
    default_parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    default_parameter_collection.get_parameter("use pde only").value = True
    
    # Set up the study
    study = cornea.parameters.parameter_collection.Study(work_dir, default_parameter_collection)
    study.range = [["sprouting probability", 1]]
    study.random_realisations = 1
    task_list = study.get_task_list()

    # Set up comms
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    status = MPI.Status()
    
    new_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
    new_rank = new_comm.Get_rank()
    
    # Work through the tasks
    shuffle(task_list)
    task_indices = get_task_indices(rank, size, task_list)
    for eachTaskIndex in task_indices:
        eachTask = task_list[eachTaskIndex]
        
        # Make the workdir
        work_dir = eachTask[0]
        print work_dir
        param_collection = eachTask[1]
        file_handler = chaste.core.OutputFileHandler(work_dir, False)
        param_collection.save(file_handler.GetOutputDirectoryFullPath() + "/input_parameters.p")
        time.sleep(1)
        
        # Launch the task
        module_location = cornea.simulations.__file__
        print "Launching task " + str(eachTaskIndex) + " on rank:", rank, " of ", size
        newercomm = new_comm.Spawn(sys.executable, args=[os.path.dirname(module_location)+"/simulation_runner.py", "-i "+ file_handler.GetRelativePath()], maxprocs=1)
        data = newercomm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        print "Completed task " + str(eachTaskIndex) + " on rank:",  rank, " of ", size
        print "Source: ", source, " Tag ", tag
        
        # Shutdown
        newercomm.Disconnect()
    new_comm.Disconnect()

if __name__ == "__main__":
    work_dir = "Python/Cornea/ParamSweep/"
    launch(work_dir)