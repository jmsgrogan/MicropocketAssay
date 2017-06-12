#!/usr/bin/env python

import sys
import os
import random
from random import shuffle
import time
from mpi4py import MPI

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
    
    # Set up comms
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    status = MPI.Status()
    
    new_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
    new_rank = new_comm.Get_rank()
    
    # Work through the tasks
    task_indices = range(200)
    for eachTaskIndex in task_indices:
        time.sleep(1)

        print "Launching task " + str(eachTaskIndex) + " on rank:", rank, " of ", size
        newercomm = new_comm.Spawn(sys.executable, args=["simulation_runner.py", "-i test.txt"], maxprocs=1)
        data = newercomm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        print "Completed task " + str(eachTaskIndex) + " on rank:",  rank, " of ", size
        print "Source: ", source, " Tag ", tag
        
        # Shutdown
        newercomm.Disconnect()
    new_comm.Disconnect()

if __name__ == "__main__":
    work_dir = "Python/Cornea/ParamSweep_WithConsumption/"
    launch(work_dir)