#!/usr/bin/env python
"""Demonstrate the task-pull paradigm for high-throughput computing
using mpi4py. Task pull is an efficient way to perform a large number of
independent tasks when there are more tasks than processors, especially
when the run times vary for each task. 
This code is over-commented for instructional purposes.
This example was contributed by Craig Finch (cfinch@ieee.org).
Inspired by http://math.acadiau.ca/ACMMaC/Rmpi/index.html
"""

import os
import time
import pickle
import chaste
import microvessel_chaste.simulation
from microvessel_chaste.utility import *
from cornea.simulations.full_activation import study, master_work_dir, study_data
from mpi4py import MPI


def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START')

# Initializations and preliminaries
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object
local_comm = comm.Split(color=rank, key=rank)
chaste.init(comm=local_comm)

if rank == 0:

    # Master process executes code below
    test_output = os.getenv('CHASTE_TEST_OUTPUT', os.getcwd())
    work_dir = test_output + "/" + master_work_dir + "/"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    pickle.dump(study_data, open(work_dir + "study_data.p", "wb"))

    tasks = study
    task_index = 0
    num_workers = size - 1
    closed_workers = 0

    while closed_workers < num_workers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # Worker is ready, so send it a task
            if task_index < len(tasks):
                comm.send(tasks[task_index], dest=source, tag=tags.START)
                task_index += 1
            else:
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            results = data
        elif tag == tags.EXIT:
            closed_workers += 1
    print "Master finishing"

else:

    # Give master time to make the work directory
    time.sleep(1)

    # Worker processes execute code below
    name = MPI.Get_processor_name()
    print "worker start"

    while True:
        comm.send(None, dest=0, tag=tags.READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        if tag == tags.START:

            # Do the work here
            print "Parent rank", comm.Get_rank()
            print "Worker rank", local_comm.Get_rank()
            print "Chaste rank", chaste.chaste.core.PetscTools.GetMyRank()

            study, run_number, domain, pc, work_dir = task

            file_handler = chaste.core.OutputFileHandler(work_dir, False)
            int_work_dir = work_dir + "/" + study["name"]
            domain_string = "/DomainType_" + domain.replace(" ", "")
            run_string = "/Run_" + str(run_number) + "/"
            local_ext = domain_string + run_string
            local_work_dir = int_work_dir + local_ext
            if "2" in domain:
                simulation = microvessel_chaste.simulation.CornealMicropocketSimulation2()
            else:
                simulation = microvessel_chaste.simulation.CornealMicropocketSimulation3()
            simulation.SetWorkDir(local_work_dir)
            for eachParameter in pc.collection.keys():
                if eachParameter not in ["DomainType", "RandomSeed", "RunNumber"]:
                    param_nam = pc.collection[eachParameter].name
                    param_value = pc.collection[eachParameter].value
                    getattr(simulation, 'Set'+param_nam)(param_value)
            domain_type = microvessel_chaste.simulation.DomainType
            domain_key = domain
            if "Finite" in domain:
                domain_key = domain[0:-7]
            simulation.SetRandomSeed(int(pc.get_parameter("RandomSeed").value))
            simulation.SetRunNumber(int(pc.get_parameter("RunNumber").value))
            simulation.SetDomainType(getattr(domain_type, domain_key.upper()))
            simulation.Run()
            pc.save(file_handler.GetOutputDirectoryFullPath() + study["name"] + "/" +
                    local_ext + "/input_parameters.p")

            result = 1
            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            break

    comm.send(None, dest=0, tag=tags.EXIT)