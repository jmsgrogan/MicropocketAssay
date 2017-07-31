#!/usr/bin/env python

""" Parallel batch run based on an example by jborschein
(https://github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py)
Modified to allow launching of multiple single process PETSc jobs via
petsc4py.
"""

import os
import time
import sys
import pickle
import chaste
import types
import microvessel_chaste.simulation
from microvessel_chaste.utility import *
from cornea.simulations.pde_vary_h import study
from cornea.simulations.pde_vary_h import master_work_dir
from cornea.simulations.pde_vary_h import study_data

n_workers = 3
n_tasks = 50
start_worker = 'worker'
usage = 'Program should be started without argument'

# Parent
if len(sys.argv) == 1:

    from mpi4py import MPI

    # Start clock
    start = MPI.Wtime()

    # Set up the problem, avoid initializing Chaste
    test_output = os.getenv('CHASTE_TEST_OUTPUT', os.getcwd())
    if not os.path.exists(test_output + "/" + master_work_dir):
        os.makedirs(test_output + "/" + master_work_dir)
    pickle.dump(study_data, open(test_output + "/" + master_work_dir + "/study_data.p", "wb"))

    # Don't spawn until this is done
    time.sleep(2)

    # Append stop sentinel for each worker
    msg_list = study + ([StopIteration] * n_workers)

    # Spawn workers
    comm = MPI.COMM_WORLD.Spawn(
        sys.executable,
        args=[sys.argv[0], start_worker],
        maxprocs=n_workers)

    # Reply to whoever asks until done
    status = MPI.Status()
    for position, msg in enumerate(msg_list):
        comm.recv(source=MPI.ANY_SOURCE, status=status)
        comm.send(obj=msg, dest=status.Get_source())

        # Simple (loop position) progress bar
        sys.stdout.write('\rProgress:  %3i \n' % (position))
        sys.stdout.flush()

    # Gather reports from workers
    reports = comm.gather(root=MPI.ROOT)

    # Final statistics
    finish = MPI.Wtime() - start
    print('\nProcessed in %.2f secs' % finish)

    # Shutdown
    comm.Disconnect()

# Worker
elif sys.argv[1] == start_worker:

    from mpi4py import MPI

    # Connect to parent
    try:
        comm = MPI.Comm.Get_parent()
        rank = comm.Get_rank()
    except:
        raise ValueError('Could not connect to parent - ' + usage)

    # Ask for work until stop sentinel
    local_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
    chaste.init(comm=local_comm)

    log = []
    for task in iter(lambda: comm.sendrecv(dest=0), StopIteration):
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
        log.append(task)

    # Collective report to parent
    comm.gather(sendobj=log, root=0)

    # Shutdown
    comm.Disconnect()

# Catch
else:
    raise ValueError(usage)