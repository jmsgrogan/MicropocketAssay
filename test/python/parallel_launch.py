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
from mpi4py import MPI

n_workers = 16
start_worker = 'worker'
usage = 'Program should be started without argument'

# Parent
if len(sys.argv) == 1:

    # Start clock
    start = MPI.Wtime()

    # Don't spawn until this is done
    time.sleep(2)

    # Append stop sentinel for each worker
    study = range(100)
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

    # Connect to parent
    try:
        comm = MPI.Comm.Get_parent()
        rank = comm.Get_rank()
    except:
        raise ValueError('Could not connect to parent - ' + usage)

    # Ask for work until stop sentinel
    local_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
    print "Worker rank", local_comm.Get_rank()

    log = []
    for task in iter(lambda: comm.sendrecv(dest=0), StopIteration):
        print "Worker task", task
        time.sleep(2)

        log.append(task)

    # Collective report to parent
    comm.gather(sendobj=log, root=0)

    # Shutdown
    comm.Disconnect()

# Catch
else:
    raise ValueError(usage)