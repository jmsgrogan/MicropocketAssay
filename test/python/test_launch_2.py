#!/usr/bin/env python
from mpi4py import MPI
import chaste
import sys
import os

rank = MPI.COMM_WORLD.Get_rank()
new_comm = MPI.COMM_WORLD.Split(color=rank, key=rank)
new_rank = new_comm.Get_rank()

# cwd=os.getcwd()
# os.mkdir(str(rank))
# directory=os.path.join(cwd,str(rank))
# print(rank,directory)
# os.chdir(directory)

print "Launch rank:", rank
new_comm.Spawn("python", args=["test_launch_script.py", "-i "+str(rank)], maxprocs=1)
print "return rank:", rank