import time
import chaste
chaste.init()
from petsc4py import PETSc
from mpi4py import MPI
from argparse import ArgumentParser

def run_simulation(value):
    
    print "Enter ", value
    print PETSc.COMM_WORLD
    cw = PETSc.COMM_WORLD.tompi4py()
    print isinstance(cw, MPI.Intracomm)
    print cw.Get_size()
    print cw.Get_rank()
    print cw.Get_parent()
    time.sleep(value)
    print "Exit ", value

if __name__=="__main__":
    
    parser = ArgumentParser()
    parser.add_argument('-i', type=int)
    args = parser.parse_args()
    value = args.i
    run_simulation(value)