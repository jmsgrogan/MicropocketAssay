import time
import petsc4py
petsc4py.init()
from petsc4py import PETSc
from mpi4py import MPI
from argparse import ArgumentParser

def run_simulation(work_dir):
    
    intracom = PETSc.COMM_WORLD.tompi4py()
    intercom = intracom.Get_parent()

    print "Entering on local: ", intracom.Get_rank(), " on global: ", intercom.Get_rank(), " at: ", work_dir
    print "Proc name: ", MPI.Get_processor_name()
    
    time.sleep(10)

    print "Exiting on local: ", intracom.Get_rank(), " on global: ", intercom.Get_rank()
    intercom.send(1, dest=0, tag=1234)
    intercom.Disconnect()

if __name__=="__main__":
    
    parser = ArgumentParser()
    parser.add_argument('-i', type=str)
    args = parser.parse_args()
    work_dir = args.i.replace(" ", "")
    run_simulation(work_dir)