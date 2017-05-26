import time
import chaste
chaste.init()
import microvessel_chaste
from petsc4py import PETSc
from mpi4py import MPI
from argparse import ArgumentParser

import parameter_collection

def run_simulation(work_dir):
    
    intracom = PETSc.COMM_WORLD.tompi4py()
    intercom = intracom.Get_parent()

    print "Entering LP: ", intracom.Get_rank(), " on GP: ", intercom.Get_rank(), " at: ", work_dir
    parameters = parameter_collection.SimulationParameterCollection()
    parameters.load(work_dir+"/input_parameters.p")
    print parameters.get_parameter("sprouting probability").value
    print "Exiting LP: ", intracom.Get_rank(), " on GP: ", intercom.Get_rank()

if __name__=="__main__":
    
    parser = ArgumentParser()
    parser.add_argument('-i', type=str)
    args = parser.parse_args()
    work_dir = args.i.replace(" ", "")
    run_simulation(work_dir)