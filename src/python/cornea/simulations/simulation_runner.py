import time
import chaste
chaste.init()
import microvessel_chaste
from petsc4py import PETSc
from mpi4py import MPI
from argparse import ArgumentParser

import cornea.parameters.parameter_collection
import cornea.simulations.base_simulation

def run_simulation(work_dir):
    
    intracom = PETSc.COMM_WORLD.tompi4py()
    intercom = intracom.Get_parent()

    print "Entering on local: ", intracom.Get_rank(), " on global: ", intercom.Get_rank(), " at: ", work_dir
    print "Proc name: ", MPI.Get_processor_name()
    local_parameters = cornea.parameters.parameter_collection.SimulationParameterCollection()
    local_parameters.load(work_dir+"/input_parameters.p")
    
    simulation = cornea.simulations.base_simulation.BaseSimulation(local_parameters, work_dir)
    
#    time.sleep(10)
    try:
        simulation.run()
    except Exception as e:
        print e.__doc__
        print e.message
    
    print "Exiting on local: ", intracom.Get_rank(), " on global: ", intercom.Get_rank()
    intercom.send(1, dest=0, tag=1234)
    intercom.Disconnect()

if __name__=="__main__":
    
    parser = ArgumentParser()
    parser.add_argument('-i', type=str)
    args = parser.parse_args()
    work_dir = args.i.replace(" ", "")
    run_simulation(work_dir)