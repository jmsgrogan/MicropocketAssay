import time
import chaste.core
chaste.init()
import microvessel_chaste
import microvessel_chaste.simulation
from petsc4py import PETSc
from mpi4py import MPI
from argparse import ArgumentParser

import cornea.parameters.parameter_collection

def run_simulation(work_dir):
    
    intracom = PETSc.COMM_WORLD.tompi4py()
    intercom = intracom.Get_parent()

    print "Entering on local: ", intracom.Get_rank(), " on global: ", intercom.Get_rank(), " at: ", work_dir
    print "Proc name: ", MPI.Get_processor_name()
    
    file_handler = chaste.core.OutputFileHandler(work_dir, False)
    local_parameters = cornea.parameters.parameter_collection.SimulationParameterCollection()
    local_parameters.load(file_handler.GetOutputDirectoryFullPath()+"/input_parameters.p")
    
    domain_type = local_parameters.get_parameter("DomainType").value
    if "2" in domain_type:
        simulation = microvessel_chaste.simulation.CornealMicropocketSimulation2()
    else: 
        simulation = microvessel_chaste.simulation.CornealMicropocketSimulation3()
    simulation.SetWorkDir(work_dir)
    local_parameters.collection["RandomSeed"].value = int(local_parameters.collection["RandomSeed"].value)
    local_parameters.collection["RunNumber"].value = int(local_parameters.collection["RunNumber"].value)
    
    for eachParameter in local_parameters.collection.keys():
        if eachParameter not in ["DomainType"]:
            param_nam = local_parameters.collection[eachParameter].name
            param_value = local_parameters.collection[eachParameter].value
            getattr(simulation, 'Set'+param_nam)(param_value)
                    
            domain_type_enum = microvessel_chaste.simulation.DomainType
            simulation.SetDomainType(getattr(domain_type_enum, domain_type.upper()))
    
#    time.sleep(10)
    try:
        simulation.Run()
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