import chaste # Core Chaste functionality
chaste.init() # Initialize MPI and PETSc

import cornea.simulations.base_simulation
import cornea.parameters.default_parameters
        
if __name__ == '__main__':
    
    work_dir = "Python/Cornea/TestSimulationFixedGradient/"
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    parameter_collection.get_parameter("use fixed gradient").value = True
    
    domain_types = ["Planar 2D", "Circle 2D", "Planar 3D", "Circle 3D", "Hemisphere 3D"]
    #domain_types = ["Planar 2D"]
    random_seeds = [1234]

    for eachDomainType in domain_types:
        run_number = 0
        parameter_collection.get_parameter("domain type").value = eachDomainType
        for eachSeed in random_seeds:
            parameter_collection.get_parameter("run number").value = run_number
            parameter_collection.get_parameter("random seed").value = eachSeed
            local_work_dir = work_dir + "/DomainType_" + eachDomainType.replace(" ", "") + "/Run_" + str(run_number)
            simulation = cornea.simulations.base_simulation.BaseSimulation(parameter_collection, local_work_dir)
            simulation.run()
            run_number += 1