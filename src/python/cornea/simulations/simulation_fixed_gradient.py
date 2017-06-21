import chaste # Core Chaste functionality
chaste.init() # Initialize MPI and PETSc

import microvessel_chaste.simulation
import cornea.parameters.default_parameters
        
if __name__ == '__main__':
    
    work_dir = "Python/Cornea/TestSimulationFixedGradient/"
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    parameter_collection.get_parameter("UseFixedGradient").value = True
    
    domain_types = ["Planar_2D", "Circle_2D", "Planar_3D", "Circle_3D", "Hemisphere_3D"]
    #domain_types = ["Planar 2D"]
    random_seeds = [1234]

    for eachDomainType in domain_types:
        run_number = 0
        parameter_collection.get_parameter("DomainType").value = eachDomainType
        for eachSeed in random_seeds:
            parameter_collection.get_parameter("RunNumber").value = run_number
            parameter_collection.get_parameter("RandomSeed").value = eachSeed
            local_work_dir = work_dir + "/DomainType_" + eachDomainType.replace(" ", "") + "/Run_" + str(run_number)
            
            if "2" in eachDomainType:
                simulation = microvessel_chaste.simulation.CornealMicropocketSimulation2()
            else 
                simulation = microvessel_chaste.simulation.CornealMicropocketSimulation3()
            simulation.SetWorkDir(local_work_dir)
            for eachParameter in parameter_collection.keys():
                if eachParameter not in ["DomaintType"]:
                    getattr(simulation, 'Set'+parameter_collection[eachParameter].name)(parameter_collection[eachParameter].value)
                    
            simulation.SetDomainType(getattr(microvessel_chaste.simulation, "DomainType." + eachDomainType))

            simulation.Run()
            run_number += 1