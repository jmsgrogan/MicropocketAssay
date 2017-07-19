import chaste # Core Chaste functionality
chaste.init() # Initialize MPI and PETSc

from microvessel_chaste.utility import *
import microvessel_chaste.simulation
import cornea.parameters.default_parameters

if __name__ == '__main__':

    work_dir = "Python/Cornea/TestSimulationFixedGradient/"
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    parameter_collection.get_parameter("UseFixedGradient").value = True
    parameter_collection.get_parameter("PelletConcentration").value = 1.e-10*mole_per_metre_cubed
    parameter_collection.get_parameter("TotalTime").value = 72.0*3600.0*seconds
    #parameter_collection.get_parameter("PersistenceAngle").value = 0.0
    #parameter_collection.get_parameter("ChemotacticStrength").value = 0.0
    parameter_collection.get_parameter("SampleSpacingX").value = 30.0e-6*metres
    parameter_collection.get_parameter("OnlyPerfusedSprout").value = True
    domain_types = ["Planar_2D", "Circle_2D", "Planar_3D", "Circle_3D", "Hemisphere"]
    domain_types = ["Hemisphere"]
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
            else: 
                simulation = microvessel_chaste.simulation.CornealMicropocketSimulation3()
            simulation.SetWorkDir(local_work_dir)
            for eachParameter in parameter_collection.collection.keys():
                if eachParameter not in ["DomainType"]:
                    param_nam = parameter_collection.collection[eachParameter].name
                    param_value = parameter_collection.collection[eachParameter].value
                    getattr(simulation, 'Set'+param_nam)(param_value)
            domain_type = microvessel_chaste.simulation.DomainType
            simulation.SetDomainType(getattr(domain_type, eachDomainType.upper()))

            simulation.Run()
            run_number += 1