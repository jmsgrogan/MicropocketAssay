import uuid
import pickle
import chaste  # Core Chaste functionality
from microvessel_chaste.utility import *
import microvessel_chaste.simulation
import cornea.parameters.default_parameters
chaste.init()  # Initialize MPI and PETSc

if __name__ == '__main__':

    studies = []
    studies.append({"name": "fixed_gradient_10",
                    "switches": {"UseFixedGradient": False,
                                 "PelletConcentration": 100.0*mole_per_metre_cubed,
                                 "PelletHeight": 0.7e-3*metres,
                                 "ChemotacticStrength": 1.0}})
#     studies.append({"name": "fixed_gradient_07",
#                     "switches": {"UseFixedGradient": False,
#                                  "PelletConcentration": 100.0*mole_per_metre_cubed,
#                                  "PelletHeight": 0.7e-3*metres,
#                                  "ChemotacticStrength": 1.0}})
#     studies.append({"name": "fixed_gradient_10",
#                     "switches": {"UseFixedGradient": True,
#                                  "PelletConcentration": 1.5e-10*mole_per_metre_cubed,
#                                  "PersistenceAngle": 0.0,
#                                  "ChemotacticStrength": 1.0,
#                                  "OnlyPerfusedSprout": True}})
#     studies.append({"name": "fixed_gradient_15",
#                     "switches": {"UseFixedGradient": True,
#                                  "PelletConcentration": 2.0e-10*mole_per_metre_cubed,
#                                  "PersistenceAngle": 0.0,
#                                  "ChemotacticStrength": 1.0,
#                                  "OnlyPerfusedSprout": True}})
#     studies.append({"name": "fixed_gradient_20",
#                     "switches": {"UseFixedGradient": True,
#                                  "PelletConcentration": 2.5e-10*mole_per_metre_cubed,
#                                  "PersistenceAngle": 0.0,
#                                  "ChemotacticStrength": 1.0,
#                                  "OnlyPerfusedSprout": True}})
    restart = {'run_id': "c3e35591-035d-4eae-a83b-e40a9422af3f",
               'study_id': 0}
    restart = None

    if restart is None:
        run_id = uuid.uuid4()
    else:
        run_id = restart["run_id"]
        studies = studies[restart["study_id"]:]

    work_dir = "Python/Cornea/Study_pde_vary_h" + str(run_id) + "/"
    random_seeds = [1234, 5678, 9101112, 5745745, 235235645]
    random_seeds = [1234]
    domain_types = ["Planar_2D", "Circle_2D", "Planar_3D",
                    "Circle_3D", "Hemisphere"]
    domain_types = ["Planar_3D", "Planar_3D_Finite"]

    study_names = [x["name"] for x in studies]
    study_data = {"random_seeds": random_seeds,
                  "domain_types": domain_types,
                  "study_names": study_names}
    file_handler = chaste.core.OutputFileHandler(work_dir, False)
    if restart is None:
        pickle.dump(study_data, open(file_handler.GetOutputDirectoryFullPath() +
                                     "/study_data.p", "wb"))

    for eachStudy in studies:
        pc = cornea.parameters.default_parameters.get_default_collection()
        int_work_dir = work_dir + "/" + eachStudy["name"]
        switches = eachStudy["switches"]
        for key, value in switches.iteritems():
            pc.get_parameter(key).value = value
        run_number = 0
        for eachSeed in random_seeds:
            for eachDomainType in domain_types:
                h = pc.get_parameter("PelletHeight").value.Convert(1.0*metres)
                v = pc.get_parameter("TipVelocity").value.Convert(1.0*metre_per_second)
                t = h/v
                pc.get_parameter("TotalTime").value = 3600.0*round(0.9*t/3600.0)*seconds
                print pc.get_parameter("TotalTime").value.Convert(seconds)
                pc.get_parameter("SampleSpacingX").value = 30.0e-6*metres
                pc.get_parameter("DomainType").value = eachDomainType
                pc.get_parameter("RunNumber").value = run_number
                pc.get_parameter("RandomSeed").value = eachSeed
                if "Finite" in eachDomainType:
                    pc.get_parameter("FinitePelletWidth").value = True
                domain_string = "/DomainType_" + eachDomainType.replace(" ", "")
                run_string = "/Run_" + str(run_number) + "/"
                local_ext = domain_string + run_string
                local_work_dir = int_work_dir + local_ext
                if "2" in eachDomainType:
                    simulation = microvessel_chaste.simulation.CornealMicropocketSimulation2()
                else: 
                    simulation = microvessel_chaste.simulation.CornealMicropocketSimulation3()
                simulation.SetWorkDir(local_work_dir)
                for eachParameter in pc.collection.keys():
                    if eachParameter not in ["DomainType"]:
                        param_nam = pc.collection[eachParameter].name
                        param_value = pc.collection[eachParameter].value
                        getattr(simulation, 'Set'+param_nam)(param_value)
                domain_type = microvessel_chaste.simulation.DomainType
                domain_key = eachDomainType
                if "Finite" in eachDomainType:
                    domain_key = eachDomainType[0:-7]
                simulation.SetDomainType(getattr(domain_type, domain_key.upper()))
                simulation.Run()
                pc.save(file_handler.GetOutputDirectoryFullPath() + eachStudy["name"] + "/" +
                        local_ext + "/input_parameters.p")
            run_number += 1