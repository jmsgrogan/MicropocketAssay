import uuid
import pickle
import chaste  # Core Chaste functionality
from microvessel_chaste.utility import *
import cornea.parameters.default_parameters

study_list = []
cp = [100.0,] # nM

for idx in range(len(cp)):
    dimless_c = cp[idx]
    study_list.append({"name": "fg_cp_"+str(int(round(dimless_c))),
                       "switches": {"UseFixedGradient": True,
                                    "PelletConcentration": dimless_c*1.e-6*mole_per_metre_cubed,
                                    "ChemotacticStrength": 1.0,
                                    "PersistenceAngle": 0.0,
                                    "SproutingProbability": (1000.0 /3600.0)*per_second,
                                    "OnlyPerfusedSprout": True,
                                    "FinitePelletWidth": True,
                                    "DoAnastamosis": False}})

run_id = uuid.uuid4()
master_work_dir = "Python/Cornea/Study_full_activation" + str(run_id) + "/"
random_seeds = [1234]
domains = ["Planar_2D", "Circle_2D", "Planar_3D", "Circle_3D", "Hemisphere"]
domains = ["Planar_2D"]
study_names = [x["name"] for x in study_list]
study_data = {"random_seeds": random_seeds,
              "domain_types": domains,
              "study_names": study_names}

study = []
counter = 0
for eachStudy in study_list:
    run_number = 0
    for eachSeed in random_seeds:
        for eachDomainType in domains:
            pc = cornea.parameters.default_parameters.get_default_collection()
            switches = eachStudy["switches"]
            for key, value in switches.iteritems():
                pc.get_parameter(key).value = value
            h = pc.get_parameter("PelletHeight").value.Convert(1.0*metres)
            v = pc.get_parameter("TipVelocity").value.Convert(1.0*metre_per_second)
            t = h/v
            pc.get_parameter("TotalTime").value = 3600.0*round(0.9*t/3600.0)*seconds
            pc.get_parameter("SampleSpacingX").value = 180.0e-6*metres
            pc.get_parameter("DomainType").value = eachDomainType
            pc.get_parameter("RunNumber").value = run_number
            pc.get_parameter("RandomSeed").value = int(eachSeed)
            if "Finite" in eachDomainType:
                pc.get_parameter("FinitePelletWidth").value = True
            study.append([eachStudy, run_number, eachDomainType, pc, master_work_dir])
        run_number += 1