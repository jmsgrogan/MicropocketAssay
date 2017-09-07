import uuid
import pickle
import chaste  # Core Chaste functionality
from microvessel_chaste.utility import *
import cornea.parameters.default_parameters

study_list = []
height = 1.0 # mm

for idx in range(1):
    dimless_height = (height - 0.1*float(idx))
    study_list.append({"name": "pde_h_"+str(int(round(dimless_height*1000.0))),
                       "switches": {"UseFixedGradient": False,
                                    "PelletConcentration": 400.0*mole_per_metre_cubed,
                                    "VegfBindingConstant": 20.0,
                                    "PelletHeight": dimless_height*1e-3*metres}})

run_id = uuid.uuid4()
master_work_dir = "Python/Cornea/Fig8_" + str(run_id) + "/"
random_seeds = [1234]
#random_seeds = [1234]
domains = ["Planar_2D", "Planar_3D", "Planar_2D_Finite", "Circle_2D",
           "Planar_3D_Finite", "Circle_3D", "Hemisphere"]

domains = ["Planar_2D_Finite", "Planar_3D_Finite"]
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
            pc.get_parameter("TotalTime").value = 360.0*round(0.9*t/3600.0)*seconds
            pc.get_parameter("SampleSpacingX").value = 180.0e-6*metres
            pc.get_parameter("PelletThickness").value = 99.0e-6*metres
            pc.get_parameter("DomainType").value = eachDomainType
            pc.get_parameter("RunNumber").value = run_number
            pc.get_parameter("RandomSeed").value = int(eachSeed)
            if "Finite" in eachDomainType:
                pc.get_parameter("FinitePelletWidth").value = True
            study.append([eachStudy, run_number, eachDomainType, pc, master_work_dir])
        run_number += 1
