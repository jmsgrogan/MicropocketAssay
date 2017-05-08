
"""Copyright (c) 2005-2017, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import numpy as np
import chaste # Core Chaste functionality
chaste.init() # Initialize MPI and PETSc
import chaste.cell_based
import microvessel_chaste.mesh # Meshing
import microvessel_chaste.population.vessel # Vessel tools
import microvessel_chaste.pde # PDE and solvers
import microvessel_chaste.simulation # Flow and angiogenesis solvers
from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience

import geometries.domains
import geometries.grid
import geometries.networks
import geometries.sample_points
import geometries.pdes
      
def initialize_reference_scales():
    
    reference_length = 1.e-6 * metre()
    reference_time = 3600.0 * second()
    reference_concentration = 1.e-6*mole_per_metre_cubed()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    BaseUnits.Instance().SetReferenceTimeScale(reference_time)
    BaseUnits.Instance().SetReferenceConcentrationScale(reference_concentration)
    return [reference_length, reference_time, reference_concentration]

def run(parameter_collection, work_dir, domain_type, random_seed, run_number):
    
    # Initialize the simulation
    domain_dimensions, angiogenesis_parameters, pde_parameters, simulation_parameters = parameter_collection
    chaste.cell_based.SimulationTime.Instance().SetStartTime(0.0)
    chaste.core.RandomNumberGenerator.Instance().Reseed(random_seed)
    file_handler = chaste.core.OutputFileHandler(work_dir+"/Run" + str(run_number)+ "/", True)
    reference_length, _, reference_concentration = initialize_reference_scales()
    
    # Set up the domain
    domain, holes = geometries.domains.get_domain(domain_type, domain_dimensions, reference_length)
    grid = geometries.grid.get_grid(domain_type, domain, domain_dimensions, holes)

    # Prescribe a linearly increasing vegf field using a function map
    pde = geometries.pdes.get_transient_pde(domain_type, domain_dimensions, pde_parameters, 
                                            reference_length, reference_concentration)
    
    solver = geometries.pdes.get_pde_solver(domain_type, pde ,grid, domain, domain_dimensions, pde_parameters, 
                                                          reference_length, reference_concentration)
    solver.SetFileHandler(file_handler)
    solver.SetWriteSolution(True)
    solver.SetTargetTimeIncrement(pde_parameters["pde time increment"]) # hours
    solver.SetUseCoupling(True)
    solver.SetWriteIntermediateSolutions(False, int(2.0/pde_parameters["pde time increment"])) # every 2 hours
    
    # Set up the limbal vessel
    network = geometries.networks.get_vessel_network(domain_type, domain_dimensions, reference_length)
    
    # Set up the angiogenesis problem
    if "2" in domain_type:
        dimension = str(2)
    else:
        dimension = str(3)
    migration_rule = getattr(microvessel_chaste.simulation, 'OffLatticeMigrationRule' + dimension)()
    migration_rule.SetDiscreteContinuumSolver(solver)
    migration_rule.SetNetwork(network)
    migration_rule.SetAttractionStrength(angiogenesis_parameters["attraction strength"])
    migration_rule.SetChemotacticStrength(angiogenesis_parameters["chemotactic strength"])
    migration_rule.SetPersistenceAngleSdv((angiogenesis_parameters["persistence angle"]/180.0)*np.pi)
    
    sprouting_rule = getattr(microvessel_chaste.simulation, 'OffLatticeSproutingRule' + dimension)()
    sprouting_rule.SetDiscreteContinuumSolver(solver)
    sprouting_rule.SetVesselNetwork(network)
    sprouting_rule.SetSproutingProbability(angiogenesis_parameters["sprouting probability"])
    sprouting_rule.SetOnlySproutIfPerfused(True)
    sprouting_rule.SetTipExclusionRadius(angiogenesis_parameters["tip exclusion radius"])
    
    angiogenesis_solver = getattr(microvessel_chaste.simulation, 'AngiogenesisSolver' + dimension)()
    angiogenesis_solver.SetVesselNetwork(network)
    angiogenesis_solver.SetMigrationRule(migration_rule)
    angiogenesis_solver.SetSproutingRule(sprouting_rule)
    angiogenesis_solver.SetOutputFileHandler(file_handler)
    angiogenesis_solver.SetBoundingDomain(domain)
    angiogenesis_solver.SetDoAnastomosis(angiogenesis_parameters["do anastomosis"])
    

    num_steps = int(simulation_parameters["total time"]/simulation_parameters["time step"])
    chaste.cell_based.SimulationTime.Instance().SetEndTimeAndNumberOfTimeSteps(simulation_parameters["total time"], 
                                                                               num_steps)
    network_writer = getattr(microvessel_chaste.population.vessel, 'VesselNetworkWriter' + dimension)()
    
    if domain_type == "Planar 2D":
        denisty_map_grid = microvessel_chaste.mesh.RegularGrid2()
        denisty_map_grid.GenerateFromPart(domain, domain_dimensions["density grid spacing"])
    elif domain_type == "Planar 3D":
        denisty_map_grid = microvessel_chaste.mesh.RegularGrid3()
        denisty_map_grid.GenerateFromPart(domain, domain_dimensions["density grid spacing"])
    else:
        denisty_map_grid = grid

    sample_lines = geometries.sample_points.get_sample_points(domain_type, domain_dimensions, reference_length)
    output_density_quantities = ["Line", "Branch", "Tip"]
    output_density_files = []
    for eachQuantity in output_density_quantities:
        output_file = open(file_handler.GetOutputDirectoryFullPath()+"Sampled_" + eachQuantity + "_density.txt", "w")
        output_file.write("Time, ")
        for eachLine in sample_lines:
            y_coord = eachLine.GetPoint(0)[1]
            output_file.write(str(y_coord)+",")        
        output_file.write("\n")
        output_density_files.append(output_file)
        
    # PDE Sampling
    pde_output_file = open(file_handler.GetOutputDirectoryFullPath()+"Sampled_PDE.txt", "w")
    pde_output_file.write("Time, ")
    for eachLine in sample_lines:
        y_coord = eachLine.GetPoint(0)[1]
        pde_output_file.write(str(y_coord)+",")        
    pde_output_file.write("\n") 
        
    old_time = chaste.cell_based.SimulationTime.Instance().GetTime()
    while not chaste.cell_based.SimulationTime.Instance().IsFinished():
        
        elapsed_time = chaste.cell_based.SimulationTime.Instance().GetTimeStepsElapsed()
        time = chaste.cell_based.SimulationTime.Instance().GetTime()
        solver.SetStartTime(old_time)
        if time==old_time:
            solver.SetEndTime(time+pde_parameters["pde time increment"])
        else:
            solver.SetEndTime(time)
        solver.SetFileName("Vegf_Solution" + str(elapsed_time))
        try:
            solver.Solve()
        except chaste.ChasteException as e:
            print e.GetMessage
        
        # Sample the pde
        pde_output_file.write(str(time) + ",")
        for eachLine in sample_lines:
            solution = solver.GetSolutionP(eachLine)
            mean = np.mean(np.array(solution))
            pde_output_file.write(str(mean) + ",")              
        pde_output_file.write("\n")
        pde_output_file.flush()
        
        density_map = getattr(microvessel_chaste.pde, 'DensityMap' + dimension)()
        density_map.SetGrid(denisty_map_grid)
        density_map.SetVesselNetwork(network)

        density_map_result = getattr(microvessel_chaste.pde, 'FunctionMap' + dimension)()
        density_map_result.SetGrid(denisty_map_grid)
        density_map_result.SetVesselNetwork(network)
        density_map_result.SetFileHandler(file_handler)
        
        for idx, eachQuantity in enumerate(output_density_quantities):
            
            density_map_result.SetFileName("/" + eachQuantity + "_Density" + str(elapsed_time))
            calculation_method = getattr(density_map, 'rGetVessel' + eachQuantity + 'Density')
            if "Planar" in domain_type:
                density_map_result.UpdateSolution(calculation_method())
            else:
                density_map_result.UpdateElementSolution(calculation_method())
            density_map_result.Write()

            # Sample the density map
            output_density_files[idx].write(str(time) + ",")
            for eachLine in sample_lines:
                solution = density_map_result.GetSolutionP(eachLine)
                mean = np.mean(np.array(solution))
                output_density_files[idx].write(str(mean) + ",")              
            output_density_files[idx].write("\n")
            output_density_files[idx].flush()

        network_writer.SetFileName(
                file_handler.GetOutputDirectoryFullPath() + "/vessel_network_" + str(elapsed_time) + ".vtp")
        network_writer.SetVesselNetwork(network)
        network_writer.Write()

        # Increment the solver and simulation time
        angiogenesis_solver.Increment()
        chaste.cell_based.SimulationTime.Instance().IncrementTimeOneStep()
        old_time = time
            
    for eachFile in output_density_files:    
        eachFile.close()
    pde_output_file.close()

    chaste.cell_based.SimulationTime.Instance().Destroy()
        
if __name__ == '__main__':
    
    domain_dimensions = {"pellet height" : 1.0e-3*metre(),
                     "cornea radius" : 1.3e-3*metre(),
                     "cornea thickness" : 100.0e-6*metre(),
                     "pellet thickness" : 40.0e-6*metre(),
                     "grid spacing" : 40.0e-6*metre(),
                     "node spacing" : 40.0e-6*metre(),
                     "limbal offset" : 200.0e-6*metre(),
                     "density grid spacing" : 40.0e-6*metre(),
                     "sample spacing x" : 20.0e-6*metre(),
                     "sample spacing y" : 20.0e-6*metre(),
                     "sample spacing z" : 20.0e-6*metre(),
                     "use pellet" : True,
                     "use finite pellet width" : True,
                     "pellet radius" : 300.0e-6*metre(),
                     }
    
    angiogenesis_parameters = {"attraction strength" : 0.0,
                         "chemotactic strength" : 0.0,
                         "persistence angle" : 5.0,
                         "sprouting probability" : 0.5 /(3600.0*second()), 
                         "tip exclusion radius" : 40.0e-6*metre(), 
                         "do anastomosis" : True, 
                         }
    
    pde_parameters = {"pellet concentration" : 0.3e-9*mole_per_metre_cubed(),
                      "vegf diffusivity" : 6.94e-11 * metre_squared_per_second(),
                      "vegf decay rate" : (-0.8/3600.0) * per_second(),
                      "pde time increment" : 0.1,
                      "pellet binding constant" : 100.0,
                      "include vessel sink" : False,
                      "vegf blood concentration" : 0.0*mole_per_metre_cubed(),
                      "vessel permeability" : (3.e-4/3600.0)*metre_per_second(),
                      "uptake rate per cell" : (4.e-22/3600.0)*mole_per_second()
                      }
    
#     simulation_parameters = {"total time" : 48, 
#                              "time step" : 0.5}
    simulation_parameters = {"total time" : 24, 
                             "time step" : 0.5}
    
    parameter_collection = [domain_dimensions, angiogenesis_parameters, pde_parameters, simulation_parameters] 
    
    work_dir = "Python/Cornea/TestSimulationPde/"
    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    domain_types = ["Circle 3D"]
    random_seeds = [1234]
    
    for eachDomainType in domain_types:
        run_number = 0
        for eachSeed in random_seeds:
            run(parameter_collection, work_dir + eachDomainType.replace(" ", ""), eachDomainType, eachSeed, run_number)
            run_number += 1
    