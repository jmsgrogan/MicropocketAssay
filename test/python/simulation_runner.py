
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

import os
import numpy as np
import vtk
import chaste # Core Chaste functionality
chaste.init() # Initialize MPI and PETSc
import chaste.cell_based
import microvessel_chaste.geometry # Geometry tools
import microvessel_chaste.mesh # Meshing
import microvessel_chaste.population.vessel # Vessel tools
import microvessel_chaste.pde # PDE and solvers
import microvessel_chaste.simulation # Flow and angiogenesis solvers
import microvessel_chaste.visualization # Visualization
from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience

      
def initialize_reference_scales():
    
    reference_length = 1.e-6 * metre()
    reference_time = 3600.0 * second()
    reference_concentration = 1.e-6*mole_per_metre_cubed()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    BaseUnits.Instance().SetReferenceTimeScale(reference_time)
    BaseUnits.Instance().SetReferenceConcentrationScale(reference_concentration)
    return [reference_length, reference_time, reference_concentration]
    
def get_domain_dimensions():
    
    domain_dimensions = {"pellet height" : 1.0e-3*metre(),
                         "cornea radius" : 1.3e-3*metre(),
                         "grid spacing" : 40.0e-6*metre(),
                         "node spacing" : 40.0e-6*metre(),
                         "limbal offset" : 200.0e-6*metre(),
                         "density grid spacing" : 40.0e-6*metre(),
                         "sample spacing x" : 20.0e-6*metre(),
                         "sample spacing y" : 20.0e-6*metre(),
                         }
    return domain_dimensions 

def get_angiogenesis_parameters():
    
    angiogenesis_parameters = {"attraction strength" : 0.0,
                         "chemotactic strength" : 0.0,
                         "persistence angle" : 0.0,
                         "sprouting probability" : 0.5 /(3600.0*second()), 
                         "tip exclusion radius" : 40.0e-6*metre(), 
                         "do anastomosis" : True, 
                         }
    return angiogenesis_parameters 

def get_vegf_parameters():
    
    vegf_parameters = {"pellet concentration" : 0.3e-9*mole_per_metre_cubed(),
                         }
    return vegf_parameters 

def get_simulation_parameters():
    
    simulation_parameters = {"total time" : 48, 
                             "time step" : 0.5}
    return simulation_parameters

def get_2d_planar_domain(domain_dimensions, reference_length):  
    
    domain = microvessel_chaste.geometry.Part2()
    domain.AddRectangle(2.0*np.pi*domain_dimensions["cornea radius"], 
                        domain_dimensions["pellet height"], 
                        microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0))
    return domain  

def get_2d_planar_vessel_network(domain_dimensions, reference_length):

    generator = microvessel_chaste.population.vessel.VesselNetworkGenerator2()
    domain_length = 2.0*np.pi*domain_dimensions["cornea radius"]
    divisions = int(float(domain_length/domain_dimensions["node spacing"])) - 2
    alignment_axis = 0 # pointing x direction
    
    network = generator.GenerateSingleVessel(domain_length, 
                                             microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 
                                                                     float(domain_dimensions["limbal offset"]/reference_length), 
                                                                     0.0, reference_length), 
                                             divisions, alignment_axis)
    for eachNode in network.GetNodes():
        eachNode.GetFlowProperties().SetPressure(1.0*pascal())
    return network

def run(work_dir):
    
    # Initialize the simulation
    chaste.cell_based.SimulationTime.Instance().SetStartTime(0.0)
    chaste.core.RandomNumberGenerator.Instance().Reseed(1234)
    file_handler = chaste.core.OutputFileHandler(work_dir+"/Run/", True)
    reference_length, reference_time, reference_concentration = initialize_reference_scales()
    domain_dimensions = get_domain_dimensions()
    
    # Set up the domain
    domain = get_2d_planar_domain(domain_dimensions, reference_length)
    grid = microvessel_chaste.mesh.RegularGrid2()
    grid.GenerateFromPart(domain, domain_dimensions["grid spacing"])

    # Prescribe a linearly increasing vegf field using a function map
    funciton_map = microvessel_chaste.pde.FunctionMap2()
    funciton_map.SetGrid(grid)
    
    vegf_parameters = get_vegf_parameters()
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        dimless_pellet_height = (domain_dimensions["pellet height"]/reference_length)
        normalized_distance = y_loc/dimless_pellet_height
        vegf_field.append(normalized_distance*(vegf_parameters["pellet concentration"]/reference_concentration))
    funciton_map.SetFileHandler(file_handler)
    funciton_map.SetFileName("Vegf_Field")
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    funciton_map.Write()
    
    # Set up the limbal vessel
    network = get_2d_planar_vessel_network(domain_dimensions, reference_length)
    
    # Set up the angiogenesis problem
    angiogenesis_parameters = get_angiogenesis_parameters()
    migration_rule = microvessel_chaste.simulation.OffLatticeMigrationRule2()
    migration_rule.SetDiscreteContinuumSolver(funciton_map)
    migration_rule.SetNetwork(network)
    migration_rule.SetAttractionStrength(angiogenesis_parameters["attraction strength"])
    migration_rule.SetChemotacticStrength(angiogenesis_parameters["chemotactic strength"])
    migration_rule.SetPersistenceAngleSdv((angiogenesis_parameters["persistence angle"]/180.0)*np.pi)
    
    sprouting_rule = microvessel_chaste.simulation.OffLatticeSproutingRule2()
    sprouting_rule.SetDiscreteContinuumSolver(funciton_map)
    sprouting_rule.SetVesselNetwork(network)
    sprouting_rule.SetSproutingProbability(angiogenesis_parameters["sprouting probability"])
    sprouting_rule.SetOnlySproutIfPerfused(True)
    sprouting_rule.SetTipExclusionRadius(angiogenesis_parameters["tip exclusion radius"])
    
    angiogenesis_solver = microvessel_chaste.simulation.AngiogenesisSolver2()
    angiogenesis_solver.SetVesselNetwork(network)
    angiogenesis_solver.SetMigrationRule(migration_rule)
    angiogenesis_solver.SetSproutingRule(sprouting_rule)
    angiogenesis_solver.SetOutputFileHandler(file_handler)
    angiogenesis_solver.SetBoundingDomain(domain)
    angiogenesis_solver.SetDoAnastomosis(angiogenesis_parameters["do anastomosis"])
    
    simulation_parameters = get_simulation_parameters()
    num_steps = int(simulation_parameters["total time"]/simulation_parameters["time step"])
    chaste.cell_based.SimulationTime.Instance().SetEndTimeAndNumberOfTimeSteps(simulation_parameters["total time"], 
                                                                               num_steps)
    network_writer = microvessel_chaste.population.vessel.VesselNetworkWriter2()
    denisty_map_grid = microvessel_chaste.mesh.RegularGrid2()
    denisty_map_grid.GenerateFromPart(domain, domain_dimensions["density grid spacing"])
    num_sample_points_x = int(2.0*np.pi*float(domain_dimensions["cornea radius"]/domain_dimensions["sample spacing x"])) + 1
    num_sample_points_y = int(float((domain_dimensions["pellet height"]-domain_dimensions["limbal offset"])/domain_dimensions["sample spacing y"])) + 1

    f = open(file_handler.GetOutputDirectoryFullPath()+"sampled_line_density.txt", "w")
    f.write("Time, ")
    
    sample_lines = []
    for idx in range(num_sample_points_y):
        sample_points = vtk.vtkPoints()
        dimless_sample_spacing_x = domain_dimensions["sample spacing x"]/reference_length
        dimless_sample_spacing_y = domain_dimensions["sample spacing y"]/reference_length
        dimless_limbal_offset = domain_dimensions["limbal offset"]/reference_length
        for jdx in range(num_sample_points_x):
            sample_points.InsertNextPoint(float(jdx)*dimless_sample_spacing_x,
                                          float(idx)*dimless_sample_spacing_y + dimless_limbal_offset,
                                          0.0)
        sample_lines.append(sample_points)
        f.write(str(float(idx)*float(dimless_sample_spacing_y))+",")        
    f.write("\n")
    
    while not chaste.cell_based.SimulationTime.Instance().IsFinished():
        
        density_map = microvessel_chaste.pde.DensityMap2()
        density_map.SetGrid(denisty_map_grid)
        density_map.SetVesselNetwork(network)

        density_map_result = microvessel_chaste.pde.FunctionMap2()
        density_map_result.SetGrid(denisty_map_grid)
        density_map_result.SetVesselNetwork(network)
        density_map_result.SetFileHandler(file_handler)
        elapsed_time = chaste.cell_based.SimulationTime.Instance().GetTimeStepsElapsed()
        density_map_result.SetFileName("/line_density" + str(elapsed_time))
        density_map_result.UpdateSolution(density_map.rGetVesselLineDensity())
        density_map_result.Write()

        # Sample the density map
        time = chaste.cell_based.SimulationTime.Instance().GetTime()
        f.write(str(time) + ",")
        for eachLine in sample_lines:
            solution = density_map_result.GetSolutionP(eachLine)
            mean = np.mean(np.array(solution))
            f.write(str(mean) + ",")              
        f.write("\n")
        f.flush()

        density_map_result.SetFileName("/tip_density" + str(elapsed_time))
        density_map_result.UpdateSolution(density_map.rGetVesselTipDensity())
        density_map_result.Write()

        density_map_result.SetFileName("/branch_density" + str(elapsed_time))
        density_map_result.UpdateSolution(density_map.rGetVesselBranchDensity())
        density_map_result.Write()

        network_writer.SetFileName(
                file_handler.GetOutputDirectoryFullPath() + "/vessel_network_" + str(elapsed_time) + ".vtp")
        network_writer.SetVesselNetwork(network)
        network_writer.Write()

        # Increment the solver and simulation time
        angiogenesis_solver.Increment()
        chaste.cell_based.SimulationTime.Instance().IncrementTimeOneStep()
               
    f.close()

    chaste.cell_based.SimulationTime.Instance().Destroy()
        
if __name__ == '__main__':
    
    work_dir = "Python/Cornea/one_dimensional"
    run(work_dir)

    