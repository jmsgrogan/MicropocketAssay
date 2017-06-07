
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

import cornea.geometries.domains
import cornea.geometries.grid
import cornea.geometries.networks
import cornea.geometries.sample_points
import cornea.geometries.pdes
import cornea.parameters.default_parameters
      

class BaseSimulation():
    
    def __init__(self, parameter_collection, work_dir):
        
        self.parameter_collection = parameter_collection
        self.work_dir = work_dir
        self.domain_type = parameter_collection.get_parameter("domain type").value
        self.random_seed = int(parameter_collection.get_parameter("random seed").value)
        self.run_number = parameter_collection.get_parameter("run number").value
        self.reference_length = 1.e-6 * metre()
        self.reference_time = 3600.0 * second()
        self.reference_concentration = 1.e-6*mole_per_metre_cubed()
        self.output_density_quantities = ["Line", "Branch", "Tip"]
        self.domain = None
        self.grid = None
        self.network = None
        self.pde_solver = None
        self.angiogenesis_solver = None
        self.sample_lines = None
        self.denisty_map_grid = None
        self.num_samples_y = 1
        self.num_samples_z = 1
        
    def initialize_reference_scales(self):
        
        BaseUnits.Instance().SetReferenceLengthScale(self.reference_length)
        BaseUnits.Instance().SetReferenceTimeScale(self.reference_time)
        BaseUnits.Instance().SetReferenceConcentrationScale(self.reference_concentration)
        
    def set_up_domain(self):
        
        domain, _ = cornea.geometries.domains.get_domain(self.domain_type, 
                                                        self.parameter_collection, 
                                                        self.reference_length)
        self.domain = domain
        
    def set_up_grid(self):
        
        self.grid = cornea.geometries.grid.get_grid(self.domain_type, self.domain, self.parameter_collection)
        
    def set_up_vessel_network(self):
        
        self.network = cornea.geometries.networks.get_vessel_network(self.domain_type, 
                                                                     self.parameter_collection, 
                                                                     self.reference_length)
        
    def set_up_pde_solver(self, file_handler):
        
        if not self.parameter_collection.get_parameter("use fixed gradient"):
            pde = cornea.geometries.pdes.get_transient_pde(self.domain_type, self.parameter_collection, 
                                            self.reference_length, self.reference_concentration)
            
            self.pde_solver = cornea.geometries.pdes.get_pde_solver(self.domain_type, pde ,self.grid, 
                                                           self.domain, self.parameter_collection, 
                                                           self.reference_length, self.reference_concentration)
            
            pde_time_increment = self.parameter_collection.get_parameter("pde time increment").value
            self.pde_solver.SetFileHandler(file_handler)
            self.pde_solver.SetWriteSolution(True)
            self.pde_solver.SetTargetTimeIncrement(pde_time_increment) # hours
            self.pde_solver.SetUseCoupling(True)
            self.pde_solver.SetWriteIntermediateSolutions(False, int(2.0/pde_time_increment)) # every 2 hours  
            pde.GetDiscreteSources()[0].SetDensityMap(self.pde_solver.GetDensityMap())
            
            if self.network is not None:
                self.pde_solver.GetDensityMap().SetVesselNetwork(self.network) 
        else:
            self.pde_solver = cornea.geometries.pdes.get_pde_fixed_gradient(self.domain_type, self.grid, 
                                                                            self.parameter_collection, 
                                                                            self.reference_length, self.reference_concentration)
            self.pde_solver.SetFileHandler(file_handler)
            self.pde_solver.SetFileName("Vegf_Field")
            self.pde_solver.Write()
                    
    def set_up_angiogenesis_solver(self, file_handler):
        
        if "2" in self.domain_type:
            dimension = str(2)
        else:
            dimension = str(3)
            
        attraction_strength = self.parameter_collection.get_parameter("attraction strength").value  
        chemotactic_strength = self.parameter_collection.get_parameter("chemotactic strength").value  
        persistence_angle = self.parameter_collection.get_parameter("persistence angle").value  
        sprouting_probability = self.parameter_collection.get_parameter("sprouting probability").value  
        tip_exclusion_radius = self.parameter_collection.get_parameter("tip exclusion radius").value 
        do_anastamosis = self.parameter_collection.get_parameter("do anastomosis").value  
        
        migration_rule = getattr(microvessel_chaste.simulation, 'OffLatticeMigrationRule' + dimension)()
        migration_rule.SetDiscreteContinuumSolver(self.pde_solver)
        migration_rule.SetNetwork(self.network)
        migration_rule.SetAttractionStrength(attraction_strength)
        migration_rule.SetChemotacticStrength(chemotactic_strength)
        migration_rule.SetPersistenceAngleSdv((persistence_angle/180.0)*np.pi)
        
        sprouting_rule = getattr(microvessel_chaste.simulation, 'OffLatticeSproutingRule' + dimension)()
        sprouting_rule.SetDiscreteContinuumSolver(self.pde_solver)
        sprouting_rule.SetVesselNetwork(self.network)
        sprouting_rule.SetSproutingProbability(sprouting_probability)
        #sprouting_rule.SetOnlySproutIfPerfused(True)
        sprouting_rule.SetTipExclusionRadius(tip_exclusion_radius)
        
        self.angiogenesis_solver = getattr(microvessel_chaste.simulation, 'AngiogenesisSolver' + dimension)()
        self.angiogenesis_solver.SetVesselNetwork(self.network)
        self.angiogenesis_solver.SetMigrationRule(migration_rule)
        self.angiogenesis_solver.SetSproutingRule(sprouting_rule)
        self.angiogenesis_solver.SetOutputFileHandler(file_handler)
        self.angiogenesis_solver.SetBoundingDomain(self.domain)
        self.angiogenesis_solver.SetDoAnastomosis(do_anastamosis)
        
    def get_sample_points(self):  
        
        self.sample_lines, self.num_samples_y = cornea.geometries.sample_points.get_sample_points(self.domain_type, 
                                                                        self.parameter_collection, 
                                                                        self.reference_length)   
        self.num_samples_z = int(len(self.sample_lines)/self.num_samples_y)  
        
    def set_up_sampling_grid(self):
        
        density_grid_spacing = self.parameter_collection.get_parameter("density grid spacing").value  
        if self.domain_type == "Planar 2D":
            self.denisty_map_grid = microvessel_chaste.mesh.RegularGrid2()
            self.denisty_map_grid.GenerateFromPart(self.domain, density_grid_spacing)
        elif self.domain_type == "Planar 3D":
            self.denisty_map_grid = microvessel_chaste.mesh.RegularGrid3()
            self.denisty_map_grid.GenerateFromPart(self.domain, density_grid_spacing)
        else:
            self.denisty_map_grid = self.grid
        
    def do_sampling(self, output_file, solver, time):  
            
        output_file.write(str(time) + ",")
        for idx in range(self.num_samples_y):
            mean = 0.0
            for jdx in range(self.num_samples_z):
                sample_index = jdx*self.num_samples_y + idx
                solution = solver.GetSolutionP(self.sample_lines[sample_index])
                mean += np.mean(np.array(solution))
            mean /= float(self.num_samples_z)
            output_file.write(str(mean) + ",")              
        output_file.write("\n")
        output_file.flush()
                
    def run(self):
        
        # Initialize the simulation
        chaste.cell_based.SimulationTime.Instance().SetStartTime(0.0)
        chaste.core.RandomNumberGenerator.Instance().Reseed(self.random_seed)
        file_handler = chaste.core.OutputFileHandler(self.work_dir, True)
        pde_only = self.parameter_collection.get_parameter("use pde only").value        
        self.parameter_collection.save(file_handler.GetOutputDirectoryFullPath() + "/adopted_parameter_collection.p")
        
        print "Running Simulation in: ", file_handler.GetOutputDirectoryFullPath()
        print "With Domain Type: ", self.domain_type
        print "With Fixed Gradient: ", self.parameter_collection.get_parameter("use fixed gradient").value
        print "With PDE Only: ", pde_only
        
        self.initialize_reference_scales()
        self.set_up_domain()
        self.set_up_grid()
        
        if not pde_only:
            self.set_up_vessel_network()
        self.set_up_pde_solver(file_handler)
        
        if not pde_only:
            self.set_up_angiogenesis_solver(file_handler)
        
        if "2" in self.domain_type:
            dimension = str(2)
        else:
            dimension = str(3)
        time_step = self.parameter_collection.get_parameter("time step").value 
        total_time = self.parameter_collection.get_parameter("total time").value 
        num_steps = int(total_time/time_step)
        chaste.cell_based.SimulationTime.Instance().SetEndTimeAndNumberOfTimeSteps(total_time, num_steps)
        
        if not pde_only:
            network_writer = getattr(microvessel_chaste.population.vessel, 'VesselNetworkWriter' + dimension)()
        
        self.get_sample_points()
        self.set_up_sampling_grid()

        if not pde_only:
            output_density_files = []
            for eachQuantity in self.output_density_quantities:
                output_file = open(file_handler.GetOutputDirectoryFullPath()+"Sampled_" + 
                                   eachQuantity + "_density.txt", "w")
                output_file.write("Time, ")
                for eachLine in self.sample_lines[0:self.num_samples_y]:
                    y_coord = eachLine.GetPoint(0)[1]
                    output_file.write(str(y_coord)+",")        
                output_file.write("\n")
                output_density_files.append(output_file)
            
        # PDE Sampling
        pde_output_file = open(file_handler.GetOutputDirectoryFullPath()+"Sampled_PDE.txt", "w")
        pde_output_file.write("Time, ")
        for eachLine in self.sample_lines[0:self.num_samples_y]:
            y_coord = eachLine.GetPoint(0)[1]
            pde_output_file.write(str(y_coord)+",")        
        pde_output_file.write("\n") 
            
        # Do the solve
        old_time = chaste.cell_based.SimulationTime.Instance().GetTime()
        while not chaste.cell_based.SimulationTime.Instance().IsFinished():

            print "Simulation Time: ", old_time, " of ", total_time
            elapsed_time = chaste.cell_based.SimulationTime.Instance().GetTimeStepsElapsed()
            time = chaste.cell_based.SimulationTime.Instance().GetTime()
            
            if not self.parameter_collection.get_parameter("use fixed gradient"):
                pde_time_increment = self.parameter_collection.get_parameter("pde time increment").value  
                self.pde_solver.SetStartTime(old_time)
                if time==old_time:
                    self.pde_solver.SetEndTime(time+pde_time_increment)
                else:
                    self.pde_solver.SetEndTime(time)
                self.pde_solver.SetFileName("Vegf_Solution" + str(elapsed_time))
                try:
                    self.pde_solver.Solve()
                except chaste.ChasteException as e:
                    print e.GetMessage
                
            # Sample the pde
            self.do_sampling(pde_output_file, self.pde_solver, time)
            
            if not pde_only:
                density_map = getattr(microvessel_chaste.pde, 'DensityMap' + dimension)()
                density_map.SetGrid(self.denisty_map_grid)
                density_map.SetVesselNetwork(self.network)
        
                density_map_result = getattr(microvessel_chaste.pde, 'FunctionMap' + dimension)()
                density_map_result.SetGrid(self.denisty_map_grid)
                density_map_result.SetVesselNetwork(self.network)
                density_map_result.SetFileHandler(file_handler)
                
                for idx, eachQuantity in enumerate(self.output_density_quantities):
                    density_map_result.SetFileName("/" + eachQuantity + "_Density" + str(elapsed_time))
                    calculation_method = getattr(density_map, 'rGetVessel' + eachQuantity + 'Density')
                    if "Planar" in self.domain_type:
                        density_map_result.UpdateSolution(calculation_method())
                    else:
                        density_map_result.UpdateElementSolution(calculation_method())
                    density_map_result.Write()
    
                    self.do_sampling(output_density_files[idx], density_map_result, time)
        
                network_writer.SetFileName(
                        file_handler.GetOutputDirectoryFullPath() + "/vessel_network_" + str(elapsed_time) + ".vtp")
                network_writer.SetVesselNetwork(self.network)
                network_writer.Write()
        
                # Increment the solver and simulation time
                self.angiogenesis_solver.Increment()
            chaste.cell_based.SimulationTime.Instance().IncrementTimeOneStep()
            old_time = time
            
        if not pde_only:                
            for eachFile in output_density_files:    
                eachFile.close()
        pde_output_file.close()
        chaste.cell_based.SimulationTime.Instance().Destroy()
        
if __name__ == '__main__':
    
    work_dir = "Python/Cornea/TestSimulationBaseDefaults/"
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    domain_types = ["Planar 2D", "Circle 2D", "Planar 3D", "Circle 3D", "Hemisphere 3D"]
    domain_types = ["Planar 2D"]
    random_seeds = [1234]
    
    for eachDomainType in domain_types:
        run_number = 0
        parameter_collection.get_parameter("domain type").value = eachDomainType
        for eachSeed in random_seeds:
            parameter_collection.get_parameter("run number").value = run_number
            parameter_collection.get_parameter("random seed").value = eachSeed
            local_work_dir = work_dir + "/DomainType_" + eachDomainType.replace(" ", "") + "/Run_" + str(run_number)
            simulation = BaseSimulation(parameter_collection, work_dir)
            simulation.run()
            run_number += 1
    