
"""Copyright (c) 2005-2016, University of Oxford.
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
import chaste.cell_based # Chaste Cell Populations
chaste.init() # Initialize MPI and PETSc
import microvessel_chaste.geometry # Geometry tools
import microvessel_chaste.mesh # Meshing
import microvessel_chaste.population.vessel # Vessel tools
import microvessel_chaste.pde # PDE and solvers
import microvessel_chaste.simulation # Flow and angiogenesis solvers
import microvessel_chaste.visualization # Visualization
from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience

           
def do_2d_run(work_dir, num_runs, rand_seed, parameter_values, sensing = False, moore = False):
    
    chaste.core.RandomNumberGenerator.Instance().Reseed(rand_seed)
    
    for run_index in range(num_runs):
        chaste.cell_based.SimulationTime.Instance().SetStartTime(0.0)
        
        file_handler = chaste.core.OutputFileHandler(work_dir+"/Run_"+str(run_index)+"/", False)

        reference_length = 1.e-6 * metre()
        reference_time = 3600.0 * second()
        BaseUnits.Instance().SetReferenceLengthScale(reference_length)
        BaseUnits.Instance().SetReferenceTimeScale(reference_time)
        BaseUnits.Instance().SetReferenceConcentrationScale(1.e-6*mole_per_metre_cubed())
        
        domain = microvessel_chaste.geometry.Part2()
        domain.AddRectangle(2.e-3*metre(), 1.24e-3*metre(), 
                            microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0))
        grid = microvessel_chaste.mesh.RegularGrid2()
        grid.GenerateFromPart(domain, parameter_values["Grid Spacing"])
    
        # network
        nodes = []
        for idx in range(100):
            nodes.append(microvessel_chaste.population.vessel.VesselNode2.Create(float(idx)*20.0, 200.0, 0.0, 
                                                                                 reference_length))
        nodes[0].GetFlowProperties().SetIsInputNode(True)
        nodes[0].GetFlowProperties().SetPressure(Owen11Parameters.mpInletPressure.GetValue("User"))
        nodes[-1].GetFlowProperties().SetIsOutputNode(True)
        nodes[-1].GetFlowProperties().SetPressure(Owen11Parameters.mpOutletPressure.GetValue("User"))
        vessel = microvessel_chaste.population.Vessel2.Create(nodes)
        network = microvessel_chaste.population.VesselNetwork2.Create()
        network.AddVessel(vessel)
    
        # pde
        vegf_pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde2_2()
        vegf_pde.SetIsotropicDiffusionConstant(parameter_values["VEGF Diffusivity"])
        vegf_pde.SetContinuumLinearInUTerm(parameter_values["VEGF Decay Rate"])
        vegf_pde.SetMultiplierValue(parameter_values["Initial VEGF Concentration"])
        vessel_vegf_sink = microvessel_chaste.pde.VesselBasedDiscreteSource2()
        vessel_vegf_sink.SetReferenceConcentration(0.0*mole_per_metre_cubed())
        vessel_vegf_sink.SetVesselPermeability(parameter_values["VEGF Permeability"])
        vessel_vegf_sink.SetReferenceHaematocrit(0.45*dimensionless())
        vegf_pde.SetPelletBindingConstant(30000.0*dimensionless())
        #vessel_vegf_sink.SetUptakeRatePerCell(-(4.e-22/3600.0)*mole_per_second())
        vegf_pde.AddDiscreteSource(vessel_vegf_sink)
        vegf_solver = microvessel_chaste.pde.FiniteDifferenceSolver2()
        vegf_solver.SetParabolicPde(vegf_pde)
        vegf_solver.SetLabel("vegf")
        vegf_solver.SetGrid(grid) 
        
        ## Flow
        large_vessel_radius = 5.0e-6 * metre()
        network.SetSegmentRadii(large_vessel_radius)
        viscosity = Owen11Parameters.mpPlasmaViscosity.GetValue("User")
        network.SetSegmentViscosity(viscosity);
        
        ## Set up the pre- and post flow calculators.
        impedance_calculator = microvessel_chaste.simulation.VesselImpedanceCalculator2()
        haematocrit_calculator = microvessel_chaste.simulation.ConstantHaematocritSolver2()
        haematocrit_calculator.SetHaematocrit(Owen11Parameters.mpInflowHaematocrit.GetValue("User"))
        wss_calculator = microvessel_chaste.simulation.WallShearStressCalculator2()
        viscosity_calculator = microvessel_chaste.simulation.ViscosityCalculator2()
    
        ## Set up and configure the structural adaptation solver.
        structural_adaptation_solver = microvessel_chaste.simulation.StructuralAdaptationSolver2()
        structural_adaptation_solver.SetTolerance(0.0001)
        structural_adaptation_solver.SetMaxIterations(2)
        structural_adaptation_solver.SetTimeIncrement(Owen11Parameters.mpVesselRadiusUpdateTimestep.GetValue("User"));
        structural_adaptation_solver.AddPreFlowSolveCalculator(impedance_calculator)
        structural_adaptation_solver.AddPostFlowSolveCalculator(haematocrit_calculator)
        structural_adaptation_solver.AddPostFlowSolveCalculator(wss_calculator)
        structural_adaptation_solver.AddPostFlowSolveCalculator(viscosity_calculator)
        
        ## Set up a regression solver.
        regression_solver = microvessel_chaste.simulation.WallShearStressBasedRegressionSolver2()
        
        ## Set up an angiogenesis solver and add sprouting and migration rules.
        
        angiogenesis_solver = microvessel_chaste.simulation.AngiogenesisSolver2()
        sprouting_rule = microvessel_chaste.simulation.Owen2011SproutingRule2()
        migration_rule = microvessel_chaste.simulation.TipAttractionLatticeBasedMigrationRule2()
        sprouting_rule.SetDiscreteContinuumSolver(vegf_solver)
        migration_rule.SetDiscreteContinuumSolver(vegf_solver)
        
        sprouting_rule.SetTipExclusionRadius(60.e-6*metre())
        sprouting_rule.SetVesselEndCutoff(0.0*metre())
        #sprouting_rule.SetSproutingProbability(7.5e-4*per_second())
        sprouting_rule.SetSproutingProbability(.75e-3*per_second())
        #sprouting_rule.SetHalfMaxVegf(0.5e-9*mole_per_metre_cubed())
        sprouting_rule.SetHalfMaxVegf(0.5e-9*mole_per_metre_cubed())
        migration_rule.SetDiscreteContinuumSolver(vegf_solver)
        migration_rule.SetUseMooreNeighbourhood(moore)
        migration_rule.SetUseTipAttraction(sensing)
        migration_rule.SetCellChemotacticParameter((5.7/3600.0)*metre_pow5_per_second_per_mole())
        
        angiogenesis_solver.SetMigrationRule(migration_rule)
        angiogenesis_solver.SetSproutingRule(sprouting_rule)
        angiogenesis_solver.SetVesselGrid(grid)
        angiogenesis_solver.SetVesselNetwork(network)

        ## The microvessel solver will manage all aspects of the vessel solve.
        microvessel_solver = microvessel_chaste.simulation.MicrovesselSolver2()
        microvessel_solver.SetVesselNetwork(network)
        microvessel_solver.SetOutputFrequency(6)
        microvessel_solver.SetOutputFileHandler(file_handler)
        microvessel_solver.AddDiscreteContinuumSolver(vegf_solver)
        microvessel_solver.SetStructuralAdaptationSolver(structural_adaptation_solver)
        microvessel_solver.SetRegressionSolver(regression_solver)
        microvessel_solver.SetAngiogenesisSolver(angiogenesis_solver)
        
        chaste.cell_based.SimulationTime.Instance().SetEndTimeAndNumberOfTimeSteps(parameter_values["End Time"], 
                                                                                   parameter_values["Time Steps"])
        microvessel_solver.Run()
        
        ## Dump the parameters to file for inspection.
        ParameterCollection.Instance().DumpToFile(file_handler.GetOutputDirectoryFullPath()+"parameter_collection.xml")

        network_recorder.line_density_file.close() 
        network_recorder.branch_density_file.close() 
        network_recorder.perfusion_density_file.close() 
        
        chaste.cell_based.SimulationTime.Instance().Destroy()
        
def do_3d_run(work_dir, num_runs, rand_seed, parameter_values, sensing = False, moore = False):
    
    chaste.core.RandomNumberGenerator.Instance().Reseed(rand_seed)
    
    for run_index in range(num_runs):
        chaste.cell_based.SimulationTime.Instance().SetStartTime(0.0)
        
        file_handler = chaste.core.OutputFileHandler(work_dir+"/Run_"+str(run_index)+"/", False)

        reference_length = 1.e-6 * metre()
        reference_time = 3600.0 * second()
        BaseUnits.Instance().SetReferenceLengthScale(reference_length)
        BaseUnits.Instance().SetReferenceTimeScale(reference_time)
        BaseUnits.Instance().SetReferenceConcentrationScale(1.e-6*mole_per_metre_cubed())
        
        domain = microvessel_chaste.geometry.Part3()
        domain.AddCuboid(2.e-3*metre(), 1.24e-3*metre(), 100.0e-6*metre(),
                            microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0))
        grid = microvessel_chaste.mesh.RegularGrid3()
        grid.GenerateFromPart(domain, parameter_values["Grid Spacing"])
    
        # network
        nodes = []
        for idx in range(100):
            nodes.append(microvessel_chaste.population.vessel.VesselNode3.Create(float(idx)*20.0, 200.0, 40.0, 
                                                                                 reference_length))
        nodes[0].GetFlowProperties().SetIsInputNode(True)
        nodes[0].GetFlowProperties().SetPressure(Owen11Parameters.mpInletPressure.GetValue("User"))
        nodes[-1].GetFlowProperties().SetIsOutputNode(True)
        nodes[-1].GetFlowProperties().SetPressure(Owen11Parameters.mpOutletPressure.GetValue("User"))
        vessel = microvessel_chaste.population.Vessel3.Create(nodes)
        network = microvessel_chaste.population.VesselNetwork3.Create()
        network.AddVessel(vessel)
    
        # pde
        vegf_pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde3_3()
        vegf_pde.SetIsotropicDiffusionConstant(parameter_values["VEGF Diffusivity"])
        vegf_pde.SetContinuumLinearInUTerm(parameter_values["VEGF Decay Rate"])
        vegf_pde.SetMultiplierValue(parameter_values["Initial VEGF Concentration"])
        vegf_pde.SetPelletBindingConstant(30000.0*dimensionless())
        vessel_vegf_sink = microvessel_chaste.pde.VesselBasedDiscreteSource3()
        vessel_vegf_sink.SetReferenceConcentration(0.0*mole_per_metre_cubed())
        vessel_vegf_sink.SetVesselPermeability(parameter_values["VEGF Permeability"])
        #vessel_vegf_sink.SetUptakeRatePerCell(-(4.e-22/3600.0)*mole_per_second())
        vessel_vegf_sink.SetReferenceHaematocrit(0.45*dimensionless())
        vegf_pde.AddDiscreteSource(vessel_vegf_sink)
        vegf_solver = microvessel_chaste.pde.FiniteDifferenceSolver3()
        vegf_solver.SetParabolicPde(vegf_pde)
        vegf_solver.SetLabel("vegf")
        vegf_solver.SetGrid(grid) 
        
        ## Flow
        large_vessel_radius = 5.0e-6 * metre()
        network.SetSegmentRadii(large_vessel_radius)
        viscosity = Owen11Parameters.mpPlasmaViscosity.GetValue("User")
        network.SetSegmentViscosity(viscosity);
        
        ## Set up the pre- and post flow calculators.
        impedance_calculator = microvessel_chaste.simulation.VesselImpedanceCalculator3()
        haematocrit_calculator = microvessel_chaste.simulation.ConstantHaematocritSolver3()
        haematocrit_calculator.SetHaematocrit(Owen11Parameters.mpInflowHaematocrit.GetValue("User"))
        wss_calculator = microvessel_chaste.simulation.WallShearStressCalculator3()
        viscosity_calculator = microvessel_chaste.simulation.ViscosityCalculator3()
    
        ## Set up and configure the structural adaptation solver.
        structural_adaptation_solver = microvessel_chaste.simulation.StructuralAdaptationSolver3()
        structural_adaptation_solver.SetTolerance(0.0001)
        structural_adaptation_solver.SetMaxIterations(2)
        structural_adaptation_solver.SetTimeIncrement(Owen11Parameters.mpVesselRadiusUpdateTimestep.GetValue("User"));
        structural_adaptation_solver.AddPreFlowSolveCalculator(impedance_calculator)
        structural_adaptation_solver.AddPostFlowSolveCalculator(haematocrit_calculator)
        structural_adaptation_solver.AddPostFlowSolveCalculator(wss_calculator)
        structural_adaptation_solver.AddPostFlowSolveCalculator(viscosity_calculator)
        
        ## Set up a regression solver.
        regression_solver = microvessel_chaste.simulation.WallShearStressBasedRegressionSolver3()
        
        ## Set up an angiogenesis solver and add sprouting and migration rules.
        
        angiogenesis_solver = microvessel_chaste.simulation.AngiogenesisSolver3()
        sprouting_rule = microvessel_chaste.simulation.Owen2011SproutingRule3()
        migration_rule = microvessel_chaste.simulation.TipAttractionLatticeBasedMigrationRule3()
        sprouting_rule.SetDiscreteContinuumSolver(vegf_solver)
        migration_rule.SetDiscreteContinuumSolver(vegf_solver)
        
        sprouting_rule.SetTipExclusionRadius(60.e-6*metre())
        sprouting_rule.SetVesselEndCutoff(0.0*metre())
        #sprouting_rule.SetSproutingProbability(7.5e-4*per_second()) # original run
        #sprouting_rule.SetSproutingProbability(4.0e-4*per_second()) # new run
        sprouting_rule.SetSproutingProbability(0.75e-3*per_second()) # lowest run
        #sprouting_rule.SetHalfMaxVegf(0.5e-9*mole_per_metre_cubed())
        sprouting_rule.SetHalfMaxVegf(0.65e-9*mole_per_metre_cubed())
        migration_rule.SetDiscreteContinuumSolver(vegf_solver)
        migration_rule.SetUseMooreNeighbourhood(moore)
        migration_rule.SetUseTipAttraction(sensing)
        migration_rule.SetCellChemotacticParameter((5.7/3600.0)*metre_pow5_per_second_per_mole())
        
        angiogenesis_solver.SetMigrationRule(migration_rule)
        angiogenesis_solver.SetSproutingRule(sprouting_rule)
        angiogenesis_solver.SetVesselGrid(grid)
        angiogenesis_solver.SetVesselNetwork(network)

        ## The microvessel solver will manage all aspects of the vessel solve.
        microvessel_solver = microvessel_chaste.simulation.MicrovesselSolver3()
        microvessel_solver.SetVesselNetwork(network)
        microvessel_solver.SetOutputFrequency(6)
        microvessel_solver.SetOutputFileHandler(file_handler)
        microvessel_solver.AddDiscreteContinuumSolver(vegf_solver)
        microvessel_solver.SetStructuralAdaptationSolver(structural_adaptation_solver)
        microvessel_solver.SetRegressionSolver(regression_solver)
        microvessel_solver.SetAngiogenesisSolver(angiogenesis_solver)
        
        chaste.cell_based.SimulationTime.Instance().SetEndTimeAndNumberOfTimeSteps(parameter_values["End Time"], 
                                                                                   parameter_values["Time Steps"])
        microvessel_solver.Run()
        
        ## Dump the parameters to file for inspection.
        ParameterCollection.Instance().DumpToFile(file_handler.GetOutputDirectoryFullPath()+"parameter_collection.xml")

        network_recorder.line_density_file.close() 
        network_recorder.branch_density_file.close() 
        network_recorder.perfusion_density_file.close() 
        
        chaste.cell_based.SimulationTime.Instance().Destroy()
        
if __name__ == '__main__':
    
    work_dir = "Python/TestTipSensingPaper_ThesisParams_Jan25"
    num_runs = 10
    rand_seed = 1234
    
    parameter_values = {"VEGF Diffusivity":6.94e-11 * metre_squared_per_second(),
                        "VEGF Decay Rate":(-0.8/3600.0) * per_second(),
                        "Initial VEGF Concentration":3.93e-1*mole_per_metre_cubed(),
                        "VEGF Permeability": (3.e-4/3600.0)*metre_per_second(),
                        "End Time": 24,#120.0,
                        "Time Steps": 144,#720,
                        "Grid Spacing": 20.0e-6*metre()}
    
    # 2D, von neumann, no sensing
    #do_2d_run(work_dir+"/2D_VN_NoSense", num_runs, rand_seed, parameter_values, sensing = False, moore = False)
    
    # 2D, von neumann, sensing
    #do_2d_run(work_dir+"/2D_VN_Sense", num_runs, rand_seed, parameter_values, sensing = True, moore = False)
    
    # 2D, moore, sensing
    #do_2d_run(work_dir+"/2D_Moore_Sense", num_runs, rand_seed, parameter_values, sensing = True, moore = True)
    
    # 2D, moore, no sensing
    do_2d_run(work_dir+"/2D_Moore_NoSense", num_runs, rand_seed, parameter_values, sensing = False, moore = True)
    
    # 3D, von neumann, no sensing
    #do_3d_run(work_dir+"/3D_VN_NoSense", num_runs, rand_seed, parameter_values, sensing = False, moore = False)
    
    # 3D, von neumann, sensing
    #do_3d_run(work_dir+"/3D_VN_Sense", num_runs, rand_seed, parameter_values, sensing = True, moore = False)
    
    # 3D, moore, sensing
    #do_3d_run(work_dir+"/3D_Moore_Sense", num_runs, rand_seed, parameter_values, sensing = True, moore = True)
    
    #rand_seed = 5678
    #do_3d_run(work_dir+"/3D_Moore_Sense_b", num_runs, rand_seed, parameter_values, sensing = True, moore = True)
    
    # 3D, moore, no sensing
    #do_3d_run(work_dir+"/3D_Moore_NoSense", num_runs, rand_seed, parameter_values, sensing = False, moore = True)

    #rand_seed = 5678
    #do_3d_run(work_dir+"/3D_Moore_NoSense_b", num_runs, rand_seed, parameter_values, sensing = False, moore = True)