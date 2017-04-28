/*

Copyright (c) 2005-2017, University of Oxford.
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

*/

#ifndef TESTONEDIMENSIONALDOMAIN_HPP_
#define TESTONEDIMENSIONALDOMAIN_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Part.hpp"
#include "UnitCollection.hpp"
#include "CoupledVegfPelletDiffusionReactionPde.hpp"
#include "CoupledLumpedSystemFiniteDifferenceSolver.hpp"
#include "VesselNetwork.hpp"
#include "VesselNetworkGenerator.hpp"
#include "DiscreteContinuumBoundaryCondition.hpp"
#include "DiscreteContinuumMesh.hpp"
#include "DiscreteContinuumMeshGenerator.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CoupledLumpedSystemFiniteElementSolver.hpp"
#include "SetUpDomainsAndPdes.hpp"
#include "FunctionMap.hpp"
#include "OffLatticeMigrationRule.hpp"
#include "OffLatticeSproutingRule.hpp"
#include "AngiogenesisSolver.hpp"
#include "DensityMap.hpp"
#include "VesselNetworkWriter.hpp"

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestOneDimensionalDomain : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void xTestTransportOnly() throw(Exception)
    {
        InitializeReferenceScales();
        units::quantity<unit::length> pellet_height = 1.0e-3*unit::metres;
        units::quantity<unit::length> cornea_radius = 1.3e-3*unit::metres;

        boost::shared_ptr<Part<2> > p_domain = Get2DPlaneDomain(pellet_height, cornea_radius);
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        p_grid->GenerateFromPart(p_domain, 40.0e-6*unit::metres);
        boost::shared_ptr<CoupledVegfPelletDiffusionReactionPde<2> > p_pde = GetPde<2>();
        p_pde->SetPelletSurfaceArea(2.0*cornea_radius*M_PI*1e-6*unit::metres);

        // Solve the finite difference problem
        CoupledLumpedSystemFiniteDifferenceSolver<2> fd_solver;
        fd_solver.SetGrid(p_grid);
        fd_solver.SetPde(p_pde);
        MAKE_PTR_ARGS(OutputFileHandler, p_fd_output_file_handler, ("TestOneDimensionalDomain/TestTransportOnly"));
        fd_solver.SetFileHandler(p_fd_output_file_handler);
        fd_solver.SetWriteSolution(true);
        fd_solver.SetTargetTimeIncrement(0.1); // hours
        fd_solver.SetStartTime(0.0);
        fd_solver.SetUseCoupling(true);
        fd_solver.SetEndTime(24.0*6); // 6 days
        fd_solver.SetWriteIntermediateSolutions(true, 10*24); // every 2 hours
        fd_solver.Solve();

        // Test the intermediate solutions
        vtkSmartPointer<vtkPoints> p_sample_points = GetSamplePointsFor2DPlaneDomain(pellet_height, cornea_radius);
        std::vector<std::pair<std::vector<double>, double> > intermediate_solutions =
                fd_solver.rGetIntermediateSolutions();

        std::ofstream outfile;
        std::string file_name = p_fd_output_file_handler->GetOutputDirectoryFullPath()+"sample_values.txt";
        outfile.open(file_name.c_str());

        outfile << "Time, ";
        for(unsigned idx=0; idx<p_sample_points->GetNumberOfPoints(); idx++)
        {
            double x_loc = p_sample_points->GetPoint(idx)[0];
            outfile << boost::lexical_cast<std::string>(x_loc) << ",";
        }
        outfile << "\n";
        for(unsigned idx=0; idx<intermediate_solutions.size();idx++)
        {
            double time = intermediate_solutions[idx].second;
            outfile << boost::lexical_cast<std::string>(time) << ",";

            fd_solver.UpdateSolution(intermediate_solutions[idx].first);
            std::vector<units::quantity<unit::concentration> > solution = fd_solver.GetConcentrations(p_sample_points);
            for(unsigned jdx=0; jdx<solution.size(); jdx++)
            {
                double c_numerical_nondim = solution[jdx]/(1.e-9*unit::mole_per_metre_cubed);
                outfile << boost::lexical_cast<std::string>(c_numerical_nondim) << ",";
            }
            outfile << "\n";
        }
        outfile.close();
    }

    void TestVesselsOnly() throw(Exception)
    {
        MAKE_PTR_ARGS(OutputFileHandler, p_output_file_handler, ("TestOneDimensionalDomain/TestVesselOnly"));
        InitializeReferenceScales();
        RandomNumberGenerator::Instance()->Reseed(12345);

        units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();
        units::quantity<unit::length> pellet_height = 1.2e-3*unit::metres;
        units::quantity<unit::length> cornea_radius = (1.0e-3)*unit::metres;

        boost::shared_ptr<Part<2> > p_domain = Get2DPlaneDomain(pellet_height, cornea_radius);
        boost::shared_ptr<RegularGrid<2> > p_grid = RegularGrid<2>::Create();
        units::quantity<unit::length> spacing(40.0*unit::microns); //um
        p_grid->GenerateFromPart(p_domain, spacing);

        // Prescribe a linearly increasing vegf field using a function map
        boost::shared_ptr<FunctionMap<2> > p_funciton_map = FunctionMap<2>::Create();
        p_funciton_map->SetGrid(p_grid);
        std::vector<units::quantity<unit::concentration> > vegf_field = std::vector<units::quantity<unit::concentration> >(p_grid->GetNumberOfPoints(), 0.0*unit::mole_per_metre_cubed);
        for (unsigned idx = 0; idx < p_grid->GetNumberOfPoints(); idx++)
        {
            vegf_field[idx] = 0.3*(p_grid->GetPoint(idx).GetLocation(reference_length)[1] / (pellet_height/reference_length))*1.e-9*unit::mole_per_metre_cubed;
        }

        p_funciton_map->SetFileHandler(p_output_file_handler);
        p_funciton_map->SetFileName("Function");
        p_funciton_map->SetLabel("vegf");
        p_funciton_map->UpdateSolution(vegf_field);
        p_funciton_map->Write();

        //Set up the limbal vessel
        VesselNetworkGenerator<2> generator;
        unsigned divisions = unsigned(2.0*M_PI*cornea_radius/spacing) - 2; // divide the vessel to coincide with grid
        unsigned alignment_axis = 0; // pointing x direction
        boost::shared_ptr<VesselNetwork<2> > p_network = generator.GenerateSingleVessel(cornea_radius*2.0*M_PI,
                                                                                        DimensionalChastePoint<2>(0.0, 5.0, 0.0, spacing),
                                                                                            divisions, alignment_axis);

        for(unsigned idx=0;idx<p_network->GetNodes().size();idx++)
        {
            p_network->GetNodes()[idx]->GetFlowProperties()->SetPressure(1.0*unit::pascals);
        }

        boost::shared_ptr<OffLatticeMigrationRule<2> > p_migration_rule = OffLatticeMigrationRule<2>::Create();
        p_migration_rule->SetDiscreteContinuumSolver(p_funciton_map);
        p_migration_rule->SetNetwork(p_network);
        p_migration_rule->SetAttractionStrength(0.0);
        p_migration_rule->SetChemotacticStrength(0.0);
        p_migration_rule->SetPersistenceAngleSdv((2.0/180.0)*M_PI);

        boost::shared_ptr<OffLatticeSproutingRule<2> > p_sprouting_rule = OffLatticeSproutingRule<2>::Create();
        p_sprouting_rule->SetDiscreteContinuumSolver(p_funciton_map);
        p_sprouting_rule->SetVesselNetwork(p_network);
        p_sprouting_rule->SetSproutingProbability(0.5 /(3600.0*unit::seconds));
        p_sprouting_rule->SetOnlySproutIfPerfused(true);
        p_sprouting_rule->SetTipExclusionRadius(40.0e-6*unit::metres);

        AngiogenesisSolver<2> angiogenesis_solver;
        angiogenesis_solver.SetVesselNetwork(p_network);
        angiogenesis_solver.SetMigrationRule(p_migration_rule);
        angiogenesis_solver.SetSproutingRule(p_sprouting_rule);
        angiogenesis_solver.SetOutputFileHandler(p_output_file_handler);
        angiogenesis_solver.SetBoundingDomain(p_domain);
        angiogenesis_solver.SetDoAnastomosis(true);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(48.0, 96);
        // Set up a vessel network writer
        boost::shared_ptr<VesselNetworkWriter<2> > p_network_writer = VesselNetworkWriter<2>::Create();

        boost::shared_ptr<RegularGrid<2> > p_denisty_map_grid = RegularGrid<2>::Create();
        units::quantity<unit::length> density_grid_spacing(80.0*unit::microns); //um
        p_denisty_map_grid->GenerateFromPart(p_domain, density_grid_spacing);

        // Set up the sample lines and output files
        units::quantity<unit::length> sample_spacing_x = 40.0e-6*unit::metres;
        units::quantity<unit::length> sample_spacing_y = 20.0e-6*unit::metres;
        unsigned num_sample_points_x = 2.0*M_PI*cornea_radius/sample_spacing_x + 1u;
        unsigned num_sample_points_y = (pellet_height-200.0e-6*unit::metres)/sample_spacing_y + 1u;

        std::ofstream outfile;
        std::string file_name = p_output_file_handler->GetOutputDirectoryFullPath()+"sampled_line_density.txt";
        outfile.open(file_name.c_str());
        outfile << "Time, ";

        std::vector<vtkSmartPointer<vtkPoints> > sample_lines;
        for(unsigned idx=0; idx<num_sample_points_y; idx++)
        {
            vtkSmartPointer<vtkPoints> p_sample_points = vtkSmartPointer<vtkPoints>::New();
            for(unsigned jdx=0; jdx<num_sample_points_x; jdx++)
            {
                p_sample_points->InsertNextPoint(double(jdx)*(sample_spacing_x/reference_length),
                        (double(idx)*sample_spacing_y)/reference_length+200.0, 0.0);
            }
            sample_lines.push_back(p_sample_points);
            outfile << boost::lexical_cast<std::string>(double(double(idx)*sample_spacing_y/reference_length)) << ",";
        }
        outfile << "\n";

        // Loop for the duration of the simulation time
        while (!SimulationTime::Instance()->IsFinished())
        {
            boost::shared_ptr<DensityMap<2> > p_density_map = DensityMap<2>::Create();
            p_density_map->SetGrid(p_denisty_map_grid);
            p_density_map->SetVesselNetwork(p_network);

            boost::shared_ptr<FunctionMap<2> > p_density_map_result = FunctionMap<2>::Create();
            p_density_map_result->SetGrid(p_denisty_map_grid);
            p_density_map_result->SetVesselNetwork(p_network);
            p_density_map_result->SetFileHandler(p_output_file_handler);
            p_density_map_result->SetFileName("/line_density"
                    + boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed()));
            p_density_map_result->UpdateSolution(p_density_map->rGetVesselLineDensity());
            p_density_map_result->Write();

            // Sample the density map
            double time = SimulationTime::Instance()->GetTime();
            outfile << boost::lexical_cast<std::string>(time) << ",";
            for(unsigned idx=0;idx<sample_lines.size();idx++)
            {
                std::vector<double> solution = p_density_map_result->GetSolution(sample_lines[idx]);
                double sum = std::accumulate(solution.begin(), solution.end(), 0.0);
                double mean = sum / solution.size();
                outfile << boost::lexical_cast<std::string>(mean) << ",";
            }
            outfile << std::endl;

            p_density_map_result->SetFileName("/tip_density"
                    + boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed()));
            p_density_map_result->UpdateSolution(p_density_map->rGetVesselTipDensity());
            p_density_map_result->Write();

            p_density_map_result->SetFileName("/branch_density"
                    + boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed()));
            p_density_map_result->UpdateSolution(p_density_map->rGetVesselBranchDensity());
            p_density_map_result->Write();

            p_network_writer->SetFileName(
                    p_output_file_handler->GetOutputDirectoryFullPath() + "/vessel_network_"
                            + boost::lexical_cast<std::string>(SimulationTime::Instance()->GetTimeStepsElapsed())
                            + ".vtp");
            p_network_writer->SetVesselNetwork(p_network);
            p_network_writer->Write();

            // Increment the solver and simulation time
            angiogenesis_solver.Increment();
            SimulationTime::Instance()->IncrementTimeOneStep();
        }
        outfile.close();
    }
};

#endif /*TESTCOUPLEDLUMPEDSYSTEMFINITEDIFFERENCESOLVER_HPP_*/
