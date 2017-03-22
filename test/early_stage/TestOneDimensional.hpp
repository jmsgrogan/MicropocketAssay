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
#include <boost/lexical_cast.hpp>
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

#include "PetscAndVtkSetupAndFinalize.hpp"

class TestOneDimensionalDomain : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestTransportOnly() throw(Exception)
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
};

#endif /*TESTCOUPLEDLUMPEDSYSTEMFINITEDIFFERENCESOLVER_HPP_*/
