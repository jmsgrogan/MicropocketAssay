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

#ifndef SETUPDOMAINSANDPDES_HPP_
#define SETUPDOMAINSANDPDES_HPP_

#include <vector>
#include <string>
#include <math.h>
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
#include "CoupledLumpedSystemFiniteElementSolver.hpp"

inline void InitializeReferenceScales()
{
    BaseUnits::Instance()->SetReferenceLengthScale(1.e-6*unit::metres);
    BaseUnits::Instance()->SetReferenceConcentrationScale(1.e-9*unit::mole_per_metre_cubed);
    BaseUnits::Instance()->SetReferenceTimeScale(3600.0*unit::seconds);
}

template <unsigned DIM>
inline boost::shared_ptr<CoupledVegfPelletDiffusionReactionPde<DIM> > GetPde()
{
    boost::shared_ptr<CoupledVegfPelletDiffusionReactionPde<DIM> > p_pde = CoupledVegfPelletDiffusionReactionPde<DIM>::Create();
    units::quantity<unit::diffusivity> vegf_diffusivity(6.94e-11 * unit::metre_squared_per_second);
    units::quantity<unit::rate> vegf_decay_rate((-0.8/3600.0) * unit::per_second);
    p_pde->SetIsotropicDiffusionConstant(vegf_diffusivity);
    p_pde->SetContinuumLinearInUTerm(vegf_decay_rate);
    units::quantity<unit::concentration> initial_vegf_concentration(3.93e-1*unit::mole_per_metre_cubed);
    p_pde->SetCurrentVegfInPellet(initial_vegf_concentration);
    p_pde->SetPelletBindingConstant(100.0);
    p_pde->SetPelletDepth(50.0e-6*unit::metres);
    p_pde->SetPelletSurfaceArea(2000e-6*unit::metres*p_pde->GetPelletDepth());
    p_pde->SetCorneaPelletPermeability(0.002*p_pde->GetCorneaPelletPermeability());
    return p_pde;
}

inline boost::shared_ptr<Part<2> > Get2DCircleDomain(units::quantity<unit::length> pellet_height,
        units::quantity<unit::length> cornea_radius, units::quantity<unit::length> pellet_radius)
{
    units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();
    units::quantity<unit::length> delta = pellet_height-cornea_radius+pellet_radius;
    boost::shared_ptr<Part<2> > p_domain = Part<2> ::Create();
    p_domain->AddCircle(cornea_radius, DimensionalChastePoint<2>(0.0, 0.0, 0.0));
    boost::shared_ptr<Polygon<2> > p_polygon = p_domain->AddCircle(pellet_radius,
            DimensionalChastePoint<2>(0.0, -delta/reference_length, 0.0, reference_length));
    p_polygon->LabelAllEdges("Inner Boundary");
    p_domain->AddHoleMarker(DimensionalChastePoint<2>(0.0, -delta/reference_length, 0.0, reference_length));
    return p_domain;
}

inline boost::shared_ptr<Part<2> > Get3DCircleDomain(units::quantity<unit::length> pellet_height,
        units::quantity<unit::length> cornea_radius, units::quantity<unit::length> pellet_radius)
{
    units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();
    units::quantity<unit::length> delta = pellet_height-cornea_radius+pellet_radius;
    boost::shared_ptr<Part<2> > p_domain = Part<2> ::Create();
    p_domain->AddCircle(cornea_radius, DimensionalChastePoint<2>(0.0, 0.0, 0.0));
    boost::shared_ptr<Polygon<2> > p_polygon = p_domain->AddCircle(pellet_radius,
            DimensionalChastePoint<2>(0.0, -delta/reference_length, 0.0, reference_length));
    p_polygon->LabelAllEdges("Inner Boundary");
    p_domain->AddHoleMarker(DimensionalChastePoint<2>(0.0, -delta/reference_length, 0.0, reference_length));
    return p_domain;
}

inline boost::shared_ptr<Part<2> > Get2DPlaneDomain(units::quantity<unit::length> pellet_height,
        units::quantity<unit::length> cornea_radius)
{
    boost::shared_ptr<Part<2> > p_domain = Part<2>::Create();
    units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();
    p_domain->AddRectangle(2.0*M_PI*cornea_radius, pellet_height, DimensionalChastePoint<2>(0.0, 0.0, 0.0));

    p_domain->LabelEdges(DimensionalChastePoint<2>(M_PI*cornea_radius/reference_length,
            pellet_height/reference_length, 0, reference_length), "Top Boundary");
    return p_domain;
}

inline boost::shared_ptr<Part<2> > Get2DPlaneDomainWithPellet(units::quantity<unit::length> pellet_height,
        units::quantity<unit::length> cornea_radius, units::quantity<unit::length> pellet_radius)
{
    // Manually build the polygon
    units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();
    units::quantity<unit::length> cornea_width = 2.0*M_PI*cornea_radius;
    units::quantity<unit::length> pellet_width = 2.0*M_PI*pellet_radius;

    std::vector<boost::shared_ptr<DimensionalChastePoint<2> > > points;
    points.push_back(DimensionalChastePoint<2>::Create(0.0, 0.0, 0.0, reference_length));
    points.push_back(DimensionalChastePoint<2>::Create(cornea_width/reference_length, 0.0, 0.0, reference_length));
    points.push_back(DimensionalChastePoint<2>::Create(cornea_width/reference_length, pellet_height/reference_length, 0.0,
            reference_length));
    points.push_back(DimensionalChastePoint<2>::Create(cornea_width/(2.0*reference_length) + pellet_width/(2.0*reference_length), pellet_height/reference_length, 0.0,
            reference_length));
    points.push_back(DimensionalChastePoint<2>::Create(cornea_width/(2.0*reference_length) - pellet_width/(2.0*reference_length), pellet_height/reference_length, 0.0,
            reference_length));
    points.push_back(DimensionalChastePoint<2>::Create(0.0, pellet_height/reference_length, 0.0, reference_length));

    boost::shared_ptr<Polygon<2> > p_polygon = Polygon<2>::Create(points);
    boost::shared_ptr<Part<2> > p_domain = Part<2>::Create();
    p_domain->AddPolygon(p_polygon);
    p_domain->LabelEdges(DimensionalChastePoint<2>(cornea_width/(2.0*reference_length),
            pellet_height/reference_length, 0, reference_length), "Top Boundary");
    return p_domain;
}

inline vtkSmartPointer<vtkPoints> GetSamplePointsFor2DPlaneDomain(units::quantity<unit::length> pellet_height,
        units::quantity<unit::length> cornea_radius)
{
    units::quantity<unit::length> cell_spacing = 40.0e-6*unit::metres;
    units::quantity<unit::length> vessel_height = 200.0e-6*unit::metres;

    unsigned num_sample_points = 2.0*M_PI*cornea_radius/cell_spacing + 1u;
    units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();

    vtkSmartPointer<vtkPoints> p_sample_points = vtkSmartPointer<vtkPoints>::New();
    for(unsigned idx=0; idx<num_sample_points; idx++)
    {
        p_sample_points->InsertNextPoint(double(idx)*(cell_spacing/reference_length),
                vessel_height/reference_length, 0.0);
    }
    return p_sample_points;
}

inline vtkSmartPointer<vtkPoints> GetSamplePointsFor2DCircleDomain(units::quantity<unit::length> pellet_height,
        units::quantity<unit::length> cornea_radius, units::quantity<unit::length> pellet_radius)
{
    vtkSmartPointer<vtkPoints> p_sample_points = vtkSmartPointer<vtkPoints>::New();
    units::quantity<unit::length> reference_length = BaseUnits::Instance()->GetReferenceLengthScale();
    units::quantity<unit::length> cell_spacing(40.0*unit::microns);
    units::quantity<unit::length> sampling_radius = cornea_radius-200e-6*unit::metres;
    unsigned num_cells = (2.0*M_PI*sampling_radius)/cell_spacing +1u;
    double sweep_angle = 2.0*M_PI/num_cells;

    for(unsigned idx=0; idx<num_cells; idx++)
    {
        double this_angle = double(idx)*sweep_angle+M_PI;
        double x_coord = (sampling_radius/reference_length)*std::sin(this_angle);
        double y_coord = (sampling_radius/reference_length)*std::cos(this_angle);
        std::cout << "angle " << this_angle << std::endl;
        std::cout << "x " << x_coord << std::endl;
        std::cout << "y " << y_coord << std::endl;
        p_sample_points->InsertNextPoint(x_coord, y_coord, 0.0);
    }
    return p_sample_points;
}

#endif /*SETUPDOMAINSANDPDES_HPP_*/
