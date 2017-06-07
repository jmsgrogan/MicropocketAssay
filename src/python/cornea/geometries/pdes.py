import numpy as np
import microvessel_chaste.pde
from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience

def get_transient_pde(domain_type, parameter_collection, reference_length, reference_concentration):
    
    if "2" in domain_type:
        pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde2_2()
    else:
        pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde3_3()
        
    vegf_diffusivity = parameter_collection.get_parameter("vegf diffusivity").value
    vegf_decay_rate = parameter_collection.get_parameter("vegf decay rate").value
    pellet_concentration = parameter_collection.get_parameter("pellet concentration").value
    pellet_binding_constant = parameter_collection.get_parameter("pellet binding constant").value
    pellet_thickness = parameter_collection.get_parameter("pellet thickness").value 
    pellet_radius = parameter_collection.get_parameter("pellet radius").value 
    vegf_blood_concentration = parameter_collection.get_parameter("vegf blood concentration").value 
    vegf_permeability = parameter_collection.get_parameter("vessel permeability").value 
    uptake_rate_per_cell = parameter_collection.get_parameter("uptake rate per cell").value
    include_sink = parameter_collection.get_parameter("include vessel sink").value
    
    pde.SetIsotropicDiffusionConstant(vegf_diffusivity)
    pde.SetContinuumLinearInUTerm(vegf_decay_rate)
    pde.SetCurrentVegfInPellet(pellet_concentration)
    pde.SetPelletBindingConstant(pellet_binding_constant*dimensionless())
    pde.SetPelletDepth(pellet_thickness)
    
    surface_area = 2.0*np.pi*pellet_radius*pellet_thickness
    surface_area += 2.0*np.pi*pellet_radius*pellet_radius
    pde.SetPelletSurfaceArea(surface_area)
    pde.SetCorneaPelletPermeability(0.002*pde.GetCorneaPelletPermeability())
    
    if include_sink:
        if "2" in domain_type:
            vessel_vegf_sink = microvessel_chaste.pde.VesselBasedDiscreteSource2()
        else:
            vessel_vegf_sink = microvessel_chaste.pde.VesselBasedDiscreteSource3()
        vessel_vegf_sink.SetReferenceConcentration(vegf_blood_concentration)
        vessel_vegf_sink.SetVesselPermeability(vegf_permeability)
        vessel_vegf_sink.SetReferenceHaematocrit(0.45*dimensionless());
        vessel_vegf_sink.SetUptakeRatePerCell(-1.0*uptake_rate_per_cell);
        pde.AddDiscreteSource(vessel_vegf_sink)  
    return pde

def get_2d_planar_solver(pde, grid, domain, parameter_collection, reference_length, reference_concentration):
    
    finite_width = parameter_collection.get_parameter("use finite pellet width").value
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    if not finite_width:
        solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteDifferenceSolver2()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.SetLabel("vegf")
    else:
        boundary_condition = microvessel_chaste.pde.DiscreteContinuumBoundaryCondition2()
        boundary_condition.SetValue(pellet_conc)
        boundary_condition.SetType(microvessel_chaste.pde.BoundaryConditionType.EDGE)
        boundary_condition.SetIsRobin(True)
        boundary_condition.SetLabel("Pellet Interface")
        boundary_condition.SetDomain(domain)
        
        solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteElementSolver2()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.SetLabel("vegf")
        solver.AddBoundaryCondition(boundary_condition) 
    return solver

def get_2d_circle_solver(pde, grid, domain, parameter_collection, reference_length, reference_concentration):
    
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    boundary_condition = microvessel_chaste.pde.DiscreteContinuumBoundaryCondition2()
    boundary_condition.SetValue(pellet_conc)
    boundary_condition.SetType(microvessel_chaste.pde.BoundaryConditionType.EDGE)
    boundary_condition.SetIsRobin(True)
    boundary_condition.SetLabel("Pellet Interface")
    boundary_condition.SetDomain(domain)
    
    solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteElementSolver2()
    solver.SetGrid(grid)
    solver.SetPde(pde)
    solver.SetLabel("vegf")
    solver.AddBoundaryCondition(boundary_condition)
    return solver
    
def get_3d_circle_solver(pde, grid, domain, parameter_collection, reference_length, reference_concentration):
    
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
        
    boundary_condition = microvessel_chaste.pde.DiscreteContinuumBoundaryCondition3()
    boundary_condition.SetValue(pellet_conc)
    boundary_condition.SetType(microvessel_chaste.pde.BoundaryConditionType.POLYGON)
    boundary_condition.SetIsRobin(True)
    boundary_condition.SetLabel("Pellet Interface")
    boundary_condition.SetDomain(domain)
    
    solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteElementSolver3()
    solver.SetGrid(grid)
    solver.SetPde(pde)
    solver.SetLabel("vegf")
    solver.AddBoundaryCondition(boundary_condition) 
    return solver
    
def get_3d_hemisphere_solver(pde, grid, domain, parameter_collection, reference_length, reference_concentration):
    
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    boundary_condition = microvessel_chaste.pde.DiscreteContinuumBoundaryCondition3()
    boundary_condition.SetValue(pellet_conc)
    boundary_condition.SetType(microvessel_chaste.pde.BoundaryConditionType.POLYGON)
    boundary_condition.SetIsRobin(True)
    boundary_condition.SetLabel("Pellet Interface")
    boundary_condition.SetDomain(domain)
    
    solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteElementSolver3()
    solver.SetGrid(grid)
    solver.SetPde(pde)
    solver.SetLabel("vegf")
    solver.AddBoundaryCondition(boundary_condition)   
    return solver     

def get_3d_planar_solver(pde, grid, domain, parameter_collection, reference_length, reference_concentration):

    finite_width = parameter_collection.get_parameter("use finite pellet width").value
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    if not finite_width:    
        solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteDifferenceSolver3()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.SetLabel("vegf")
    else:
        boundary_condition = microvessel_chaste.pde.DiscreteContinuumBoundaryCondition3()
        boundary_condition.SetValue(pellet_conc)
        boundary_condition.SetType(microvessel_chaste.pde.BoundaryConditionType.POLYGON)
        boundary_condition.SetIsRobin(True)
        boundary_condition.SetLabel("Pellet Interface")
        boundary_condition.SetDomain(domain)
        
        solver = microvessel_chaste.pde.CoupledLumpedSystemFiniteElementSolver3()
        solver.SetGrid(grid)
        solver.SetPde(pde)
        solver.SetLabel("vegf")
        solver.AddBoundaryCondition(boundary_condition)         
    return solver

def get_2d_planar_fixed_gradient(grid, parameter_collection, reference_length, reference_concentration):
    
    pellet_height = parameter_collection.get_parameter("pellet height").value
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    funciton_map = microvessel_chaste.pde.FunctionMap2()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        dimless_pellet_height = (pellet_height/reference_length)
        normalized_distance = y_loc/dimless_pellet_height
        vegf_field.append(normalized_distance*(pellet_conc/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map
    
def get_3d_planar_fixed_gradient(grid, parameter_collection, reference_length, reference_concentration):
    
    pellet_height = parameter_collection.get_parameter("pellet height").value
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    funciton_map = microvessel_chaste.pde.FunctionMap3()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        dimless_pellet_height = (pellet_height/reference_length)
        normalized_distance = y_loc/dimless_pellet_height
        vegf_field.append(normalized_distance*(pellet_conc/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map

def get_2d_circle_fixed_gradient(grid, parameter_collection, reference_length, reference_concentration):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    funciton_map = microvessel_chaste.pde.FunctionMap2()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        x_loc = grid.GetPoint(idx).GetLocation(reference_length)[0]
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        radius = np.sqrt(x_loc*x_loc + y_loc*y_loc)
        dimless_pellet_height = (cornea_radius/reference_length)
        normalized_distance = float(radius/dimless_pellet_height)
        vegf_field.append((1.0-normalized_distance)*float(pellet_conc/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map

def get_3d_circle_fixed_gradient(grid, parameter_collection, reference_length, reference_concentration):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    funciton_map = microvessel_chaste.pde.FunctionMap3()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        x_loc = grid.GetPoint(idx).GetLocation(reference_length)[0]
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        radius = np.sqrt(x_loc*x_loc + y_loc*y_loc)
        dimless_pellet_height = (cornea_radius/reference_length)
        normalized_distance = float(radius/dimless_pellet_height)
        vegf_field.append((1.0-normalized_distance)*float(pellet_conc/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map

def get_3d_hemisphere_fixed_gradient(grid, parameter_collection, reference_length, reference_concentration):
    
    pellet_conc = parameter_collection.get_parameter("pellet concentration").value
    
    funciton_map = microvessel_chaste.pde.FunctionMap3()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        x_loc = grid.GetPoint(idx).GetLocation(reference_length)[0]
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        z_loc = grid.GetPoint(idx).GetLocation(reference_length)[2]
        radius = np.sqrt(x_loc*x_loc + y_loc*y_loc)
        angle = np.arctan(z_loc/radius)
        frac = angle/(np.pi/2.0)
        vegf_field.append(frac*float(pellet_conc/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map
    
def get_pde_fixed_gradient(domain_type, grid, parameter_collection, 
                       reference_length, reference_concentration):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_fixed_gradient(grid, parameter_collection, 
                       reference_length, reference_concentration)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_fixed_gradient(grid, parameter_collection, 
                       reference_length, reference_concentration)
      
    elif domain_type == "Circle 2D":
        return get_2d_circle_fixed_gradient(grid, parameter_collection, 
                       reference_length, reference_concentration)
      
    elif domain_type == "Circle 3D":
        return get_3d_circle_fixed_gradient(grid, parameter_collection, 
                       reference_length, reference_concentration)
        
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_fixed_gradient(grid, parameter_collection, 
                       reference_length, reference_concentration)
        
def get_pde_solver(domain_type, pde, grid, domain, parameter_collection, 
                       reference_length, reference_concentration):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_solver(pde, grid, domain, parameter_collection, 
                       reference_length, reference_concentration)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_solver(pde, grid, domain, parameter_collection, 
                       reference_length, reference_concentration)
      
    elif domain_type == "Circle 2D":
        return get_2d_circle_solver(pde, grid, domain, parameter_collection, 
                       reference_length, reference_concentration)
      
    elif domain_type == "Circle 3D":
        return get_3d_circle_solver(pde, grid, domain, parameter_collection, 
                       reference_length, reference_concentration)
        
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_solver(pde, grid, domain, parameter_collection, 
                       reference_length, reference_concentration)