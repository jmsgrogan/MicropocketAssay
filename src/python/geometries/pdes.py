import numpy as np
import microvessel_chaste.pde

def get_transient_pde(domain_type, domain_dimensions, pde_parameters, reference_length, reference_concentration):
    
    if "2" in domain_type:
        pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde2()
    else:
        pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde3()
        
    pde.SetIsotropicDiffusionConstant(pde_parameters["vegf diffusivity"])
    pde.SetContinuumLinearInUTerm(pde_parameters["vegf decay rate"])
    pde.SetCurrentVegfInPellet(pde_parameters["pellet concentration"])
    pde.SetPelletBindingConstant(pde_parameters["pellet binding constant"])
    pde.SetPelletDepth(domain_dimensions["pellet thickness"])
    
    surface_area = 2.0*np.pi*domain_dimensions["pellet radius"]*domain_dimensions["pellet thickness"]
    surface_area += 2.0*np.pi*domain_dimensions["pellet radius"]*domain_dimensions["pellet radius"]
    
    pde.SetPelletSurfaceArea(surface_area)
    pde.SetCorneaPelletPermeability(0.002*pde.GetCorneaPelletPermeability())
    return pde

def get_2d_planar_fixed_gradient(grid, domain_dimensions, pde_parameters, reference_length, reference_concentration):
    
    funciton_map = microvessel_chaste.pde.FunctionMap2()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        dimless_pellet_height = (domain_dimensions["pellet height"]/reference_length)
        normalized_distance = y_loc/dimless_pellet_height
        vegf_field.append(normalized_distance*(pde_parameters["pellet concentration"]/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map
    
def get_3d_planar_fixed_gradient(grid, domain_dimensions, pde_parameters, reference_length, reference_concentration):
    
    funciton_map = microvessel_chaste.pde.FunctionMap3()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        dimless_pellet_height = (domain_dimensions["pellet height"]/reference_length)
        normalized_distance = y_loc/dimless_pellet_height
        vegf_field.append(normalized_distance*(pde_parameters["pellet concentration"]/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map

def get_2d_circle_fixed_gradient(grid, domain_dimensions, pde_parameters, reference_length, reference_concentration):
    
    funciton_map = microvessel_chaste.pde.FunctionMap2()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        x_loc = grid.GetPoint(idx).GetLocation(reference_length)[0]
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        radius = np.sqrt(x_loc*x_loc + y_loc*y_loc)
        dimless_pellet_height = (domain_dimensions["cornea radius"]/reference_length)
        normalized_distance = float(radius/dimless_pellet_height)
        vegf_field.append((1.0-normalized_distance)*float(pde_parameters["pellet concentration"]/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map

def get_3d_circle_fixed_gradient(grid, domain_dimensions, pde_parameters, reference_length, reference_concentration):
    
    funciton_map = microvessel_chaste.pde.FunctionMap3()
    funciton_map.SetGrid(grid)
    
    vegf_field = []
    for idx in range(grid.GetNumberOfPoints()):
        x_loc = grid.GetPoint(idx).GetLocation(reference_length)[0]
        y_loc = grid.GetPoint(idx).GetLocation(reference_length)[1]
        radius = np.sqrt(x_loc*x_loc + y_loc*y_loc)
        dimless_pellet_height = (domain_dimensions["cornea radius"]/reference_length)
        normalized_distance = float(radius/dimless_pellet_height)
        vegf_field.append((1.0-normalized_distance)*float(pde_parameters["pellet concentration"]/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map

def get_3d_hemisphere_fixed_gradient(grid, domain_dimensions, pde_parameters, reference_length, reference_concentration):
    
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
        vegf_field.append(frac*float(pde_parameters["pellet concentration"]/reference_concentration))
    funciton_map.SetLabel("vegf")
    funciton_map.UpdateSolution(vegf_field)
    return funciton_map
    
def get_pde_fixed_gradient(domain_type, grid, domain_dimensions, pde_parameters, 
                       reference_length, reference_concentration):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_fixed_gradient(grid, domain_dimensions, pde_parameters, 
                       reference_length, reference_concentration)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_fixed_gradient(grid, domain_dimensions, pde_parameters, 
                       reference_length, reference_concentration)
      
    elif domain_type == "Circle 2D":
        return get_2d_circle_fixed_gradient(grid, domain_dimensions, pde_parameters, 
                       reference_length, reference_concentration)
      
    elif domain_type == "Circle 3D":
        return get_3d_circle_fixed_gradient(grid, domain_dimensions, pde_parameters, 
                       reference_length, reference_concentration)
        
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_fixed_gradient(grid, domain_dimensions, pde_parameters, 
                       reference_length, reference_concentration)