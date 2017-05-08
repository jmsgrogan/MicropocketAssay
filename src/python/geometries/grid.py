import numpy as np
import microvessel_chaste.mesh
from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience

def get_2d_planar_grid(domain, domain_dimensions):  
    
    if not domain_dimensions["use finite pellet width"]:
        grid = microvessel_chaste.mesh.RegularGrid2()
        grid.GenerateFromPart(domain, domain_dimensions["grid spacing"])
    else:
        generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
        generator.SetDomain(domain)
        generator.SetMaxElementArea(5e4*(1.e-18*metre_cubed()))
        generator.Update()  
        grid = generator.GetMesh()     
    return grid

def get_3d_planar_grid(domain, domain_dimensions):  

    if not domain_dimensions["use finite pellet width"]:    
        grid = microvessel_chaste.mesh.RegularGrid3()
        grid.GenerateFromPart(domain, domain_dimensions["grid spacing"])
    else:
        generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
        generator.SetDomain(domain)
        generator.SetMaxElementArea(5e4*(1.e-18*metre_cubed()))
        generator.Update()  
        grid = generator.GetMesh()   
    return grid

def get_2d_circle_grid(domain, domain_dimensions, holes = None):  
    
    generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
    generator.SetDomain(domain)
    generator.SetMaxElementArea(5e4*(1.e-18*metre_cubed()))
    
    if holes is not None and len(holes) is not 0:
        generator.SetHoles(holes)
    generator.Update()
    return generator.GetMesh()

def get_3d_circle_grid(domain, domain_dimensions, holes = None):  
    
    generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
    generator.SetDomain(domain)
    generator.SetMaxElementArea(5e4*(1.e-18*metre_cubed()))
    if holes is not None and len(holes) is not 0:
        generator.SetHoles(holes)
    generator.Update()
    return generator.GetMesh()

def get_3d_hemisphere_grid(domain, domain_dimensions, holes = None):  
    
    generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
    generator.SetDomain(domain)
    generator.SetMaxElementArea(5e4*(1.e-18*metre_cubed()))
    if holes is not None and len(holes) is not 0:
        generator.SetHoles(holes)
    generator.Update()
    mesh = generator.GetMesh()

    if holes is not None and len(holes) is not 0:
        rotation_angle = np.pi/25.0
        domain.RotateAboutAxis((0, 1, 0), -rotation_angle)
        mesh.Rotate((0, 1, 0), -4.0*rotation_angle)  
    return mesh

def get_grid(domain_type, domain, domain_dimensions, holes = None):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_grid(domain, domain_dimensions)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_grid(domain, domain_dimensions)
    
    elif domain_type == "Circle 2D":
        return get_2d_circle_grid(domain, domain_dimensions, holes)
    
    elif domain_type == "Circle 3D":
        return get_3d_circle_grid(domain, domain_dimensions, holes)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_grid(domain, domain_dimensions, holes)