import microvessel_chaste.mesh

def get_2d_planar_grid(domain, parameter_collection):  
    
    use_finite_pellet = parameter_collection.get_parameter("use finite pellet width").value
    grid_spacing = parameter_collection.get_parameter("grid spacing").value
    mesh_area = parameter_collection.get_parameter("element area 2d").value
    
    if not use_finite_pellet:
        grid = microvessel_chaste.mesh.RegularGrid2()
        grid.GenerateFromPart(domain, grid_spacing)
    else:
        generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
        generator.SetDomain(domain)
        generator.SetMaxElementArea(mesh_area)
        generator.Update()  
        grid = generator.GetMesh()     
    return grid

def get_3d_planar_grid(domain, parameter_collection):  
    
    use_finite_pellet = parameter_collection.get_parameter("use finite pellet width").value
    grid_spacing = parameter_collection.get_parameter("grid spacing").value
    mesh_area = parameter_collection.get_parameter("element area 2d").value

    if not use_finite_pellet:    
        grid = microvessel_chaste.mesh.RegularGrid3()
        grid.GenerateFromPart(domain, grid_spacing)
    else:
        generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
        generator.SetDomain(domain)
        generator.SetMaxElementArea(mesh_area)
        generator.Update()  
        grid = generator.GetMesh()   
    return grid

def get_2d_circle_grid(domain, parameter_collection, holes = None):  
    
    mesh_area = parameter_collection.get_parameter("element area 2d").value
    
    generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
    generator.SetDomain(domain)
    generator.SetMaxElementArea(mesh_area)
    
    if holes is not None and len(holes) is not 0:
        generator.SetHoles(holes)
    generator.Update()
    return generator.GetMesh()

def get_3d_circle_grid(domain, parameter_collection, holes = None): 
    
    mesh_area = parameter_collection.get_parameter("element area 3d").value 
    
    generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
    generator.SetDomain(domain)
    generator.SetMaxElementArea(mesh_area)
    if holes is not None and len(holes) is not 0:
        generator.SetHoles(holes)
    generator.Update()
    return generator.GetMesh()

def get_3d_hemisphere_grid(domain, parameter_collection, holes = None): 
    
    mesh_area = parameter_collection.get_parameter("element area 3d").value  
    
    generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
    generator.SetDomain(domain)
    generator.SetMaxElementArea(mesh_area)
    if holes is not None and len(holes) is not 0:
        generator.SetHoles(holes)
    generator.Update()
    mesh = generator.GetMesh()

    return mesh

def get_grid(domain_type, domain, parameter_collection, holes = None):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_grid(domain, parameter_collection)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_grid(domain, parameter_collection)
    
    elif domain_type == "Circle 2D":
        return get_2d_circle_grid(domain, parameter_collection, holes)
    
    elif domain_type == "Circle 3D":
        return get_3d_circle_grid(domain, parameter_collection, holes)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_grid(domain, parameter_collection, holes)