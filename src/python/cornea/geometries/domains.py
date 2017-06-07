"""
A collection of simulation domains
"""

import numpy as np
import microvessel_chaste.geometry
import microvessel_chaste.mesh

def get_2d_planar_domain(parameter_collection, reference_length):  
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    use_pellet = parameter_collection.get_parameter("use pellet").value
    
    domain = microvessel_chaste.geometry.Part2()
    domain_width = 2.0*np.pi*cornea_radius
    domain_height = pellet_height
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint2()
    if not use_pellet:
        domain.AddRectangle(domain_width, 
                            domain_height, 
                            microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0))
        
        domain.AddAttributeToEdgeIfFound(microvessel_chaste.mesh.DimensionalChastePoint2(np.pi*cornea_radius/reference_length,
                pellet_height/reference_length, 0, reference_length), "Pellet Interface", 1.0)
    else:
        points = []
        points.append(microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0, reference_length));
        points.append(microvessel_chaste.mesh.DimensionalChastePoint2(domain_width/reference_length, 0.0, 0.0, reference_length));
        points.append(microvessel_chaste.mesh.DimensionalChastePoint2(domain_width/reference_length, domain_height/reference_length, 0.0,
                reference_length))
        points.append(microvessel_chaste.mesh.DimensionalChastePoint2(domain_width/(2.0*reference_length) + domain_height/(2.0*reference_length), domain_height/reference_length, 0.0,
                reference_length))
        points.append(microvessel_chaste.mesh.DimensionalChastePoint2(domain_width/(2.0*reference_length) - domain_height/(2.0*reference_length), domain_height/reference_length, 0.0,
                reference_length))
        points.append(microvessel_chaste.mesh.DimensionalChastePoint2(0.0, domain_height/reference_length, 0.0, reference_length));

        polygon = microvessel_chaste.geometry.Polygon2(points)
        domain.AddPolygon(polygon)
        domain.AddAttributeToEdgeIfFound(microvessel_chaste.mesh.DimensionalChastePoint2(domain_width/(2.0*reference_length),
                domain_height/reference_length, 0, reference_length), "Pellet Interface", 1.0);
    return domain, holes

def get_3d_planar_domain(parameter_collection, reference_length):  
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    cornea_thickness = parameter_collection.get_parameter("cornea thickness").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    pellet_radius = parameter_collection.get_parameter("pellet radius").value
    pellet_thickness = parameter_collection.get_parameter("pellet thickness").value
    use_pellet = parameter_collection.get_parameter("use pellet").value
    
    domain = microvessel_chaste.geometry.Part3()
    domain_width = 2.0*np.pi*cornea_radius
    domain_height = pellet_height
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint3()
    
    if not use_pellet:
        domain.AddCuboid(domain_width, domain_height, cornea_thickness,
                            microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0))
        for eachFacet in domain.GetFacets():
            distance = eachFacet.GetCentroid().GetDistance(microvessel_chaste.mesh.DimensionalChastePoint3(domain_width/(2.0*reference_length), 
                                                                                       domain_height/reference_length, 
                                                                                       cornea_thickness/(2.0*reference_length),
                                                                                       reference_length))
            if float(distance/reference_length) < 1e-3:
                eachFacet.GetPolygons()[0].AddAttribute("Pellet Interface", 1.0)
    else:
        domain.AddCuboid(domain_width, 
                            domain_height, 
                            cornea_thickness,
                            microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0))
        
        gap = (cornea_thickness - pellet_thickness)/(2.0*reference_length)
        points = []
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width-pellet_radius)/(2.0*reference_length), 
                                                                      domain_height/reference_length, gap, reference_length));
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width+pellet_radius)/(2.0*reference_length), 
                                                                      domain_height/reference_length, gap, reference_length));
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width+pellet_radius)/(2.0*reference_length), 
                                                                      domain_height/reference_length, cornea_thickness/reference_length-gap,reference_length))
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width-pellet_radius)/(2.0*reference_length), 
                                                                      domain_height/reference_length, cornea_thickness/reference_length-gap,reference_length))
        polygon = microvessel_chaste.geometry.Polygon3(points)  
        polygon.AddAttribute("Pellet Interface", 1.0)
        for eachFacet in domain.GetFacets():
            distance = eachFacet.GetCentroid().GetDistance(microvessel_chaste.mesh.DimensionalChastePoint3(domain_width/(2.0*reference_length), 
                                                                                       domain_height/reference_length, 
                                                                                       cornea_thickness/(2.0*reference_length),
                                                                                       reference_length))
            if float(distance/reference_length) < 1e-3:
                domain.AddPolygon(polygon, False, eachFacet)
    return domain, holes

def get_2d_circle_domain(parameter_collection, reference_length): 
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    pellet_radius = parameter_collection.get_parameter("pellet radius").value
    use_pellet = parameter_collection.get_parameter("use pellet").value 
    
    delta = pellet_height-cornea_radius+pellet_radius
    domain = microvessel_chaste.geometry.Part2()
    domain.AddCircle(cornea_radius, 
                       microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0), 24)
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint2()
    if use_pellet:
        polygon = domain.AddCircle(pellet_radius,
                microvessel_chaste.mesh.DimensionalChastePoint2(0.0, -1.0*delta/reference_length, 0.0, reference_length), 24)
        polygon.AddAttributeToAllEdges("Pellet Interface", 1.0)
        domain.AddHoleMarker(microvessel_chaste.mesh.DimensionalChastePoint2(0.0, -1.0*delta/reference_length, 0.0, reference_length))
        holes.append(microvessel_chaste.mesh.DimensionalChastePoint2(0.0, -1.0*delta/reference_length, 0.0, reference_length))
    return domain, holes

def get_3d_circle_domain(parameter_collection, reference_length):  
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    cornea_thickness = parameter_collection.get_parameter("cornea thickness").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    pellet_radius = parameter_collection.get_parameter("pellet radius").value
    pellet_thickness = parameter_collection.get_parameter("pellet thickness").value
    use_pellet = parameter_collection.get_parameter("use pellet").value
    
    delta = pellet_height-cornea_radius+pellet_radius
    domain = microvessel_chaste.geometry.Part3()
    circle = domain.AddCircle(cornea_radius, 
                       microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0), 24)
    domain.Extrude(circle, cornea_thickness)
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint3()
    
    if use_pellet:
        gap = (cornea_thickness - pellet_thickness)/(2.0*reference_length)
        
        pellet = microvessel_chaste.geometry.Part3()
        circle = pellet.AddCircle(pellet_radius,
                microvessel_chaste.mesh.DimensionalChastePoint3(0.0, -1.0*delta/reference_length, 0.0, reference_length), 24)
        pellet.Extrude(circle, pellet_thickness)
        pellet.Translate(microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, gap, reference_length))
        polygons = pellet.GetPolygons() 
        half_height = cornea_thickness/(2.0*reference_length)
        domain.AddHoleMarker(microvessel_chaste.mesh.DimensionalChastePoint3(0.0, -1.0*delta/reference_length, half_height, reference_length))
        holes.append(microvessel_chaste.mesh.DimensionalChastePoint3(0.0, -1.0*delta/reference_length, half_height, reference_length))
        domain.AppendPart(pellet)
        for eachPolygon in polygons:
            eachPolygon.AddAttribute("Pellet Interface", 1.0)
    return domain, holes

def get_3d_hemisphere_domain(parameter_collection, reference_length): 
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    cornea_thickness = parameter_collection.get_parameter("cornea thickness").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    pellet_radius = parameter_collection.get_parameter("pellet radius").value
    pellet_thickness = parameter_collection.get_parameter("pellet thickness").value
    use_pellet = parameter_collection.get_parameter("use pellet").value 
    
    generator = microvessel_chaste.geometry.MappableGridGenerator()
    num_divisions_x = 20
    num_divisions_y = 20
    azimuth_angle = 1.0 * np.pi
    polar_angle = 0.999 * np.pi
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint3()
    domain = generator.GenerateHemisphere(cornea_radius, 
                                                     cornea_thickness,
                                                     num_divisions_x, num_divisions_y, azimuth_angle, polar_angle);

    if use_pellet:
        pellet_domain = microvessel_chaste.geometry.Part3()
    
        gap = (cornea_thickness - pellet_thickness)/(2.0*reference_length)/4.0
        base = cornea_radius/reference_length + gap - cornea_thickness/reference_length
        pellet_domain.AddCylinder(pellet_radius, 
                                  pellet_thickness, 
                                  microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, base, reference_length))
    
        # Rotate the part
        polygons = pellet_domain.GetPolygons()
        rotation_angle = np.arccos(float((pellet_height+pellet_radius)/cornea_radius))
        pellet_centre = microvessel_chaste.mesh.DimensionalChastePoint3(
                0.0, 0.0, base + pellet_thickness/(2.0*reference_length), reference_length)
        pellet_domain.RotateAboutAxis((0, 1, 0), rotation_angle)
        pellet_centre.RotateAboutAxis((0, 1, 0), rotation_angle);
    
        # Add the pellet domain to the cornea
        domain.AppendPart(pellet_domain)
        domain.AddHoleMarker(pellet_centre)
        for eachPolygon in polygons:
            eachPolygon.AddAttribute("Pellet Interface", 1.0)
    return domain, holes

def get_domain(domain_type, parameter_collection, reference_length):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_domain(parameter_collection, reference_length)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_domain(parameter_collection, reference_length)
    
    elif domain_type == "Circle 2D":
        return get_2d_circle_domain(parameter_collection, reference_length)
    
    elif domain_type == "Circle 3D":
        return get_3d_circle_domain(parameter_collection, reference_length)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_domain(parameter_collection, reference_length)
    
if __name__ == '__main__':
    
    import chaste.core
    import cornea.parameters.default_parameters
    
    work_dir = "Python/Cornea/TestDomains/"
    reference_length = 1.e-6 * metre()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    
    # Write the domain
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), True)
        print file_handler.GetOutputDirectoryFullPath()
        domain, holes = get_domain(eachDomainType, parameter_collection, reference_length)
        domain.Write(file_handler.GetOutputDirectoryFullPath() + "part.vtp",
                     microvessel_chaste.geometry.GeometryFormat.VTP, True)
        if "2" in eachDomainType:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter2()
            generator.SetMaxElementArea(parameter_collection.get_parameter("element area 2d").value)
        else:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter3()
            generator.SetMaxElementArea(parameter_collection.get_parameter("element area 3d").value)
        generator.SetDomain(domain)
        generator.Update()
        mesh_writer.SetMesh(generator.GetMesh(), True)
        mesh_writer.SetFileName(file_handler.GetOutputDirectoryFullPath() + "mesh")
        mesh_writer.Write()
        
    parameter_collection.get_parameter("use pellet").value = True
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), False)
        print file_handler.GetOutputDirectoryFullPath()
        domain, holes = get_domain(eachDomainType, parameter_collection, reference_length)
        domain.Write(file_handler.GetOutputDirectoryFullPath() + "part_with_pellet.vtp",
                     microvessel_chaste.geometry.GeometryFormat.VTP, True)
        
        if "2" in eachDomainType:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter2()
            generator.SetMaxElementArea(parameter_collection.get_parameter("element area 2d").value)
        else:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter3()
            generator.SetMaxElementArea(parameter_collection.get_parameter("element area 3d").value) 
        generator.SetDomain(domain)
        generator.SetHoles(holes)
        generator.Update()
        mesh = generator.GetMesh()
        mesh_writer.SetMesh(mesh, True)
        mesh_writer.SetFileName(file_handler.GetOutputDirectoryFullPath() + "mesh_with_pellet")
        mesh_writer.Write()
        