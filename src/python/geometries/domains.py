"""
A collection of simulation domains
"""

import numpy as np
import microvessel_chaste.geometry
import microvessel_chaste.mesh

def get_2d_planar_domain(domain_dimensions, reference_length):  
    
    domain = microvessel_chaste.geometry.Part2()
    domain_width = 2.0*np.pi*domain_dimensions["cornea radius"]
    domain_height = domain_dimensions["pellet height"]
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint2()
    
    if not domain_dimensions["use pellet"]:
        domain.AddRectangle(domain_width, 
                            domain_height, 
                            microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0))
        
        domain.AddAttributeToEdgeIfFound(microvessel_chaste.mesh.DimensionalChastePoint2(np.pi*domain_dimensions["cornea radius"]/reference_length,
                domain_dimensions["pellet height"]/reference_length, 0, reference_length), "Pellet Interface", 1.0)
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

def get_3d_planar_domain(domain_dimensions, reference_length):  
    
    domain = microvessel_chaste.geometry.Part3()
    domain_width = 2.0*np.pi*domain_dimensions["cornea radius"]
    domain_height = domain_dimensions["pellet height"]
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint3()
    
    if not domain_dimensions["use pellet"]:
        domain.AddCuboid(domain_width, domain_height, domain_dimensions["cornea thickness"],
                            microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0))
        for eachFacet in domain.GetFacets():
            distance = eachFacet.GetCentroid().GetDistance(microvessel_chaste.mesh.DimensionalChastePoint3(domain_width/(2.0*reference_length), 
                                                                                       domain_height/reference_length, 
                                                                                       domain_dimensions["cornea thickness"]/(2.0*reference_length),
                                                                                       reference_length))
            if float(distance/reference_length) < 1e-3:
                eachFacet.GetPolygons()[0].AddAttribute("Pellet Interface", 1.0)
    else:
        domain.AddCuboid(domain_width, 
                            domain_height, 
                            domain_dimensions["cornea thickness"],
                            microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0))
        
        gap = (domain_dimensions["cornea thickness"] - domain_dimensions["pellet thickness"])/(2.0*reference_length)
        
        points = []
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width-domain_dimensions["pellet radius"])/(2.0*reference_length), 
                                                                      domain_height/reference_length, gap, reference_length));
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width+domain_dimensions["pellet radius"])/(2.0*reference_length), 
                                                                      domain_height/reference_length, gap, reference_length));
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width+domain_dimensions["pellet radius"])/(2.0*reference_length), 
                                                                      domain_height/reference_length, domain_dimensions["cornea thickness"]/reference_length-gap,reference_length))
        points.append(microvessel_chaste.mesh.DimensionalChastePoint3((domain_width-domain_dimensions["pellet radius"])/(2.0*reference_length), 
                                                                      domain_height/reference_length, domain_dimensions["cornea thickness"]/reference_length-gap,reference_length))
        polygon = microvessel_chaste.geometry.Polygon3(points)  
        polygon.AddAttribute("Pellet Interface", 1.0)
        for eachFacet in domain.GetFacets():
            distance = eachFacet.GetCentroid().GetDistance(microvessel_chaste.mesh.DimensionalChastePoint3(domain_width/(2.0*reference_length), 
                                                                                       domain_height/reference_length, 
                                                                                       domain_dimensions["cornea thickness"]/(2.0*reference_length),
                                                                                       reference_length))
            if float(distance/reference_length) < 1e-3:
                domain.AddPolygon(polygon, False, eachFacet)
    return domain, holes

def get_2d_circle_domain(domain_dimensions, reference_length):  
    
    delta = domain_dimensions["pellet height"]-domain_dimensions["cornea radius"]+domain_dimensions["pellet radius"]
    domain = microvessel_chaste.geometry.Part2()
    domain.AddCircle(domain_dimensions["cornea radius"], 
                       microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 0.0, 0.0), 24)
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint2()
    if domain_dimensions["use pellet"]:
        polygon = domain.AddCircle(domain_dimensions["pellet radius"],
                microvessel_chaste.mesh.DimensionalChastePoint2(0.0, -1.0*delta/reference_length, 0.0, reference_length), 24)
        polygon.AddAttributeToAllEdges("Pellet Interface", 1.0)
        domain.AddHoleMarker(microvessel_chaste.mesh.DimensionalChastePoint2(0.0, -1.0*delta/reference_length, 0.0, reference_length))
        holes.append(microvessel_chaste.mesh.DimensionalChastePoint2(0.0, -1.0*delta/reference_length, 0.0, reference_length))
    return domain, holes

def get_3d_circle_domain(domain_dimensions, reference_length):  
    
    delta = domain_dimensions["pellet height"]-domain_dimensions["cornea radius"]+domain_dimensions["pellet radius"]
    domain = microvessel_chaste.geometry.Part3()
    circle = domain.AddCircle(domain_dimensions["cornea radius"], 
                       microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, 0.0), 24)
    domain.Extrude(circle, domain_dimensions["cornea thickness"])
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint3()
    
    if domain_dimensions["use pellet"]:
        gap = (domain_dimensions["cornea thickness"] - domain_dimensions["pellet thickness"])/(2.0*reference_length)
        
        pellet = microvessel_chaste.geometry.Part3()
        circle = pellet.AddCircle(domain_dimensions["pellet radius"],
                microvessel_chaste.mesh.DimensionalChastePoint3(0.0, -1.0*delta/reference_length, 0.0, reference_length), 24)
        pellet.Extrude(circle, domain_dimensions["pellet thickness"])
        pellet.Translate(microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, gap, reference_length))
        polygons = pellet.GetPolygons() 
        half_height = domain_dimensions["cornea thickness"]/(2.0*reference_length)
        domain.AddHoleMarker(microvessel_chaste.mesh.DimensionalChastePoint3(0.0, -1.0*delta/reference_length, half_height, reference_length))
        holes.append(microvessel_chaste.mesh.DimensionalChastePoint3(0.0, -1.0*delta/reference_length, half_height, reference_length))
        domain.AppendPart(pellet)
        for eachPolygon in polygons:
            eachPolygon.AddAttribute("Pellet Interface", 1.0)
    
    return domain, holes

def get_3d_hemisphere_domain(domain_dimensions, reference_length):  
    
    generator = microvessel_chaste.geometry.MappableGridGenerator()
    num_divisions_x = 20
    num_divisions_y = 20
    azimuth_angle = 1.0 * np.pi
    polar_angle = 0.999 * np.pi
    holes = microvessel_chaste.mesh.VectorDimensionalChastePoint3()
    domain = generator.GenerateHemisphere(domain_dimensions["cornea radius"], 
                                                     domain_dimensions["cornea thickness"],
                                                     num_divisions_x, num_divisions_y, azimuth_angle, polar_angle);

    if domain_dimensions["use pellet"]:
        pellet_domain = microvessel_chaste.geometry.Part3()
    
        gap = (domain_dimensions["cornea thickness"] - domain_dimensions["pellet thickness"])/(2.0*reference_length)/4.0
        base = domain_dimensions["cornea radius"]/reference_length + gap - domain_dimensions["cornea thickness"]/reference_length
        pellet_domain.AddCylinder(domain_dimensions["pellet radius"], 
                                  domain_dimensions["pellet thickness"], 
                                  microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 0.0, base, reference_length))
    
        # Rotate the part
        polygons = pellet_domain.GetPolygons()
        rotation_angle = np.pi/25.0
        domain.RotateAboutAxis((0, 1, 0), rotation_angle)
        #p_vegf_domain->RotateAboutAxis(rotation_axis, rotation_angle);
        #pellet_centre.RotateAboutAxis(rotation_axis, rotation_angle);
    
        # Add the pellet domain to the cornea
        domain.AppendPart(pellet_domain)
        domain.AddHoleMarker(
            microvessel_chaste.mesh.DimensionalChastePoint3(
                0.0, 0.0, base + domain_dimensions["pellet thickness"]/(2.0*reference_length), reference_length))
        holes.append(microvessel_chaste.mesh.DimensionalChastePoint3(
                0.0, 0.0, base + domain_dimensions["pellet thickness"]/(2.0*reference_length), reference_length))
        
        for eachPolygon in polygons:
            eachPolygon.AddAttribute("Pellet Interface", 1.0)
    return domain, holes

def get_domain(domain_type, domain_dimensions, reference_length):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_domain(domain_dimensions, reference_length)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_domain(domain_dimensions, reference_length)
    
    elif domain_type == "Circle 2D":
        return get_2d_circle_domain(domain_dimensions, reference_length)
    
    elif domain_type == "Circle 3D":
        return get_3d_circle_domain(domain_dimensions, reference_length)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_domain(domain_dimensions, reference_length)
    
if __name__ == '__main__':
    
    import chaste.core
    from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience
    
    work_dir = "Python/Cornea/TestDomains/"
    reference_length = 1.e-6 * metre()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    domain_dimensions = {"pellet height" : 1.0e-3*metre(),
                     "cornea radius" : 1.3e-3*metre(),
                     "cornea thickness" : 100.0e-6*metre(),
                     "pellet thickness" : 40.0e-6*metre(),
                     "grid spacing" : 40.0e-6*metre(),
                     "node spacing" : 40.0e-6*metre(),
                     "limbal offset" : 200.0e-6*metre(),
                     "density grid spacing" : 40.0e-6*metre(),
                     "sample spacing x" : 20.0e-6*metre(),
                     "sample spacing y" : 20.0e-6*metre(),
                     "use pellet" : False,
                     "pellet radius" : 300.0e-6*metre(),
                     }
    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    
    # Write the domain
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), True)
        domain, holes = get_domain(eachDomainType, domain_dimensions, reference_length)
        domain.Write(file_handler.GetOutputDirectoryFullPath() + "part.vtp",
                     microvessel_chaste.geometry.GeometryFormat.VTP, True)
         
        if "2" in eachDomainType:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter2()
        else:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter3()
        generator.SetDomain(domain)
        generator.SetMaxElementArea(1e4*(1.e-18*metre_cubed()))
        generator.Update()
        mesh_writer.SetMesh(generator.GetMesh())
        mesh_writer.SetFileName(file_handler.GetOutputDirectoryFullPath() + "mesh")
        mesh_writer.Write()
        
    domain_dimensions["use pellet"] = True
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), False)
        domain, holes = get_domain(eachDomainType, domain_dimensions, reference_length)
        domain.Write(file_handler.GetOutputDirectoryFullPath() + "part_with_pellet.vtp",
                     microvessel_chaste.geometry.GeometryFormat.VTP, True)
        
        if "2" in eachDomainType:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator2_2()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter2()
        else:
            generator = microvessel_chaste.mesh.DiscreteContinuumMeshGenerator3_3()
            mesh_writer = microvessel_chaste.mesh.MultiFormatMeshWriter3()
        generator.SetMaxElementArea(1e4*(1.e-18*metre_cubed()))   

        generator.SetDomain(domain)
        generator.SetHoles(holes)
        generator.Update()
        mesh = generator.GetMesh()
        if "Hemisphere 3D" in eachDomainType:
            rotation_angle = np.pi/25.0
            domain.RotateAboutAxis((0, 1, 0), -rotation_angle)
            mesh.Rotate((0, 1, 0), -4.0*rotation_angle)  
        mesh_writer.SetMesh(mesh)
        mesh_writer.SetFileName(file_handler.GetOutputDirectoryFullPath() + "mesh_with_pellet")
        mesh_writer.Write()
        