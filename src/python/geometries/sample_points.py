import numpy as np
import vtk

def get_2d_planar_sample_points(domain_dimensions, reference_length):

    domain_width = 2.0*np.pi*domain_dimensions["cornea radius"]
    num_sample_points_x = int(float(domain_width/domain_dimensions["sample spacing x"])) + 1
    num_sample_points_y = int(float((domain_dimensions["pellet height"]-domain_dimensions["limbal offset"])/domain_dimensions["sample spacing y"])) + 1

    sample_lines = []
    for idx in range(num_sample_points_y):
        sample_points = vtk.vtkPoints()
        dimless_sample_spacing_x = domain_dimensions["sample spacing x"]/reference_length
        dimless_sample_spacing_y = domain_dimensions["sample spacing y"]/reference_length
        dimless_limbal_offset = domain_dimensions["limbal offset"]/reference_length
        for jdx in range(num_sample_points_x):
            sample_points.InsertNextPoint(float(jdx)*dimless_sample_spacing_x,
                                          float(idx)*dimless_sample_spacing_y + dimless_limbal_offset,
                                          0.0)
        sample_lines.append(sample_points)
    return sample_lines

def get_3d_planar_sample_points(domain_dimensions, reference_length):

    domain_width = 2.0*np.pi*domain_dimensions["cornea radius"]
    num_sample_points_x = int(float(domain_width/domain_dimensions["sample spacing x"])) + 1
    num_sample_points_y = int(float((domain_dimensions["pellet height"]-domain_dimensions["limbal offset"])/domain_dimensions["sample spacing y"])) + 1
    num_sample_points_z = int(float(domain_dimensions["cornea thickness"]/domain_dimensions["sample spacing z"])) + 1
    
    sample_lines = []
    for kdx in range(num_sample_points_z):
        for idx in range(num_sample_points_y):
            sample_points = vtk.vtkPoints()
            dimless_sample_spacing_x = domain_dimensions["sample spacing x"]/reference_length
            dimless_sample_spacing_y = domain_dimensions["sample spacing y"]/reference_length
            dimless_sample_spacing_z = domain_dimensions["sample spacing z"]/reference_length
            dimless_limbal_offset = domain_dimensions["limbal offset"]/reference_length
            for jdx in range(num_sample_points_x):
                sample_points.InsertNextPoint(float(jdx)*dimless_sample_spacing_x,
                                              float(idx)*dimless_sample_spacing_y + dimless_limbal_offset,
                                              float(kdx)*dimless_sample_spacing_z)
            sample_lines.append(sample_points)
    return sample_lines

def get_2d_circle_sample_points(domain_dimensions, reference_length):

    domain_width = 2.0*np.pi*domain_dimensions["cornea radius"]
    num_sample_points_x = int(float(domain_width/domain_dimensions["sample spacing x"])) + 1
    num_sample_points_y = int(float(domain_dimensions["pellet height"]/domain_dimensions["sample spacing y"])) + 1

    sample_lines = []

    for idx in range(num_sample_points_y):
        sampling_radius = domain_dimensions["cornea radius"]-domain_dimensions["limbal offset"]-float(idx)*domain_dimensions["sample spacing y"]
        num_nodes = int(float((2.0*np.pi*sampling_radius)/domain_dimensions["node spacing"])) + 1
        sweep_angle = 2.0*np.pi/num_nodes
        sample_points = vtk.vtkPoints()
        for jdx in range(num_sample_points_x):
            this_angle = float(jdx)*sweep_angle+np.pi
            x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
            y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
            sample_points.InsertNextPoint(x_coord,y_coord, 0.0)
        sample_lines.append(sample_points)
    return sample_lines

def get_3d_circle_sample_points(domain_dimensions, reference_length):

    domain_width = 2.0*np.pi*domain_dimensions["cornea radius"]
    num_sample_points_x = int(float(domain_width/domain_dimensions["sample spacing x"])) + 1
    num_sample_points_y = int(float(domain_dimensions["pellet height"]/domain_dimensions["sample spacing y"])) + 1
    num_sample_points_z = int(float(domain_dimensions["cornea thickness"]/domain_dimensions["sample spacing z"])) + 1
    sample_lines = []

    for kdx in range(num_sample_points_z):
        for idx in range(num_sample_points_y):
            sampling_radius = domain_dimensions["cornea radius"]-domain_dimensions["limbal offset"]-float(idx)*domain_dimensions["sample spacing y"]
            num_nodes = int(float((2.0*np.pi*sampling_radius)/domain_dimensions["node spacing"])) + 1
            sweep_angle = 2.0*np.pi/num_nodes
            sample_points = vtk.vtkPoints()
            dimless_sample_spacing_z = domain_dimensions["sample spacing z"]/reference_length
            for jdx in range(num_sample_points_x):
                this_angle = float(jdx)*sweep_angle+np.pi
                x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
                y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
                
                sample_points.InsertNextPoint(x_coord,y_coord, float(kdx)*dimless_sample_spacing_z)
            sample_lines.append(sample_points)
    return sample_lines

def get_3d_hemisphere_sample_points(domain_dimensions, reference_length):

    pellet_angle = np.arcsin(float(domain_dimensions["pellet height"]/domain_dimensions["cornea radius"]))
    offset_angle = np.arcsin(float(domain_dimensions["limbal offset"]/domain_dimensions["cornea radius"]))
    
    y_extent = domain_dimensions["cornea radius"]*(pellet_angle-offset_angle)
    num_sample_points_y = int(float(y_extent/domain_dimensions["sample spacing y"])) + 1

    sample_lines = []
    for idx in range(num_sample_points_y):
        current_angle = offset_angle + float(idx)*(pellet_angle-offset_angle)/float(num_sample_points_y)
        current_radius = (domain_dimensions["cornea radius"]-domain_dimensions["cornea thickness"]/2.0)*np.cos(current_angle)
        current_height = (domain_dimensions["cornea radius"]-domain_dimensions["cornea thickness"]/2.0)*np.sin(current_angle)
        num_nodes = int(float((2.0*np.pi*current_radius)/domain_dimensions["node spacing"])) + 1
        sweep_angle = 2.0*np.pi/num_nodes
        sample_points = vtk.vtkPoints()
        for jdx in range(num_nodes):
            this_angle = float(jdx)*sweep_angle+np.pi
            x_coord = (current_radius/reference_length)*np.sin(this_angle)
            y_coord = (current_radius/reference_length)*np.cos(this_angle)
            sample_points.InsertNextPoint(x_coord,y_coord, current_height/reference_length)
        sample_lines.append(sample_points)
    return sample_lines

def get_sample_points(domain_type, domain_dimensions, reference_length):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_sample_points(domain_dimensions, reference_length)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_sample_points(domain_dimensions, reference_length)
       
    elif domain_type == "Circle 2D":
        return get_2d_circle_sample_points(domain_dimensions, reference_length)
       
    elif domain_type == "Circle 3D":
        return get_3d_circle_sample_points(domain_dimensions, reference_length)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_sample_points(domain_dimensions, reference_length)

if __name__ == '__main__':
    
    import chaste.core
    from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience
    
    work_dir = "Python/Cornea/TestSamplePoints/"
    reference_length = 1.e-6 * metre()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    domain_dimensions = {"pellet height" : 1.0e-3*metre(),
                     "cornea radius" : 1.3e-3*metre(),
                     "cornea thickness" : 80.0e-6*metre(),
                     "pellet thickness" : 60.0e-6*metre(),
                     "grid spacing" : 40.0e-6*metre(),
                     "node spacing" : 40.0e-6*metre(),
                     "limbal offset" : 200.0e-6*metre(),
                     "density grid spacing" : 40.0e-6*metre(),
                     "sample spacing x" : 40.0e-6*metre(),
                     "sample spacing y" : 40.0e-6*metre(),
                     "sample spacing z" : 40.0e-6*metre(),
                     "use pellet" : False,
                     "pellet radius" : 300.0e-6*metre(),
                     }
    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    
    # Write the domain
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), True)
        sample_points = get_sample_points(eachDomainType, domain_dimensions, reference_length)
        polydata = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        for eachLine in sample_points:
            for idx in range(eachLine.GetNumberOfPoints()):
                points.InsertNextPoint(eachLine.GetPoint(idx))
        polydata.SetPoints(points)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file_handler.GetOutputDirectoryFullPath() + "sample.vtp")
        writer.SetInputData(polydata)
        writer.Write()