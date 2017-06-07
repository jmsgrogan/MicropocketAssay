import numpy as np
import vtk

def get_2d_planar_sample_points(parameter_collection, reference_length):

    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    sample_spacing_x = parameter_collection.get_parameter("sample spacing x").value
    sample_spacing_y = parameter_collection.get_parameter("sample spacing y").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    domain_width = 2.0*np.pi*cornea_radius
    num_sample_points_x = int(float(domain_width/sample_spacing_x)) + 1
    num_sample_points_y = int(float((pellet_height-limbal_offset)/sample_spacing_y)) + 1

    sample_lines = []
    for idx in range(num_sample_points_y):
        sample_points = vtk.vtkPoints()
        dimless_sample_spacing_x = sample_spacing_x/reference_length
        dimless_sample_spacing_y = sample_spacing_y/reference_length
        dimless_limbal_offset = limbal_offset/reference_length
        for jdx in range(num_sample_points_x):
            sample_points.InsertNextPoint(float(jdx)*dimless_sample_spacing_x,
                                          float(idx)*dimless_sample_spacing_y + dimless_limbal_offset,
                                          0.0)
        sample_lines.append(sample_points)
    return [sample_lines, num_sample_points_y]

def get_3d_planar_sample_points(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    conrea_thickness = parameter_collection.get_parameter("cornea thickness").value
    sample_spacing_x = parameter_collection.get_parameter("sample spacing x").value
    sample_spacing_y = parameter_collection.get_parameter("sample spacing y").value
    sample_spacing_z = parameter_collection.get_parameter("sample spacing z").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    domain_width = 2.0*np.pi*cornea_radius
    num_sample_points_x = int(float(domain_width/sample_spacing_x)) + 1
    num_sample_points_y = int(float((pellet_height-limbal_offset)/sample_spacing_y)) + 1
    num_sample_points_z = int(float(conrea_thickness/sample_spacing_z)) + 1
    
    sample_lines = []
    for kdx in range(num_sample_points_z):
        for idx in range(num_sample_points_y):
            sample_points = vtk.vtkPoints()
            dimless_sample_spacing_x = sample_spacing_x/reference_length
            dimless_sample_spacing_y = sample_spacing_y/reference_length
            dimless_sample_spacing_z = sample_spacing_z/reference_length
            dimless_limbal_offset = limbal_offset/reference_length
            for jdx in range(num_sample_points_x):
                sample_points.InsertNextPoint(float(jdx)*dimless_sample_spacing_x,
                                              float(idx)*dimless_sample_spacing_y + dimless_limbal_offset,
                                              float(kdx)*dimless_sample_spacing_z)
            sample_lines.append(sample_points)
    return [sample_lines, num_sample_points_y]

def get_2d_circle_sample_points(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    node_spacing = parameter_collection.get_parameter("node spacing").value
    sample_spacing_x = parameter_collection.get_parameter("sample spacing x").value
    sample_spacing_y = parameter_collection.get_parameter("sample spacing y").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    domain_width = 2.0*np.pi*cornea_radius
    num_sample_points_x = int(float(domain_width/sample_spacing_x)) + 1
    num_sample_points_y = int(float(pellet_height/sample_spacing_y)) + 1

    sample_lines = []

    for idx in range(num_sample_points_y):
        sampling_radius = cornea_radius-limbal_offset-float(idx)*sample_spacing_y
        num_nodes = int(float((2.0*np.pi*sampling_radius)/node_spacing)) + 1
        sweep_angle = 2.0*np.pi/num_nodes
        sample_points = vtk.vtkPoints()
        for jdx in range(num_sample_points_x):
            this_angle = float(jdx)*sweep_angle+np.pi
            x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
            y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
            sample_points.InsertNextPoint(x_coord,y_coord, 0.0)
        sample_lines.append(sample_points)
    return [sample_lines, num_sample_points_y]

def get_3d_circle_sample_points(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    conrea_thickness = parameter_collection.get_parameter("cornea thickness").value
    sample_spacing_x = parameter_collection.get_parameter("sample spacing x").value
    sample_spacing_y = parameter_collection.get_parameter("sample spacing y").value
    sample_spacing_z = parameter_collection.get_parameter("sample spacing z").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value
    node_spacing = parameter_collection.get_parameter("node spacing").value

    domain_width = 2.0*np.pi*cornea_radius
    num_sample_points_x = int(float(domain_width/sample_spacing_x)) + 1
    num_sample_points_y = int(float(pellet_height/sample_spacing_y)) + 1
    num_sample_points_z = int(float(conrea_thickness/sample_spacing_z)) + 1
    sample_lines = []

    for kdx in range(num_sample_points_z):
        for idx in range(num_sample_points_y):
            sampling_radius = cornea_radius-limbal_offset-float(idx)*sample_spacing_y
            num_nodes = int(float((2.0*np.pi*sampling_radius)/node_spacing)) + 1
            sweep_angle = 2.0*np.pi/num_nodes
            sample_points = vtk.vtkPoints()
            dimless_sample_spacing_z = sample_spacing_z/reference_length
            for jdx in range(num_sample_points_x):
                this_angle = float(jdx)*sweep_angle+np.pi
                x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
                y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
                sample_points.InsertNextPoint(x_coord,y_coord, float(kdx)*dimless_sample_spacing_z)
            sample_lines.append(sample_points)
    return [sample_lines, num_sample_points_y] 

def get_3d_hemisphere_sample_points(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    pellet_height = parameter_collection.get_parameter("pellet height").value
    conrea_thickness = parameter_collection.get_parameter("cornea thickness").value
    sample_spacing_y = parameter_collection.get_parameter("sample spacing y").value
    sample_spacing_z = parameter_collection.get_parameter("sample spacing z").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value
    node_spacing = parameter_collection.get_parameter("node spacing").value

    pellet_angle = np.arcsin(float(pellet_height/cornea_radius))
    offset_angle = np.arcsin(float(limbal_offset/cornea_radius))
    
    y_extent = cornea_radius*(pellet_angle-offset_angle)
    num_sample_points_y = int(float(y_extent/sample_spacing_y)) + 1
    num_sample_points_z = int(float(conrea_thickness/sample_spacing_z)) + 1

    sample_lines = []
    for kdx in range(num_sample_points_z):
        for idx in range(num_sample_points_y):
            current_angle = offset_angle + float(idx)*(pellet_angle-offset_angle)/float(num_sample_points_y)
            sample_offset_from_outside = float(kdx)*sample_spacing_z
            current_radius = (cornea_radius-sample_offset_from_outside)*np.cos(current_angle)
            current_height = (cornea_radius-sample_offset_from_outside)*np.sin(current_angle)
            num_nodes = int(float((2.0*np.pi*current_radius)/node_spacing)) + 1
            sweep_angle = 2.0*np.pi/num_nodes
            sample_points = vtk.vtkPoints()
            for jdx in range(num_nodes):
                this_angle = float(jdx)*sweep_angle+np.pi
                x_coord = (current_radius/reference_length)*np.sin(this_angle)
                y_coord = (current_radius/reference_length)*np.cos(this_angle)
                sample_points.InsertNextPoint(x_coord,y_coord, current_height/reference_length)
            sample_lines.append(sample_points)
    return [sample_lines, num_sample_points_y]

def get_sample_points(domain_type, parameter_collection, reference_length):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_sample_points(parameter_collection, reference_length)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_sample_points(parameter_collection, reference_length)
       
    elif domain_type == "Circle 2D":
        return get_2d_circle_sample_points(parameter_collection, reference_length)
       
    elif domain_type == "Circle 3D":
        return get_3d_circle_sample_points(parameter_collection, reference_length)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_sample_points(parameter_collection, reference_length)

if __name__ == '__main__':
    
    import chaste.core
    import cornea.parameters.default_parameters
    from microvessel_chaste.utility import * # Dimensional analysis: bring in all units for convenience
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    
    work_dir = "Python/Cornea/TestSamplePoints/"
    reference_length = 1.e-6 * metre()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    
    # Write the domain
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), True)
        sample_lines, num_sample_points_y = get_sample_points(eachDomainType, parameter_collection, reference_length)
        print eachDomainType, " num_z ", len(sample_lines)/num_sample_points_y
        polydata = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        for eachLine in sample_lines:
            for idx in range(eachLine.GetNumberOfPoints()):
                points.InsertNextPoint(eachLine.GetPoint(idx))
        polydata.SetPoints(points)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file_handler.GetOutputDirectoryFullPath() + "sample.vtp")
        writer.SetInputData(polydata)
        writer.Write()