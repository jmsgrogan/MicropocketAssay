"""
A collection of simulation vessel networks
"""

import numpy as np
import microvessel_chaste.mesh
import microvessel_chaste.population.vessel
from microvessel_chaste.utility import *

def get_2d_planar_vessel_network(parameter_collection, reference_length):

    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    node_spacing = parameter_collection.get_parameter("node spacing").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value
    
    generator = microvessel_chaste.population.vessel.VesselNetworkGenerator2()
    domain_length = 2.0*np.pi*cornea_radius
    divisions = int(float(domain_length/node_spacing)) - 2
    alignment_axis = 0 # pointing x direction
    network = generator.GenerateSingleVessel(domain_length, 
                                             microvessel_chaste.mesh.DimensionalChastePoint2(0.0, 
                                                                     float(limbal_offset/reference_length), 
                                                                     0.0, reference_length), 
                                             divisions, alignment_axis)
    for eachNode in network.GetNodes():
        eachNode.GetFlowProperties().SetPressure(1.0*pascal())
    return network

def get_3d_planar_vessel_network(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    cornea_thickness = parameter_collection.get_parameter("cornea thickness").value
    node_spacing = parameter_collection.get_parameter("node spacing").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    generator = microvessel_chaste.population.vessel.VesselNetworkGenerator3()
    domain_length = 2.0*np.pi*cornea_radius
    divisions = int(float(domain_length/node_spacing)) - 2
    alignment_axis = 0 # pointing x direction
    
    mid_point = cornea_thickness/(2.0*reference_length)
    network = generator.GenerateSingleVessel(domain_length, 
                                             microvessel_chaste.mesh.DimensionalChastePoint3(0.0, 
                                                                     float(limbal_offset/reference_length), 
                                                                     mid_point, reference_length), 
                                             divisions, alignment_axis)
    for eachNode in network.GetNodes():
        eachNode.GetFlowProperties().SetPressure(1.0*pascal())
    return network

def get_2d_circle_vessel_network(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    node_spacing = parameter_collection.get_parameter("node spacing").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    sampling_radius = cornea_radius-limbal_offset
    num_nodes = int(float((2.0*np.pi*sampling_radius)/node_spacing)) + 1
    sweep_angle = 2.0*np.pi/float(num_nodes)

    network = microvessel_chaste.population.vessel.VesselNetwork2.Create()
    nodes = []
    for idx in range(num_nodes):
        this_angle = float(idx)*sweep_angle+np.pi
        x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
        y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
        nodes.append(microvessel_chaste.population.vessel.VesselNode2.Create(
            microvessel_chaste.mesh.DimensionalChastePoint2(x_coord, y_coord, 0.0, reference_length)))

    for idx in range(len(nodes)):
        network.AddVessel(microvessel_chaste.population.vessel.Vessel2.Create(nodes[idx-1], nodes[idx]))
    network.AddVessel(microvessel_chaste.population.vessel.Vessel2.Create(nodes[len(nodes)-1], nodes[0]))
    
    for idx in range(len(nodes)):
        network.GetNodes()[idx].GetFlowProperties().SetPressure(1.0*pascal())        
    return network

def get_3d_circle_vessel_network(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    cornea_thickness = parameter_collection.get_parameter("cornea thickness").value
    node_spacing = parameter_collection.get_parameter("node spacing").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    sampling_radius = cornea_radius-limbal_offset
    num_nodes = int(float((2.0*np.pi*sampling_radius)/node_spacing)) + 1
    sweep_angle = 2.0*np.pi/float(num_nodes)
    mid_point = cornea_thickness/(2.0*reference_length)

    network = microvessel_chaste.population.vessel.VesselNetwork3.Create()
    nodes = []
    for idx in range(num_nodes):
        this_angle = float(idx)*sweep_angle+np.pi
        x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
        y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
        nodes.append(microvessel_chaste.population.vessel.VesselNode3.Create(
            microvessel_chaste.mesh.DimensionalChastePoint3(x_coord, y_coord, mid_point, reference_length)))

    for idx in range(len(nodes)):
        network.AddVessel(microvessel_chaste.population.vessel.Vessel3.Create(nodes[idx-1], nodes[idx]))
    network.AddVessel(microvessel_chaste.population.vessel.Vessel3.Create(nodes[len(nodes)-1], nodes[0]))
    
    for idx in range(len(nodes)):
        network.GetNodes()[idx].GetFlowProperties().SetPressure(1.0*pascal())        
    return network

def get_3d_hemisphere_vessel_network(parameter_collection, reference_length):
    
    cornea_radius = parameter_collection.get_parameter("cornea radius").value
    cornea_thickness = parameter_collection.get_parameter("cornea thickness").value
    node_spacing = parameter_collection.get_parameter("node spacing").value
    limbal_offset = parameter_collection.get_parameter("limbal offset").value

    dimless_cornea_mid_radius = float((cornea_radius-cornea_thickness/2.0)/reference_length)
    dimless_offset = float(limbal_offset/reference_length)
    sampling_radius = np.sqrt(dimless_cornea_mid_radius**2-dimless_offset**2)*reference_length
    num_nodes = int(float((2.0*np.pi*sampling_radius)/node_spacing)) + 1
    sweep_angle = 2.0*np.pi/float(num_nodes)

    network = microvessel_chaste.population.vessel.VesselNetwork3.Create()
    nodes = []
    for idx in range(num_nodes):
        this_angle = float(idx)*sweep_angle+np.pi
        x_coord = (sampling_radius/reference_length)*np.sin(this_angle)
        y_coord = (sampling_radius/reference_length)*np.cos(this_angle)
        z_coord = dimless_offset
        nodes.append(microvessel_chaste.population.vessel.VesselNode3.Create(
            microvessel_chaste.mesh.DimensionalChastePoint3(x_coord, y_coord, z_coord, reference_length)))    
    
    for idx in range(len(nodes)):
        network.AddVessel(microvessel_chaste.population.vessel.Vessel3.Create(nodes[idx-1], nodes[idx]))
    network.AddVessel(microvessel_chaste.population.vessel.Vessel3.Create(nodes[len(nodes)-1], nodes[0]))    
    
    for idx in range(len(nodes)):
        network.GetNodes()[idx].GetFlowProperties().SetPressure(1.0*pascal())    
    return network

def get_vessel_network(domain_type, parameter_collection, reference_length):
    
    if domain_type == "Planar 2D":
        return get_2d_planar_vessel_network(parameter_collection, reference_length)
    
    elif domain_type == "Planar 3D":
        return get_3d_planar_vessel_network(parameter_collection, reference_length)
     
    elif domain_type == "Circle 2D":
        return get_2d_circle_vessel_network(parameter_collection, reference_length)
     
    elif domain_type == "Circle 3D":
        return get_3d_circle_vessel_network(parameter_collection, reference_length)
    
    elif domain_type == "Hemisphere 3D":
        return get_3d_hemisphere_vessel_network(parameter_collection, reference_length)

if __name__ == '__main__':
    
    import chaste.core
    import cornea.parameters.default_parameters
    
    work_dir = "Python/Cornea/TestNetworks/"
    reference_length = 1.e-6 * metre()
    BaseUnits.Instance().SetReferenceLengthScale(reference_length)
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
        
    domain_types = ["Planar 2D", "Planar 3D", "Circle 2D", "Circle 3D", "Hemisphere 3D"]
    
    # Write the domain
    for eachDomainType in domain_types:
        file_handler = chaste.core.OutputFileHandler(work_dir + "/" + eachDomainType.replace(" ", ""), True)
        network = get_vessel_network(eachDomainType, parameter_collection, reference_length)
        network.Write(file_handler.GetOutputDirectoryFullPath() + "network.vtp")