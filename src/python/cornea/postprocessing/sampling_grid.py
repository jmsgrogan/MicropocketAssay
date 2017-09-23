import vtk
import numpy as np
from microvessel_chaste.utility import *


def GetDensityMetrics(file_path, domain, pc):

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.Update()

    grid = reader.GetOutput()
    line_density = grid.GetCellData().GetArray("Line Density")

    line_density_volume = 0.0
    total_volume = 0.0
    offset = pc.get_parameter("LimbalOffset").value
    offset = offset.Convert(1.0e-6*metres)
    height = pc.get_parameter("PelletHeight").value
    height = height.Convert(1.0e-6*metres)

    for idx in range(grid.GetNumberOfCells()):
        points = grid.GetCell(idx).GetPoints()
        loc = np.zeros(3)
        for jdx in range(points.GetNumberOfPoints()):
            if jdx == 0:
                loc = np.array(points.GetPoint(jdx))
            else:
                loc += np.array(points.GetPoint(jdx))
        loc /= float(points.GetNumberOfPoints())
        cell_volume = abs(vtk.vtkMeshQuality.HexVolume(grid.GetCell(idx)))
        total_volume += cell_volume
        cell_line_density = line_density.GetTuple1(idx)
        if cell_line_density >= 1.e-9:
            line_density_volume += cell_volume
    line_fraction = line_density_volume/total_volume
    return [line_fraction, 1.0]


def DoLineSampling(file_path, domain, pc):

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.Update()

    sample_spacing = 60.0
    height = pc.get_parameter("PelletHeight").value
    height = height.Convert(1.0e-6*metres)
    radius = pc.get_parameter("CorneaRadius").value
    radius = radius.Convert(1.0e-6*metres)
    width = 2.0*np.pi*radius

    num_samples = int(height/sample_spacing)
    points = vtk.vtkPoints()
    for idx in range(num_samples):
        points.InsertNextPoint(width/2.0, float(idx*sample_spacing), 0.0)
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
        
    probe = vtk.vtkProbeFilter()
    probe.SetSourceData(reader.GetOutput())
    probe.SetInputData(poly)
    probe.Update()
    
    results = probe.GetOutput().GetPointData().GetArray("Line Density")
    for idx in range(results.GetNumberOfTuples()):
        print "y", points.GetPoint(idx)[1], " rho ", results.GetTuple1(idx)
