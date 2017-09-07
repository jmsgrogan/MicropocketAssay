import vtk
import numpy as np

path = "/home/grogan/test.vtu"

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(path)
reader.Update()

grid = reader.GetOutput()
line_density = grid.GetCellData().GetArray("Line Density")
tip_density = grid.GetCellData().GetArray("Tip Density")

line_density_volume = 0.0
tip_density_volume = 0.0
total_volume = 0.0

for idx in range(grid.GetNumberOfCells()):
    points = grid.GetCell(idx).GetPoints()
    loc = np.zeros(3)
    for jdx in range(points.GetNumberOfPoints()):
        if jdx == 0:
            loc = np.array(points.GetPoint(jdx))
        else:
            loc += np.array(points.GetPoint(jdx))
    loc/=float(points.GetNumberOfPoints())
    cell_volume = abs(vtk.vtkMeshQuality.HexVolume(grid.GetCell(idx)))
    total_volume += cell_volume
    cell_line_density = line_density.GetTuple1(idx)
    if cell_line_density>=1.e-9:
        line_density_volume += cell_volume
    cell_tip_density = tip_density.GetTuple1(idx)
    if cell_tip_density>=1.e-9:
        tip_density_volume += cell_volume
        
print line_density_volume/total_volume
print tip_density_volume/total_volume