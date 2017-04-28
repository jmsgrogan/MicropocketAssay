import vtk

array1 = vtk.vtkDoubleArray()
array1.SetNumberOfTuples(2)
array1.SetTuple1(0, 1.0)
array1.SetTuple1(1, 2.0)

holder = vtk.vtkPointData()
holder.AddArray(array1)