"""
Using the `programmable_filter` in Paraview, it's possible to convert
input axes grids, e.g. from units of meters to km
Below are some functions that can be put into the script of a programmable
filter to do these
"""

"""
CONVERT COORDINATES AND UNITS

Convert the coordinate XYZ from units of m to km and normalize origin to 0,
Change the units of the output from units of m/s to km/s
I honestly have no idea how I made this work, the documentation for these 
filters is apalling. Wtf paraview.
"""
from paraview import vtk

pdi = self.GetInput()
pdo = self.GetOutput()

numPoints = pdi.GetNumberOfPoints()
newPoints=vtk.vtkPoints()
min_x, _, min_y, _, _, _ = pdi.GetBounds()
for i in range(0, numPoints):
    coords = pdi.GetPoint(i)
    x, y, z = coords[:3]
    x -= min_x
    y -= min_y

    # Convert to units of km
    x *= 1E-3
    y *= 1E-3
    z *= 1E-3
    newPoints.InsertPoint(i, x, y, z)

# Set the new coordinate system
pdo.SetPoints(newPoints)

# Create a new array that will hold the converted values
ivals = pdi.GetPointData().GetScalars()

ca = vtk.vtkFloatArray()
ca.SetName(ivals.GetName())
ca.SetNumberOfComponents(1)
ca.SetNumberOfTuples(numPoints)

# Copy the values over element by element and convert
for i in range(0, numPoints):
  ca.SetValue(i, ivals.GetValue(i) * 1E-3)

# Set the new values
pdo.GetPointData().AddArray(ca)


"""
CONVERT COORDINATES CARRY UNITS

Convert the coordinate XYZ from units of m to km and normalize origin to 0,
carry over the unit values but don't conver those
"""
from paraview import vtk

pdi = self.GetInput()
pdo = self.GetOutput()

numPoints = pdi.GetNumberOfPoints()
newPoints=vtk.vtkPoints()
min_x, _, min_y, _, _, _ = pdi.GetBounds()
for i in range(0, numPoints):
    coords = pdi.GetPoint(i)
    x, y, z = coords[:3]
    x -= min_x
    y -= min_y

    # Convert to units of km
    x *= 1E-3
    y *= 1E-3
    z *= 1E-3
    newPoints.InsertPoint(i, x, y, z)

# Set the new coordinate system
pdo.SetPoints(newPoints)

# Create a new array that will hold the converted values
ivals = pdi.GetPointData().GetScalars()

ca = vtk.vtkFloatArray()
ca.SetName(ivals.GetName())
ca.SetNumberOfComponents(1)
ca.SetNumberOfTuples(numPoints)

# Copy the values over element by element and convert
for i in range(0, numPoints):
  ca.SetValue(i, ivals.GetValue(i))

# Set the new values
pdo.GetPointData().AddArray(ca)


"""
CONVERT COORDINATES

Convert the coordinate XYZ from units of m to km
"""
from paraview import vtk

pdi = self.GetInput()
pdo = self.GetOutput()

numPoints = pdi.GetNumberOfPoints()
newPoints=vtk.vtkPoints()
for i in range(0, numPoints):
    coords = pdi.GetPoint(i)
    x, y, z = coords[:3]
    x *= 1E-3
    y *= 1E-3
    z *= 1E-3
    newPoints.InsertPoint(i, x, y, z)
pdo.SetPoints(newPoints)

"""
CONVERT UNITS

Convert units e.g. from m/s to km/s
"""
from paraview import vtk

pdi = self.GetInput()
pdo = self.GetOutput()
ivals = pdi.GetPointData().GetScalars()
numPoints = pdi.GetNumberOfPoints()
newPoints=vtk.vtkPoints()
ca = vtk.vtkFloatArray()
ca.SetName(ivals.GetName())
ca.SetNumberOfComponents(1)
ca.SetNumberOfTuples(numPoints)

# Copy the values over element by element and convert
for i in range(0, numPoints):
  ca.SetValue(i, ivals.GetValue(i) * 1E-3)

# Set the new values
pdo.GetPointData().AddArray(ca)

"""
CONVERT UNITS

Convert units from absolute to deviation from mean
"""
from paraview import vtk

pdi = self.GetInput()
pdo = self.GetOutput()
ivals = pdi.GetPointData().GetScalars()
numPoints = pdi.GetNumberOfPoints()
newPoints=vtk.vtkPoints()
ca = vtk.vtkFloatArray()
ca.SetName(ivals.GetName())
ca.SetNumberOfComponents(1)
ca.SetNumberOfTuples(numPoints)

# Copy the values over element by element and convert
mean_ivals = 0
for i in range(0, numPoints):
    mean_ivals += ivals.GetValue(i)
mean_ivals /= numPoints

for i in range(0, numPoints):
  ca.SetValue(i, ivals.GetValue(i) - mean_ivals)

# Set the new values
pdo.GetPointData().AddArray(ca)

"""
Using the Python calculator, you can change the units of your value 
using the following command
"""
inputs[0].PointData['vs']*1E-3

"""
Convert isocontour to Z coordinates so you can color by depth
"""
inputs[0].Points[:, 2]
