"""
To be run in the Paraview Python shell to automatically generate a model
visualization with accoutrements and take screenshots. 
"""
import os

# Here we define some preset parameters that should stay constant among certain
# file types
PRESETS = {
        "model_vp": {"cbar_title": "Vp [m/s]",
                     "colormap": "Rainbow Desaturated",
                     "invert_cmap": True,
                     },
        "model_vs": {"cbar_title": "Vs [m/s]",
                     "colormap": "Rainbow Desaturated",
                     "invert_cmap": True,
                     },
        "gradient_vp_kernel": {"cbar_title": "Vp Gradient [m^-2 s^2]",
                               "colormap": "Cool to Warm (Extended)",
                               "invert_cmap": False,
                               },

        "gradient_vs_kernel": {"cbar_title": "Vs Gradient [m^-2 s^2]",
                               "colormap": "Cool to Warm (Extended)",
                               "invert_cmap": False,
                               },
        "update_vp": {"cbar_title": "Vp Net Update [%]",
                      "colormap": "Blue Orange (divergent)",
                      "invert_cmap": False,
                      },
        "update_vs": {"cbar_title": "Vs Net Update [%]",
                      "colormap": "Blue Orange (divergent)",
                      "invert_cmap": False,
                      },
        "ratio_poissons": {"cbar_title": "Poisson's Ratio",
                           "colormap": "Blue - Green - Orange",
                           "invert_cmap": False,
                           },
        "ratio_vpvs": {"cbar_title": "Vp/Vs Ratio",
                       "colormap": "",
                       "invert_cmap": False,
                       },
        }

# File ID's specifying the chosen files
state_fid = "/Users/Chow/Documents/academic/vuw/forest/forest.pvsm"
data_fid = ("/Users/Chow/Documents/academic/vuw/forest/birch/vtk/model/"
            "update_0005_vs.vtk")
slices = [2, 4, 6, 8, 10, 15, 20, 25, 30, 50]

# Presets defined by labels assigned to the file names, ignoring number label
preset = os.path.splitext(os.path.basename(data_fid))[0].split("_")
preset = "_".join(preset[:1] + preset[2:])
assert(preset in PRESETS), f"{preset} does not match preset keys"
preset = PRESETS[preset]

# These values are defined by presets or by the user
colormap = preset["colormap"]
cbar_title = preset["cbar_title"]
invert_cmap = preset["invert_cmap"]

# Load an active state that has common traits, e.g. src-rcv glyphs, coast
servermanager.LoadState(state_fid)

# Plot surface trace
vtk = OpenDataFile(data_fid)
quantity = vtk.PointData.GetArray(0)  # e.g. model_init_vp
RenameSource("surface", vtk)
Show(vtk)
ResetCamera()

# Create annotation
renderView = GetActiveViewOrCreate("RenderView")
text = Text()
text.Text = "surface"
textDisplay = Show(text, renderView)
textDisplay.Position = [0.525, 0.05]
Render()

def slice_plane(depth):
    """
    Slice plane through the volume, defined mid script to use defined variables
    """
    renderView = GetActiveViewOrCreate("RenderView")
    slice_vtk = Slice(vtk)
    slice_vtk.SliceType = "Plane"
    slice_vtk.SliceType.Normal = [0., 0., 1.]
    origin = slice_vtk.SliceType.Origin
    origin = [origin[0], origin[1], abs(depth) * -1E3]
    slice_vtk.SliceType.Origin = origin
    sliceDisplay = Show(slice_vtk, renderView)
    RenameSource(f"{depth}km", slice_vtk)
    Hide(slice_vtk, renderView)

    return slice_vtk

# Move camera around
camera = GetActiveCamera()
camera.Dolly(1.25)  # Move closer

# Change the colormap
vsLUT = GetColorTransferFunction(quantity)
vsLUT.ApplyPreset(colormap, True)
vsLUT.RescaleTransferFunction(vmin, vmax)
if invert_cmap:
    vsLUT.InvertTransferFunction()

# Create the colorbar
cbar = GetScalarBar(vsLUT)
cbar.Title = cbar_title
cbar.ComponentTitle = ""
cbar.AutoOrient = 0
cbar.Orientation = "Vertical"
cbar.Position = [0.3523, 0.7522]

planes = []
for slice_ in slices:
    planes.append(slice_plane(slice_))


