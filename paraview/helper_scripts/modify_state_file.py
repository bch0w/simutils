"""
Used to modify an existing ParaView State File (.pvsm) to change input file name
and colorbar scale
"""
import os


def replace_parts(file, file_out, replacements):
    """Replace parts of a file with new values"""
    # Read in original file
    with open(file, "r") as f:
        lines = f.readlines()

    # Replace lines in file
    for i, line in enumerate(lines[:]):
        # Replace file name and tag
        for key, val in replacements.items():
            if key in line:
                lines[i] = line.replace(str(key), str(val))

    with open(file_out, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    # Explanation of `replacements`
    # FILENAME: name of the .vtk file
    # PARAMETER: name of the parameter inside the .vtk file (e.g., alpha_kernel)
    # COLORBAR_LABEL: how you want the colorbar to be tagged
    # SCALE: the max value of the colorbar 
    state_file = "kernel_slice_xz_plane.pvsm"
    template_file = "template_kernel_slice_xz_plane.pvsm"

    # Create a Template file
    if not os.path.exists(template_file):
        print("creating template file")
    
        # Template will have DUMMY values to make it easier to replace things
        # These original values come from the manual creation of the original
        # state file
        replacements = {"S_beta.vtk": "FILENAME",
                        "beta_kernel": "PARAMETER",
                        "S Beta Kernel": "COLORBAR_LABEL",
                        "1.5e-12": "SCALE",
                        }
        replace_parts(state_file, template_file, replacements)

    # Modify for P-Alpha kernel
    print("modifying for p-alpha kernel")
    file_out = "p_alpha_kernel_slice_xz_plane.pvsm"
    replacements = {"FILENAME": "P_alpha.vtk",
                    "PARAMETER": "alpha_kernel",
                    "COLORBAR_LABEL": "P Alpha Kernel",
                    "SCALE": 0.9E-13
                    }
    replace_parts(template_file, file_out, replacements)

    # Modify for SS-Beta kernel
    print("modifying for ss-beta kernel")
    file_out = "ss_beta_kernel_slice_xz_plane.pvsm"
    replacements = {"FILENAME": "SS_beta.vtk",
                    "PARAMETER": "beta_kernel",
                    "COLORBAR_LABEL": "SS Beta Kernel",
                    "SCALE": 1.5E-12
                    }
    replace_parts(template_file, file_out, replacements)

                

