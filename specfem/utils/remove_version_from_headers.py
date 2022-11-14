"""
Removes Version numbers from SPECFEM source code. Goes through the entire
source code and drops the version number from the header.

Super hacky to work with multiple languages and header formats
"""
import os
from glob import glob


DRYRUN = False

def drop_version(fid):
    """drops version numbers, makes some assumptions about how the header is
    formatted"""
    changed = False
    with open(fid, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines[:]):
        if "S P E C F E M 2 D" in line.upper():
            # Some files have a space at the start of the header
            if line.startswith(" "):
                space = " "
            else:
                space = ""
            if "!" in line:
                comment = "!"
            elif "#" in line:
                comment = "#"
            else:
                print(f"comment? {line}")
                a = 1/0
            lines[i] =   f"{space}{comment}                          S p e c f e m 2 D\n"
            lines[i+1] = f"{space}{comment}                          -----------------\n"
            changed = True
            break
    if changed:
        if not DRYRUN:
            with open((fid), "w") as f:
                f.writelines(lines)
        print(f"{fid} changed")
    else:
        print(f"{fid} NOT changed")



if __name__ == "__main__":
    repo = "./"

    i = 0
    for dir_ in glob(os.path.join(repo, "*")):
        for fid in glob(os.path.join(dir_, "*")):
            # assuming only one subdirectory
            if os.path.isdir(fid):
                for fid_ in glob(os.path.join(fid, "*")):
                    if os.path.isdir(fid_):
                        continue
                    else:
                        drop_version(fid_)
            else:
                drop_version(fid)
                i += 1



