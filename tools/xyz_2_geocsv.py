"""
Convert the .xyz files required for SPECFEM3D's external tomography file to
a GeoCSV file which contains more human-readable information. The GeoCSV file
can then be converted directly to netCDF using the IRIS emc-tools.
Also can handle multiple 

NOTE:
    iris emc-tools/GeoCSV_2_netCDF_3D.py 
        requires 3 header columns: 
        latitude_column, longitude_column, depth_column
        latitude can also be x or X
        longitude can also be y or Y
        X depth must be depth because part of their code fails if it's z
        header attributes such as 'author_name' must start with 'global', i.e.,
            # global_author_name: Bryant Chow

Relevant links:
    https://github.com/iris-edu/emc-tools
    http://geows.ds.iris.edu/documents/GeoCSV.pdf
"""
import os
import sys
import datetime
import numpy as np


class Field(dict):
    """
    GeoCSV requires the field to have a few parameters, this class is a simple
    dictionary that ensures all these attributes are met. Also allows expansion
    to new attributes that are automatically filled in as the Field is parsed.
    """
    def __init__(self, std_name, long_name, unit, grid=False):
        """
        Set field attributes required for GeoCSV file
        
        :type std_name: str
        :param std_name: standard or short name for the variable, e.g. vp
        :type long_name: str
        :param long_name: long name for the variable, e.g. P-velocity
        :type unit: str
        :param unit: standard unit for the variable, e.g. m/s
        :type grid: bool
        :param grid: if True, a new header line will be added describing the
            discretization of the variable. Useful only for uniformly gridded
            variables like the coordinate system, where x is divided into equal
            dx values. Defaults to False
        """
        self.std_name = std_name
        self.long_name = long_name
        self.unit = unit
        self.grid = grid

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class Converter():
    """
    Class to control the conversion from XYZ to GeoCSV
    """
    def __init__(self, fields=None, delimiter=",", header_lines=4, fmt="%.1f",
                 path_out=os.getcwd(), delimiter_in=None):
        """
        Set some constant parameters that the converter uses for internal checks

        :type fields: list of Fields
        :param fields: the different data types that are included in the XYZ
        :type delimiter: str
        :param delimiter: how to separate variables when writing to a new file
            see `delimiter_in` if the input and output files should have 
            different delimiters
        :type header_lines: int
        :param header_lines: tells the Converter how many lines to skip in the 
            input .xyz file so that it only reads data
        :type fmt: str or list of str
        :param fmt: string formatter for writing data to GeoCSV file. Piped into
            the numpy.savetxt() function so can be one blanket formatter or
            individual formats for each variable
        :type path_out: str
        :param path_out: path to save the output .csv files. If None, set to
            current working directory
        :type delimiter_in: str
        :param delimiter_in: If the input file has a different delimiter than
            the output file, set here. If set to NoneType, will be the same as
            `delimiter`
        """
        print("Initiating Converter for .xyz to GeoCSV")
        self.fields = fields or []
        self.delimiter = delimiter
        self.header_lines = header_lines
        self.fmt = fmt
        self.fstr = fmt[1:]  # for f-string formatting, strip '%'
        self.path_out = path_out
        self.delimiter_in = delimiter_in or delimiter

        # Initiate empty variables to be filled by other functions
        self.fid = None
        self.fid_out = None
        self.f_out = None
        self.data = None

    def append(self, **kwargs):
        """
        Append new Fields to the internal field list

        .. note::
            kwargs should match the required attributes of the Field class
        """
        self.fields.append(Field(**kwargs))

    def set_fids(self, fid, fid_out=None):
        """
        Point the converter to the input .xyz file, it will then auto set an
        output fid for the GeoCSV or the user can set it themselves.

        :type fid: str
        :param fid: path and name of the input .xyz file
        :type fid_out: str
        :param fid_out: optional path to the output .xyz file. If None then will
            be set automatically based on 'fid'
        """
        if os.path.exists(fid):
            self.fid = fid
            # Auto set the output fid to the be the same thing but with .csv
            if fid_out is None:
                fid_ = os.path.basename(fid)
                assert(".xyz" in fid_), f"input file requires .xyz  extension"
                fid_new = fid_.replace(".xyz", ".csv")
                fid_out = os.path.join(self.path_out, fid_new)

            self.fid_out = fid_out
            self.f_out = open(self.fid_out, "w")
        else:
            sys.exit(f"{fid} does n1Got exist and must")

        print(f"\tfile to read: {fid}")
        print(f"\toutput file will be: {fid_out}")

    def read_xyz(self, fid=None, swap_idx=None):
        """
        Read the .xyz file and store the values as a dict object

        .. warning::
            This does not talk to the append() function, which sets the names
            for each of the data columns. You must be sure that if you're 
            swapping indices, that you set the fields according to the newly
            swapped fields

        :type fid: str
        :param fid: read the input xyz file with numpy loadtxt
        :type swap_idx: list of lists
        :param swap_idx: swap two indices in list order. used because IRIS wants
            the y-column first (latitude) but my data is formatted in x, y, ...
            So I include swap_idx=[[0, 1]], to swap the first two indices
            and get y first. 
        """
        print(f"\treading input .xyz file")
        if fid is None:
            fid = self.fid
        self.data = np.loadtxt(fid, skiprows=self.header_lines)
        if swap_idx:
            for swap in swap_idx:
                print(f"\tswapping data arrays {swap[0]} and {swap[1]}")
                self.data[:, swap] = self.data[:, swap[::-1]]


    def convert_data(self, convert):
        """
        IRIS EMC standard format is to have units in 'km' or 'km/s' for velocity
        and the have depth positive down (i.e., deeper depths are larger,
        positive values). NZAtom was originally formatted in units of 'm' and 
        'm/s' and with depth positive up. Simply scale and flip some signs to 
        fix.

        :type convert: dict
        :param convert: dictionary that tells this function how to convert stuff
            keys should match `field.std_name`, corresponding value should be 
            an int or float by which the data is multiplied. 

            for example
            convert = {"vp": 1E-3}  # multiplies vp by 1E-3, for m/s -> km/s
        """
        print(f"\tconverting data for {len(convert)} data columns")

        # Match data indices with value names
        names = [_.std_name for _ in self.fields]

        for name, cnvt in convert.items():
            idx = names.index(name)
            og_val = self.data[0, idx]
            self.data[:, idx] *= cnvt

            # Print a check statement to make sure this is working
            print(f"\t\t{name}: {og_val} -> {self.data[0, idx]}")

    def convert_coordinates(self, epsg_in, epsg_out, x_in, y_in, 
                            choice="append", insert=None, replace=None):
        """
        Convert coordinates from one coordinate system to something different.
        Useful for NZAtom as people probably want Lat/Lon but the model
        is it UMT60S Cartesian

        .. rubric::
            # Simply append the new coordinates the end
            Converter.convert_coordinates(epsg_in, epsg_out, x_in, x_out, 
                                          choice=="append")
            # Insert the coordinates into columns 3 and 4 of the data array
            Converter.convert_coordinates(epsg_in, epsg_out, x_in, x_out, 
                                          choice=="insert", insert=3)
            # Replace columns 2 and 3 with the new coordinates
            Converter.convert_coordinates(epsg_in, epsg_out, x_in, x_out, 
                                          choice=="replace", replace=(3,4))
        :type choice: str
        :param choice: 'append': place new coords at end of fields, 
            'insert': insert new coords at location insert
            'replace': replace given tuple of locations 
        :type insert: int
        :param insert: place in the fields list to insert new coordinates
        :type replace: tuple of int
        :param replace: locations of the field data to be replaced
        """
        print(f"\tconverting coordinates from EPSG {epsg_in} to {epsg_out}")
        print(f"\tConverter will {choice} new coordinates to data array")
        if self.data is None:
            self.read_xyz()

        from pyproj import CRS, Transformer

        crs_in = CRS.from_epsg(epsg_in)
        crs_out = CRS.from_epsg(epsg_out)
        transform = Transformer.from_crs(crs_in, crs_out)

        x_out, y_out = transform.transform(x_in, y_in)
        # Choose how these new coordinates are inserted into the fields argument
        if choice == "append":
            self.data.append(x_out)
            self.data.append(y_out)
        elif choice == "insert":
            self.data = np.insert(self.data, insert, x_out, axis=1)
            self.data = np.insert(self.data, insert + 1, y_out, axis=1)
        elif choice == "replace":
            self.data[replace[0]] = x_out
            self.data[replace[1]] = y_out

    def write_header(self, f=None, prepend=None, append=None):
        """
        Write the header of the output file using field information. The IRIS
        emc tools requires a specific header structure that it uses to parse 
        data itself so this function tries to mimic that

        .. note::
            Header lines MUST be formatted with leading '#' and trailing '\n'
            e.g. 
            # This is a header line, blah blah blah \n

        :type f: _io.TextIOWrapper
        :param f: optional open text file, if None then set automatically
        :type prepend: dict
        :param prepend: optional lines to write after the first header
            declaration. Will be formatted with a leading '#' and a trailing
            '\n' to match the header format
        :type append: dict
        :param append: optional lines to write after the final standard header
            declaration. Will be formatted with a leading '#' and a trailing
            '\n' to match the header format
        """
        print("\twriting header information")
        if f is None:
            f = self.f_out

        f.write("# dataset: GeoCSV2.0\n")

        # Write some standard header information
        f.write(f"# created: {datetime.datetime.utcnow()} UTC\n")
        nc_file = os.path.basename(f.name).replace("csv", "nc")
        f.write(f"# netCDF_file: {nc_file} \n")
        f.write(f"# delimiter: {self.delimiter}\n")
        if prepend is not None:
            for key, val in prepend.items():
                f.write(f"# {key}: {val}\n")


        assert(len(self.fields) == len(self.data[0])), \
                "Data and field length mismatch, maybe convert coords first"

        # Write geospatial data as global header information, required by IRIS
        # i.e., '# global_geospatial_lat_min: ...'
        cmt = "global_geospatial"
        for name in ["latitude", "longitude", "depth"]:
            for i, field in enumerate(self.fields):
                if field["long_name"] == name:
                    if name == "depth":
                        name = "vertical"
                    else:
                        name = field.std_name
                    f.write(f"# {cmt}_{name}_min: "
                            f"{self.data[:, i].min(): {self.fstr}}\n")
                    f.write(f"# {cmt}_{name}_max: "
                            f"{self.data[:, i].max(): {self.fstr}}\n")
                    f.write(f"# {cmt}_{name}_units: {field.unit}\n")
        
                    # Need to define sign of the vertical axis, we're assuming
                    # that the model is much larger in one direction than the
                    # other and using that to determine which way is up 
                    if name == "vertical":
                        vert_data = self.data[:, i]
                        if abs(vert_data.min()) > abs(vert_data.max()):
                            vert_positive = "up"
                        else:
                            vert_positive = "down"
                        f.write(f"# {cmt}_vertical_positive: {vert_positive}\n")

        # Write header information REQUIRED by GeoCSV and IRIS emc-tools
        for i, field in enumerate(self.fields):
            name = field.std_name
            # Spatial dimensions are independent. All other params are dep on 
            # the spatial dimensions x, y, z
            if name in ["x", "y", "z", "lat", "lon"]:
                dim = 1
            else:
                dim = 3

            # Allows different column name w.r.t actual name, e.g. x -> long
            f.write(f"# {name}_dimensions: {dim}\n")
            f.write(f"# {name}_column: {field.std_name}\n")
            f.write(f"# {name}_long_name: {field.long_name}\n")
            f.write(f"# {name}_units: {field.unit}\n")
            f.write(f"# {name}_min: {self.data[:, i].min(): {self.fstr}}\n")
            f.write(f"# {name}_max: {self.data[:, i].max(): {self.fstr}}\n")
            if field.grid:
                unique = np.unique(self.data[:, i])
                dx = abs(unique[1] - unique[0])
                f.write(f"# d{name}: {dx}\n")

            # Need to define depth positive
            if name == "z":
                f.write(f"# {name}_positive: {vert_positive}\n")

        if append is not None:
            f.write(append)

    def write_data(self, data=None, f=None):
        """
        Write the data from the .xyz file into the GeoCSV file after the header
        with the originally specified delimiter

        :type data: np.ndarray
        :param data: optional data array to write, if None then automatically 
            read from file or internal attributes
        :type f: _io.TextIOWrapper
        :param f: optional open text file to write to. If None then auto write
            to internal open file
        """
        print("\twriting data to file")
        if f is None:
            f = self.f_out

        # Make sure data is loaded up from the file
        if data is None:
            if self.data is None:
                self.read_xyz()
            data = self.data

        # Write the first uncommented header line which are column header vals
        for i, field in enumerate(self.fields):
            if i < len(self.fields) - 1:
                f.write(f"{field.std_name}{self.delimiter}")
            else:
                f.write(f"{field.std_name}\n")

        # Use Numpy to dump data with correct delimiter
        np.savetxt(f, data, fmt=self.fmt, delimiter=self.delimiter) 

    def close(self):
        """
        Finalize the conversion process
        """
        self.f_out.close()
        print("finished")


def convert_main(input_files, output_files, path_out="./", prepend=None, 
                 append=None, convert_dict=None):
    """
    Main convert script, which takes input, output, and header append files
    but keeps the main processing the same
    """
    # Start em up boys
    conv = Converter(delimiter_in=",", delimiter="|", fmt="%.3f", 
                     path_out=path_out, header_lines=0)

    # Set the values defined by the xyz file. 
    # NOTE: Order matters here, must match the order in which the data appears 
    # in each data line of the input .xyz file! e.g., here it is:
    # x,y,z,lat,lon,vp,vs,rho,qp,qs
    # NOTE: Units are set to how they will appear in the output file, not how
    # they're set in the input file. The function convert_data() may affect
    # what the final units are
    conv.append(std_name="y", long_name="y_axis_utm", unit="km") 
    conv.append(std_name="x", long_name="x_axis_utm", unit="km")
    conv.append(std_name="z", long_name="depth", unit="km")

    # Having lat/lon headers means we will need to run convert_coordinates()
    # BEFORE writing header information. These values are not in the original
    # .xyz file
    conv.append(std_name="lat", long_name="latitude", unit="degrees_north")
    conv.append(std_name="lon", long_name="longitude", unit="degrees_east")

    # Standard seismic tomography model data
    conv.append(std_name="vp", long_name="p_velocity", unit="km/s")
    conv.append(std_name="vs", long_name="s_velocity", unit="km/s")
    conv.append(std_name="rho", long_name="density", unit="kg/m^3")
    conv.append(std_name="qp", long_name="p_attenuation", unit="count")
    conv.append(std_name="qs", long_name="s_attenuation", unit="count")

    # Local paths to the tomography .xyz files
    for fid, fid_out in zip(input_files, output_files):
        conv.set_fids(fid, fid_out)
        conv.read_xyz(swap_idx=[[0, 1]])

        # Convert the UTM60S coordinates to Lat/Lon and insert into the fields
        conv.convert_coordinates(epsg_in=32760, epsg_out=4326,
                                 x_in=conv.data[:,1], y_in=conv.data[:,0],
                                 choice="insert", insert=3)

        # Convert the actual data (e.g., units m -> km). Note 'unit' field is
        # not touched so user must ensure that output units are set correctly
        # during the append() stage
        conv.convert_data(convert=convert_dict)

        # Write all to file
        conv.write_header(prepend=prepend, append=append)
        conv.write_data()
        conv.close()


def convert_nzatom_north_xyz_2_geocsv():
    """
    MAIN function for convering NZAtom_north to GeoCSV file that is formatted to
    match the IRIS emc-tools converters. Can be used as a template for future
    conversions but some of the header information and maybe data information
    will require adjustments.

    * Turning 'grid' off on coordinate data because IRIS doesn't want that
    """
    # !!! Set filenames/filepaths below
    path_in = "./"
    input_files = ["tomography_model_mantle.xyz",
                   "tomography_model_crust.xyz",
                   "tomography_model_shallow.xyz"]

    # output files are auto-controlled by the header information
    path_out = "./"

    # Define header information required by EMC
    header = {
    "global_title": "New Zealand Adjoint Tomography Model - "
                    "North Island (NZ_ATOM_NORTH)",
    "global_id": "nz_atom_north_chow_etal_2021_vp+vs",
    "global_data_revision": "r0.1",
    "global_Conventions": "CF-1.0",
    "global_Metadata_Conventions": "Unidata Dataset Discovery v1.0",
    "global_summary": 
        "NZ_ATOM_NORTH is a 3D velocity model of the North Island of New "
        "Zealand derived using earthquake-based adjoint tomography. The "
        "starting model is defined as the ray-based NZ-Wide2.2 Velocity "
        "Model from Eberhart-Phillips et al. (2021). To derive this "
        "velocity model, we iteratively improved fits between earthquake "
        "obserations from New Zealand-based broadband seismometers and "
        "synthetically generated waveforms from spectral element "
        "simulations. The waveform bandpass of interest is 4-30s. This "
        "velocity model defines the following parameters: Vp (km/s), "
        "Vs (km/s), density (kg/m^3), Qp, and Qs. Only Vp and Vs "
        "are updated during the inversion. The reamining quantities are "
        "defined by the starting/reference velocity model",
    "global_keywords": "adjoint, seismic, earthquake, tomography",
    "global_attribution": "DOI:10.1029/2021JB022865",
    "global_author_name": "Bryant Chow",
    "global_author_contact": "bhchow@alaska.edu",
    "global_author_institution": "University of Alaska Fairbanks",
    "global_repository_name": "EMC",
    "global_repository_institution": "IRIS DMC",
    "global_repository_pid": "doi:10.17611/dp/emc.2021.nzatomnnorthvpvs.1",
    "global_reference_ellipsoid": "WGS 84",
    "global_geodetic_datum": "UTM 60S / EPSG 32760",
    # "global_unit_of_measure": "km",
    # "global_center_coordinates": f"{conv.data[:, 0].mean(): {conv.fstr}}, " 
    #                              f"{conv.data[:, 1].mean(): {conv.fstr}}",
            }

    # Generate full path to filenames for reading/writing
    input_files = [os.path.join(path_in, _) for _ in input_files]
    output_files = [
       f"{header['global_id']}-mantle.{header['global_data_revision']}-n4.csv",
       f"{header['global_id']}-crust.{header['global_data_revision']}-n4.csv",
       f"{header['global_id']}-shallow.{header['global_data_revision']}-n4.csv",
       ]
    output_files = [os.path.join(path_out, _) for _ in output_files]

    convert_dict = {"x": 1E-3, "y": 1E-3, "z": -1E-3, "vp": 1E-3, 
                    "vs": 1E-3}

    convert_main(input_files=input_files, output_files=output_files, 
                 prepend=header, append=None, convert_dict=convert_dict)


def convert_nzwide_north_xyz_2_geocsv():
    """
    MAIN function for converting initial ref model to GeoCSV file, formatted to
    match the IRIS emc-tools converters
    """
    path_in = "./"
    input_files = ["tomography_model_mantle_m00.xyz", 
                   "tomography_model_crust_m00.xyz",
                   "tomography_model_shallow_m00.xyz"]
    input_files = [os.path.join(path_in, _) for _ in input_files]

    path_out = "./"
    output_files = ["ref-model-nzwide2p2-mantle.r0.0-n4.csv",
                    "ref-model-nzwide2p2-crust.r0.0-n4.csv",
                    "ref-model-nzwide2p2-shallow.r0.0-n4.csv"]
    output_files = [os.path.join(path_out, _) for _ in output_files]

    header = {
    "global_title": "New Zealand Wide Velocity Model v2.2 - "
                    "North Island",
    "global_id": "nz_wide2p2_eberhart_phillips_etal_2021",
    "global_data_revision": "r0.0",
    "global_Conventions": "CF-1.0",
    "global_Metadata_Convetions": "Unidata Dataset Discovery v1.0",
    "global_summary": 
        "NZ-Wide2.2 Velocity Model created by Eberhart-Phillips et al. "
        "(2021). This provided reference model has been modified from its "
        "original format for use in our adjoint tomography inversion. "
        "Modifications include: rotation to the UTM-60S coordinate system, "
        "interpolation to a regular grid, and subsequent regularization. "
        "This reference model is parameterized in terms of: "
        "Vp (km/s), Vs (km/s), density (kg/m^3), Qp, and Qs.",
    "global_reference": "Eberhart-Phillips et al. (2021)",
    "global_attribution": "DOI:10.5281/zenodo.3779523",
    "global_repository_name": "EMC",
    "global_repository_institution": "IRIS DMC",
    "global_repository_pid": "doi:10.17611/dp/emc.2021.nzatomnnorthvpvs.1",
    "global_reference_ellipsoid": "WGS 84",
    "global_geodetic_datum": "UTM 60S / EPSG 32760",
            }

    convert_dict = {"x": 1E-3, "y": 1E-3, "z": -1E-3, "vp": 1E-3, 
                    "vs": 1E-3}

    convert_main(input_files=input_files, output_files=output_files, 
                 prepend=header, append=None, convert_dict=convert_dict)


if __name__ == "__main__":
    try:
        if sys.argv[1] == "atom":
            print("CONVERTING NZATOM XYZ VELOCITY MODEL TO GEOCSV")
            convert_nzatom_north_xyz_2_geocsv()
        elif sys.argv[1] == "ref":
            print("CONVERTING NZATOM XYZ VELOCITY MODEL TO GEOCSV")
            convert_nzwide_north_xyz_2_geocsv()
    except IndexError:
        print("argument must be 'atom' or 'ref'")
