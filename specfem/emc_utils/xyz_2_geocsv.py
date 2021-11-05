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
        depth must be depth because part of their code fails if it's z

Relevant links:
    https://github.com/iris-edu/emc-tools
    http://geows.ds.iris.edu/documents/GeoCSV.pdf
"""
import os
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
                 path_out=os.getcwd()):
        """
        Set some constant parameters that the converter uses for internal checks

        :type fields: list of Fields
        :param fields: the different data types that are included in the XYZ
        :type delimiter: str
        :param delimiter: how to separate variables in the file
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
        """
        print("Initiating Converter for .xyz to GeoCSV")
        self.fields = fields or []
        self.delimiter = delimiter
        self.header_lines = header_lines
        self.fmt = fmt
        self.path_out = path_out

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
            sys.exit(f"{fid} does not exist and must")

        print(f"\tfile to read: {fid}")
        print(f"\toutput file will be: {fid_out}")

    def read_xyz(self, fid=None):
        """
        Read the .xyz file and store the values as a dict object


        :type fid: str
        :param fid: read the input xyz file with numpy loadtxt
        """
        print(f"\treading input .xyz file")
        if fid is None:
            fid = self.fid

        self.data = np.loadtxt(fid, skiprows=self.header_lines)

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
        :type prepend: str
        :param prepend: optional lines to write after the first header
            declaration. MUST be formatted with a leading '#' and a trailing
            '\n' to match the header format
        :type append: str
        :param append: optional lines to write after the final standard header
            declaration. MUST be formatted with a leading '#' and a trailing
            '\n' to match the header format
        """
        print("\twriting header information")
        if f is None:
            f = self.f_out

        f.write("# dataset: GeoCSV 2.0\n")
        if prepend is not None:
            f.write(prepend)

        # Write some standard information
        f.write(f"# created: {datetime.datetime.utcnow()} UTC\n")
        f.write(f"# delimiter: {self.delimiter}\n")

        assert(len(self.fields) == len(self.data[0])), \
                "Data and field length mismatch, maybe convert coords first"

        # Write header information REQUIRED by GeoCSV and IRIS emc-tools
        for i, field in enumerate(self.fields):
            name = field.std_name
            f.write(f"# {name}_column: {name}\n")
            f.write(f"# {name}_long_name: {field.long_name}\n")
            f.write(f"# {name}_units: {field.unit}\n")
            f.write(f"# {name}_min: {self.data[:, i].min()}\n")
            f.write(f"# {name}_max: {self.data[:, i].max()}\n")
            if field.grid:
                unique = np.unique(self.data[:, i])
                dx = abs(unique[1] - unique[0])
                f.write(f"# d{name}: {dx}\n")

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

def convert_nzatom_north_xyz_2_geocsv():
    """
    MAIN function for convering NZAtom_north to GeoCSV file that is formatted to
    match the IRIS emc-tools converters. Can be used as a template for future
    conversions but some of the header information and maybe data information
    will require adjustments.
    """
    # !!! Set the files here
    path_in = "./"
    input_files = ["tomography_model_mantle.xyz"]#, 
                   #"tomography_model_crust.xyz",
                   #"tomography_model_shallow.xyz"]
    input_files = [os.path.join(path_in, _) for _ in input_files]
    path_out = "./"


    # Start em up boys
    conv = Converter(delimiter=",", fmt="%.3f", path_out=path_out)

    # Set the values defined by the xyz file. 
    # NOTE: Order matters here, must match the order in which the data appears 
    # in each data line of the input .xyz file! e.g., here it is:
    # x,y,z,lat,lon,vp,vs,rho,qp,qs
    conv.append(std_name="x", long_name="x_axis_utm", unit="m", grid=True)
    conv.append(std_name="y", long_name="y_axis_utm", unit="m", grid=True)
    conv.append(std_name="depth", long_name="z_axis_utm", unit="m", grid=True)

    # Having lat/lon headers means we will need to run convert_coordinates()
    # BEFORE writing header information. These values are not in the original
    # .xyz file
    conv.append(std_name="lat", long_name="latitude", unit="degrees_north")
    conv.append(std_name="lon", long_name="longitude", unit="degrees_east")

    # Standard seismic tomography model data
    conv.append(std_name="vp", long_name="p_velocity", unit="m/s")
    conv.append(std_name="vs", long_name="s_velocity", unit="m/s")
    conv.append(std_name="rho", long_name="density", unit="kg/m^3")
    conv.append(std_name="qp", long_name="p_attenuation", unit="count")
    conv.append(std_name="qs", long_name="s_attenuation", unit="count")

    # Local paths to the tomography .xyz files
    for fid in input_files:
        conv.set_fids(fid)
        conv.read_xyz()

        # Convert the UTM60S coordinates to Lat/Lon and insert into the fields
        conv.convert_coordinates(epsg_in=32760, epsg_out=4326,
                                 x_in=conv.data[:,0], y_in=conv.data[:,1],
                                 choice="insert", insert=3)

        conv.write_header(
                prepend="# title: New Zealand Adjoint Tomography Model - "
                        "North Island (NZ_ATOM_NORTH)\n"
                        "# id: nz_atom_north_chow_etal_2021_vp+vs\n"
                        "# author_name: Bryant Chow\n"
                        "# author_contact: bryant.chow@vuw.ac.nz\n"
                        "# attribution: DOI:10.1002/essoar.10507657.1.\n"
                        "# reference_ellipsoid: WGS 84\n"
                        "# geodetic_datum: UTM 60S / EPSG 32760\n"
                        "# unit_of_measure: m\n"
                        "# center_coordinates: 495732.01, 5572241.58\n"
                        "# vertical_positive: up\n",
                append=None
                        )

        conv.write_data()
        conv.close()


def convert_nzwide_north_xyz_2_geocsv():
    """
    MAIN function for convering NZAtom_north to GeoCSV file that is formatted to
    match the IRIS emc-tools converters. Can be used as a template for future
    conversions but some of the header information and maybe data information
    will require adjustments.
    """
    # !!! Set the files here
    path_in = "/Users/Chow/Documents/academic/vuw/data/tomo_files/nznorth19"
    input_files = ["tomography_model_mantle.xyz", 
                   "tomography_model_crust.xyz",
                   "tomography_model_shallow.xyz"]

    input_files = [os.path.join(path_in, _) for _ in input_files]
    path_out = ("/Users/Chow/Documents/academic/vuw/data/tomo_files/nzatom/"
                "initial")

    # Start em up boys
    conv = Converter(delimiter=",", fmt="%.3f", path_out=path_out)

    # Set the values defined by the xyz file. 
    # NOTE: Order matters here, must match the order in which the data appears 
    # in each data line of the input .xyz file! e.g., here it is:
    # x,y,z,lat,lon,vp,vs,rho,qp,qs
    conv.append(std_name="x", long_name="x_axis_utm", unit="m", grid=True)
    conv.append(std_name="y", long_name="y_axis_utm", unit="m", grid=True)
    conv.append(std_name="depth", long_name="z_axis_utm", unit="m", grid=True)

    # Having lat/lon headers means we will need to run convert_coordinates()
    # BEFORE writing header information. These values are not in the original
    # .xyz file
    conv.append(std_name="lat", long_name="latitude", unit="degrees_north")
    conv.append(std_name="lon", long_name="longitude", unit="degrees_east")

    # Standard seismic tomography model data
    conv.append(std_name="vp", long_name="p_velocity", unit="m/s")
    conv.append(std_name="vs", long_name="s_velocity", unit="m/s")
    conv.append(std_name="rho", long_name="density", unit="kg/m^3")
    conv.append(std_name="qp", long_name="p_attenuation", unit="count")
    conv.append(std_name="qs", long_name="s_attenuation", unit="count")

    # Local paths to the tomography .xyz files
    for fid in input_files:
        conv.set_fids(fid)
        conv.read_xyz()

        # Convert the UTM60S coordinates to Lat/Lon and insert into the fields
        conv.convert_coordinates(epsg_in=32760, epsg_out=4326,
                                 x_in=conv.data[:,0], y_in=conv.data[:,1],
                                 choice="insert", insert=3)

        conv.write_header(
                prepend="# title: New Zealand Wide Velocity Model v2.2- "
                        "North Island\n"
                        "# id: nz_wide2p2_eberhart_phillips_etal_2020\n"
                        "# author_name: Donna Eberhart-Phillips et al.\n"
                        "# attribution: DOI:10.5281/zenodo.3779523\n"
                        "# reference_ellipsoid: WGS 84\n"
                        "# geodetic_datum: UTM 60S / EPSG 32760\n"
                        "# unit_of_measure: m\n"
                        "# center_coordinates: 495732.01, 5572241.58\n"
                        "# vertical_positive: up\n",
                append=None
                        )

        conv.write_data()
        conv.close()


if __name__ == "__main__":
    convert_nzwide_north_xyz_2_geocsv()
    # convert_nzatom_north_xyz_2_geocsv()
