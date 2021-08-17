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
    dictionary that ensures all these attributes are met.
    """
    def __init__(self, std_name, long_name, unit):
        self.std_name = std_name
        self.long_name = long_name
        self.unit = unit

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class Converter():
    """
    Class to control the conversion from XYZ to GeoCSV
    """
    def __init__(self, fields=None, delimiter=",", header_lines=4, 
                 isotropic_data_nvals=6, anisotropic_data_nvals=8,
                 fmt="%.1f"):
        """
        Set some constant parameters that the converter uses for internal checks

        :type fields: list of Fields
        :param fields: the different data types that are included in the XYZ
        """
        self.fields = fields or []
        self.delimiter = delimiter
        self.header_lines = header_lines
        self.isotropic_data_nvals = isotropic_data_nvals
        self.anisotropic_data_nvals = anisotropic_data_nvals
        self.fmt = fmt

        # Initiate empty variables to be filled by other functions
        self.fid = None
        self.fid_out = None
        self.f_out = None
        self.data = None
        self.is_isotropic = None

    def append(self, **kwargs):
        """
        Append new Fields to the internal field list
        Kwargs should match the required attributes of the Field class
        """
        self.fields.append(Field(**kwargs))

    def set_fids(self, fid, fid_out=None):
        """
        Point the converter to the input .xyz file, it will then auto set an
        output fid for the GeoCSV or the user can set it themselves.
        """
        if os.path.exists(fid):
            self.fid = fid
            # Auto set the output fid to the be the same thing but with .csv
            if fid_out is None:
                fid_ = os.path.basename(fid)
                assert(".xyz" in fid_), f"input file requires .xyz  extension"
                fid_new = fid_.replace(".xyz", ".csv")
                fid_out = os.path.join(os.getcwd(), fid_new)

            self.fid_out = fid_out
            self.f_out = open(self.fid_out, "w")
        else:
            sys.exit(f"{fid} does not exist and must")

    def read_xyz(self, fid=None):
        """
        Read the .xyz file and store the values as a dict object
        """
        if fid is None:
            fid = self.fid

        self.data = np.loadtxt(fid, skiprows=self.header_lines)

        if len(self.data[0]) == self.isotropic_data_nvals:
            self.is_isotropic = True
        elif len(self.data[0]) == self.anisotropic_data_nvals:
            self.is_isotropic = False
        else:
            sys.exit("Data nvals does not match the expected quantity")

    def convert_coordinates(self, epsg_in, epsg_out, x_in, y_in, 
                            choice="append", insert=None, replace=None):
        """
        Convert coordinates from one coordinate system to something different.
        Useful for NZAtom as people probably want Lat/Lon but the model
        is it UMT60S Cartesian

        :type choice: str
        :param choice: 'append': place new coords at end of fields, 
            'insert': insert new coords at location insert
            'replace': replace given tuple of locations 
        :type insert: int
        :param insert: place in the fields list to insert new coordinates
        :type replace: tuple of int
        :param replace: locations of the field data to be replaced
        """
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

    def write_header(self, f=None, title=None):
        """
        Write the header of the output file using field information. The IRIS
        emc tools requires a specific header structure that it uses to parse 
        data itself so this function tries to mimic that
        """
        # def parse_fields(self, f, name):
        #     f.write(f"# field_{name}: ")
        #     for i, field in enumerate(self.fields):
        #         if i < len(self.fields) - 1:
        #             f.write(f"{field[name]}{self.delimiter}")
        #         else:
        #             f.write(f"{field[name]}\n")

        # write_field(self, f, "unit")
        # write_field(self, f, "dtype")
        # write_field(self, f, "long_name")

        if f is None:
            f = self.f_out

        f.write("# dataset: GeoCSV 2.0\n")
        if title is not None:
            f.write(f"# title: {title}\n")
        f.write(f"# created: {datetime.datetime.utcnow()} UTC\n")
        f.write(f"# delimiter: {self.delimiter}\n")

        assert(len(self.fields) == len(self.data[0])), \
                "Data and field length mismatch, maybe convert coords first"

        for i, field in enumerate(self.fields):
            name = field.std_name
            f.write(f"# {name}_column: {name}\n")
            f.write(f"# {name}_long_name: {field.long_name}\n")
            f.write(f"# {name}_units: {field.unit}\n")
            f.write(f"# {name}_min: {self.data[:, i].min()}\n")
            f.write(f"# {name}_max: {self.data[:, i].max()}\n")


    def append_header(self, line, f=None):
        """
        Append any extra lines to the header that does not come standard with
        the field information, e.g. reference ellipsoid.
        Automatically starts header lines with '#' and ends with newline
        """
        if f is None:
            f = self.f_out
        f.write(f"# {line}\n")

    def write_data(self, data=None, f=None):
        """
        Write the data read from the xyz file
        """
        if f is None:
            f = self.f_out

        # Make sure data is loaded up from the file
        if data is None:
            if self.data is None:
                self.read_xyz()
            data = self.data

        # Write the first header line 
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


if __name__ == "__main__":
    # Load and run the converter
    conv = Converter(delimiter=",", fmt="%.3f")

    # Set the values defined by the xyz file. Order matters here, must match
    # the order in which the data appears in each data line!
    conv.append(std_name="x", long_name="x_axis_utm", unit="m")
    conv.append(std_name="y", long_name="y_axis_utm", unit="m")
    conv.append(std_name="depth", long_name="z_axis_utm", unit="m")

    # Having lat/lon headers means you need to run convert_coordinates()
    conv.append(std_name="lat", long_name="latitude", unit="degrees_north")
    conv.append(std_name="lon", long_name="longitude", unit="degrees_east")

    conv.append(std_name="vp", long_name="p_velocity", unit="m/s")
    conv.append(std_name="vs", long_name="s_velocity", unit="m/s")
    conv.append(std_name="rho", long_name="density", unit="kg/m^3")

    conv.append(std_name="qp", long_name="p_attenuation", unit="count")
    conv.append(std_name="qs", long_name="s_attenuation", unit="count")
    
    path_ = "/Users/Chow/Documents/academic/vuw/data/tomo_files/birch_m11"
    for tag in ["mantle"]:#, "crust", "shallow"]:
        fid = os.path.join(path_, f"tomography_model_{tag}.xyz")

        conv.set_fids(fid)
        conv.read_xyz()

        # Convert the UTM60S coordinates to Lat/Lon and insert into the fields
        conv.convert_coordinates(epsg_in=32760, epsg_out=4326, 
                                 x_in=conv.data[:,0], y_in=conv.data[:,1],
                                 choice="insert", insert=3)

        conv.write_header(title="New Zealand Adjoint Tomography Model - "
                                "North Island (NZ_ATOM_NORTH)")

        # Add more header information here, e.g. reference ellipsoid, 
        conv.append_header(f"id: nz_atom_north_chow_etal_2021_vp+vs")
        conv.append_header(f"author_name: Bryant Chow")
        conv.append_header(f"author_contact: bryant.chow@vuw.ac.nz")
        conv.append_header(f"attribution: ???")
        conv.append_header(f"reference_ellipsoid: WGS 84")
        conv.append_header(f"geodetic_datum: UTM 60S / EPSG 32760")
        conv.append_header(f"unit_of_measure: m")
        conv.append_header(f"center_coordinates: 495732.01, 5572241.58")
        conv.append_header(f"vertical_positive: up")


        conv.write_data()
        conv.close()


