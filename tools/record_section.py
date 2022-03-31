#!/usr/bin/evn python3
"""
Record Section plotting tool for seismic waveforms (observed and synthetic)

This is a refactor of Pysep's Python utility `plotw_rs`, a record section
plotting script. The intent of this script is to plot multiple time series'
based on source-receiver characteristics (i.e., src-rcv distance, backazimuth).

.. note:: Code History
    Written by Carl Tape (11/21/2011) and Yun Wang (11/2011) in Matlab
    Translated to Python by Nealy Sims (1/2021)
    Upgraded by Aakash Gupta (9/2021)
    Refactored by Bryant Chow (3/2022)
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
from matplotlib.ticker import MultipleLocator
from obspy.geodetics import (kilometers2degrees, degrees2kilometers,
                             gps2dist_azimuth)


def myround(x, base=5, choice="near"):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base
    Copied from Pyatoa

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == "near":
        roundout = int(base * round(float(x)/base))
    elif choice == "down":
        roundout = int(base * np.floor(float(x)/base))
    elif choice == "up":
        roundout = int(base * np.ceil(float(x)/base))
    return roundout


class Dict(dict):
    """Easy dictionary overload for nicer get/set attribute characteristics"""
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class RecordSection:
    """
    Record section plotting tool which takes ObsPy streams and preprocesses and
    sorts source-receiver pairs based on User input, and produces record section
    waveform figures
    """
    def __init__(self, st, st_syn=None, sort_by="default", scale_amp=None,
                 time_shift_s=None, move_out=None, time_marker=None,
                 min_period_s=None, max_period_s=None, preprocess="st",
                 max_traces_per_rs=50, integrate=0, xlim_s=None,
                 components="ZRTNE12",
                 y_axis_spacing=1, sort_by_azimuth_start_deg=0.,
                 distance_units="km", geometric_spreading_factor=0.5,
                 geometric_spreading_k_val=None, figsize=(9, 11), show=True,
                 save="./record_section.png", overwrite=False,
                 trace_kwargs=None, **kwargs):
        """
        Set the default record section plotting parameters and enforce types
        Run some internal parameter derivation functions by manipulating input
        data and parameters.

        :type st: obspy.core.stream.Stream
        :param st: Stream objects containing observed time series to be plotted
            on the record section. Can contain any number of traces
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream objects containing synthetic time series to be
            plotted on the record section. Must contain the same number of
                traces as `st`
        :type sort_by: str
        :param sort_by: How to sort the Y-axis of the record section, available:
            - None: Not set, don't sort, just iterate directly through Stream
            - 'alphabetical': sort alphabetically
            - 'azimuth': sort by source-receiver azimuth (deg) with constant
                vertical spacing
            - 'backazimuth': sort by source-receiver backazimuth (deg) with
                constant vertical spacing. Requires `sort_by_azimuth_start_deg`
            - 'distance': sort by source-receiver distance (km) with constant
                vertical spacing. Requires `sort_by_azimuth_start_deg` AND
                `distance_units`
            - 'abs_distance': absolute vertical spacing of waveforms defined by
                source-receiver distance (km). Requires `distance_units`
            - 'abs_azimuth': absolute vertical spacing of waveforms defined
                by source-receiver azimuth (deg). Requires
                `sort_by_azimuth_start_deg`
            - 'abs_backazimuth': absolute vertical spacing of waveforms by
                source-receiver backazimuth (deg).
            - '*_r': Add a '_r' to any of the values about to REVERSE the sort,
                e.g., alphabetical_r sort will go Z->A
        :type scale_amp: list
        :param scale_amp: scale amplitude of waveforms by available:
            - None: Not set, no amplitude scaling, waveforms shown raw
            - 'normalize': scale each trace by the maximum amplitude,
                i.e., > a /= max(abs(a))  # where 'a' is time series amplitudes
            - 'geometric_spreading': scale amplitudes globally by predicting the
                expected geometric spreading amplitude reduction and correcting
                for this factor. Requires `geometric_spreading_factor`, optional
                `geometric_spreading_k_val`
        :type time_shift_s: float OR list of float
        :param time_shift_s: apply static time shift to waveforms, two options:
            1. float (e.g., -10.2), will shift ALL waveforms by
                that number (i.e., -10.2 second time shift applied)
            2. list (e.g., [5., -2., ... 11.2]), will apply individual time
                shifts to EACH trace in the stream. The length of this list MUST
                match the number of traces in your input stream.
        :type move_out: float
        :param move_out: Optional. A velocity value that will be used to
            calculate move out, which will time shift seismograms based on
            their source receiver distance. This parameter will be ADDED
            to time_shift_s (both float and list), if it is provided.
            Should be in units of `distance_units`/s
        :type time_marker: list of float
        :param time_marker: absolute time markers. If this parameter is set,
            `time_shift_s` MUST be set as Option 1 (single float)
        :type min_period_s: float
        :param min_period_s: minimum filter period in seconds
        :type max_period_s: float
        :param max_period_s: maximum filter period in seconds
        :type preprocess: str
        :param preprocess: choose which data to preprocess, options are:
            - 'st': process waveforms in st
            - 'st_syn': process waveforms in st_syn. st still must be given
            - 'both': process waveforms in both st and st_syn
            - None: do not run preprocessing
        :type max_traces_per_rs: int
        :param max_traces_per_rs: maximum number of traces to show on a single
            record section plot
        :type xlim_s: list of float
        :param xlim_s: [start, stop] in units of time, seconds, to set the
            xlimits of the figure
        :type components: str
        :param components: a sequence of strings representing acceptable
            components from the data. Also determines the order these are shown
            EVEN when sorted by other variables. For example, components=='ZR'
            would only display Z and R components, and Z components would be
            should BEFORE R components for the SAME station.
        :type integrate: int
        :param integrate: apply integration `integrate` times on all traces.
            acceptable values [-inf, inf], where positive values are integration
            and negative values are differentiation
            e.g., if integrate == 2,  will integrate each trace twice.
            or    if integrate == -1, will differentiate once
            or    if integrate == 0,  do nothing
        :type y_axis_spacing: float
        :param y_axis_spacing: spacing between adjacent seismograms applied to
            Y-axis on relative (not absolute) scales.
            Set to 1 for default behavior
        :type sort_by_azimuth_start_deg: float
        :param sort_by_azimuth_start_deg: azimuthal angle for the top record,
            applied to any Y-axis sorting which has 'azimuth' in its name.
            Set to 0 for default behavior
        :type distance_units: str
        :param distance_units: Y-axis units epicentral distance record sections
            'km': kilometers on the sphere
            'deg': degrees on the sphere
            'km_utm': kilometers on flat plane, UTM coordinate system
        :type geometric_spreading_factor: float
        :param geometric_spreading_factor: geometrical spreading factor when
            using the `scale_amp` parameter. Defaults to 0.5 for surface waves.
            Use values of 0.5 to 1.0 for regional surface waves
        :type geometric_spreading_k_val: float
        :param geometric_spreading_k_val: K value used to scale the geometric
            spreading factor (TODO figure out what this actually is)
        :type figsize: tuple of float
        :param figsize: size the of the figure, passed into plt.subplots()
         :type show: bool
        :param show: show the figure as a graphical output
        :type save: str
        :param save: path to save output figure, will create the parent
            directory if it doesn't exist. If None, will not save.
        :type overwrite: bool
        :param overwrite: if the path defined by `save` exists, will overwrite
            the existing figure
        :type kwargs: dict
        :param kwargs: keyword arguments, expected as a dictionary, will be
            passed around to plotting functions.
        :rtype: Dict
        :return: a Dictionary object containing all the user set parameters
        :raises AssertionError: if any parameters are set incorrectly
        """
        # User defined parameters, do some type-setting
        self.st = st.copy()
        try:
            self.st_syn = st_syn.copy()
        except AttributeError:
            self.st_syn = None

        # Y-Axis sorting parameters
        self.sort_by = sort_by.lower()
        self.y_axis_spacing = float(y_axis_spacing)
        self.sort_by_azimuth_start_deg = float(sort_by_azimuth_start_deg)
        self.components = str(components)

        # Amplitude scaling parameters
        self.scale_amp = scale_amp
        self.geometric_spreading_factor = float(geometric_spreading_factor)
        self.geometric_spreading_k_val = geometric_spreading_k_val

        # Time shift parameters
        self.move_out = move_out
        self.time_shift_s = time_shift_s
        self.time_marker = time_marker

        # Filtering parameters
        self.min_period_s = min_period_s
        self.max_period_s = max_period_s
        self.preprocess = preprocess
        self.max_traces_per_rs = int(max_traces_per_rs)
        self.integrate = int(integrate)

        # Plotting parameters
        self.xlim_s = xlim_s
        self.distance_units = distance_units.lower()
        self.figsize = figsize
        self.show = bool(show)
        self.save = save
        self.overwrite = bool(overwrite)
        self.trace_kwargs = Dict(trace_kwargs)
        self.kwargs = Dict(kwargs)

        # Run checks to ensure that all the parameters are set properly
        # And get some new, initial parameters that are required for other
        # 'get' functions
        self.check_parameters()
        self.stats = self.get_srcrcv_stats()
        self.distances, self.azimuths, self.backazimuths = \
            self.get_srcrcv_dist_az_baz()

        # Internally used parameters that will be filled out by other functions
        self.f = None
        self.ax = None
        self.idx = []
        self.station_ids = []
        self.max_amplitudes = []
        self.amplitude_scaling = []
        self.y_axis = []
        self.sorted_idx = []

    def check_parameters(self):
        """
        Check that parameters are set properly and in line with how they
        are expected by the program

        .. note::
            Not using assertions here because we want all incorrect parameters
            to be evaluated together and displayed at once, that way the user
            doesn't have to run this function multiple times to figure out how
            to set their parameters correct

        :raises AssertionError: If any parameters are not set as expected by
            plotw_rs functionality
        """
        print("checking parameter acceptability")

        # Used to keep track of which parameters failed and in what way
        err = Dict()

        # Check to make sure there is data
        if not bool(self.st):
            err.st = f"stream has {len(self.st)} traces, no data"

        # Check that stream has SAC headers if we want to sort by dist or (b)az.
        # Pysep should have created these
        if self.sort_by != "default" and \
                any(_ in self.sort_by for _ in ["azimuth", "distance"]):
            _idx = []
            for i, tr in enumerate(self.st):
                if not hasattr(tr.stats, "sac"):
                    _idx.append(i)
            if _idx:
                err.st = (f"{len(_idx)} traces have no SAC header, plotw_rs "
                          f"expects SAC headers for sorting. Trace indexes "
                          f"are: {_idx}")

        if self.st_syn is not None:
            if len(self.st) != len(self.st_syn):
                err.st_syn = f"length must match `st` (which is {len(self.st)})"

        # Check the `sort_by` sorting parameter options
        acceptable_sort_by = ["default", "azimuth", "backazimuth",
                              "distance", "alphabetical", "abs_azimuth",
                              "abs_distance"]
        # Allow reverse sorts
        acceptable_sort_by += [f"{_}_r" for _ in acceptable_sort_by]
        if self.sort_by not in acceptable_sort_by:
            err.sort_by = f"must be in {acceptable_sort_by}"

        if "azimuth" in self.sort_by:
            if not (0 <= self.sort_by_azimuth_start_deg <= 360):
                err.sort_by_azimuth_start_deg = f"0 < azi < 360"

        acceptable_distance_units = ["km", "km_utm", "deg"]
        if ("distance" in self.sort_by) and \
                (self.distance_units not in acceptable_distance_units):
            err.sort_by_azimuth_start_deg = \
                f"must be in {acceptable_distance_units}"

        if self.scale_amp is not None:
            acceptable_scale_amp = ["normalize", "geometric_spreading"]
            if self.scale_amp not in acceptable_scale_amp:
                err.scale_amp = f"must be in {acceptable_scale_amp}"

        if self.time_shift_s is not None:
            acceptable_time_shift_s = [1, len(self.st)]
            if len(self.time_shift_s) not in acceptable_time_shift_s:
                err.time_shift_s = f"must be in {acceptable_time_shift_s}"

        if self.time_marker is not None:
            if (self.time_shift_s is not None) and (self.time_shift_s != 1):
                err.time_marker = ("time_shift_s cannot be variable when "
                                        "using this option")

        if self.min_period_s is not None and self.max_period_s is not None:
            if self.min_period_s >= self.max_period_s:
                err.min_period_s = "must be less than `max_period_s`"

        if self.preprocess is not None:
            acceptable_preprocess = ["both", "st"]
            if self.preprocess not in acceptable_preprocess:
                err.preprocess = f"must be in {acceptable_preprocess}"

        if self.xlim_s is not None:
            if len(self.xlim_s) != 2:
                err.xlim_s = f"must be of length 2, [start, stop]"
            elif self.xlim_s[0] > self.xlim_s[1]:
                err.xlim_s = f"start time must be less than stop time"

        if len(self.figsize) != 2:
            err.figsize = "must be tuple defining (horizontal, vertical) extent"

        if os.path.exists(self.save) and not self.overwrite:
            err.save = f"path {self.save} already exists, will not overwrite"

        _dirname = os.path.dirname(self.save)
        if not os.path.exists(_dirname):
            print(f"creating output directory {_dirname}")
            os.makedirs(_dirname)

        if err:
            out = "ERROR - Parameter errors, please make following changes:\n"
            out += "\n".join([f"\t{key}: {val}" for key, val in err.items()])
            print(out)
            sys.exit(-1)

    def get_skip_idx(self):
        """
        Get a list of any traces that don't adhere to user-defined boundaries
        such as dist, az, baz, id, or component matches. Don't actually remove
        the traces from the stream but rather just collect indices we can use
        to skip when plotting

        TODO add distance, azi and backazi criteria
        :rtype: np.array
        :return: returns an indexing list which can be used to skip over
            traces that don't adhere to certain criteria
        """
        skip_idx = []
        for idx in self.idx:
            tr = self.st[idx]
            # Component-wise removal
            if tr.stats.component not in self.components:
                skip_idx.append(idx)
            # !!! Add more here
        print(f"criteria check will remove "
              f"{len(skip_idx)}/{len(self.st)} traces")

        return np.array(skip_idx)

    def get_parameters(self):
        """
        Calculate parameters in a specific order and based on the user-defined
        information. Definitions of parameters are below:

        Calculated Parameters
        ::
            np.array idx:
                a linear indexing of all the traces in the stream
            np.array station_ids:
                an ordered list of station ids, used to get station names
                that match the index defined in `idx`
            np.array max_amplitudes:
                abs max amplitudes of each trace, used for normalization
            np.array amplitude_scaling:
                An array to scale amplitudes based on user choices
            np.array y_axis:
                Y-Axis values based on sorting algorithm, used for plotting
            np.array distances:
                source-receiver distances in `distance_units` units
            np.array azimuths:
                source-receiver azimuths in degrees
            np.array backazimuths:
                source-receiver backazimuths in degrees
            np.array sorted_idx:
                sorted indexing on all the traces of the stream based on the
                chosen sorting algorithm
        """
        self.idx = np.arange(0, len(self.st), 1)
        self.station_ids = np.array([tr.get_id() for tr in self.st])
        self.max_amplitudes = np.array([max(abs(tr.data)) for tr in self.st])

        # Get array-like factors which control the sort order, amplitudes, etc.
        sorted_idx = self.get_sorted_idx()
        skip_idx = self.get_skip_idx()
        # Remove skip indexes from sorted index to get the final ordered
        # list of traces to plot
        self.sorted_idx = np.array([_ for _ in sorted_idx if _ not in skip_idx])

        self.y_axis = self.get_y_axis_positions()
        self.amplitude_scaling = self.get_amplitude_scaling()
        self.time_shift_s = self.get_time_shifts()  # OVERWRITES user input

    def get_time_shifts(self):
        """
        Very simple function which allows float inputs for time shifts and
        ensures that time shifts are always per-trace arrays

        :rtype: np.array
        :return: a stream-lengthed array of time shifts that can be applied
            per trace
        """
        # No user input means time shifts will be 0, so nothing happens
        time_shift_arr = np.zeros(len(self.st))
        if self.time_shift_s is not None:
            # User inputs a static time shift
            if isinstance(self.time_shift_s, (int, float)):
                time_shift_arr += self.time_shift_s
            # User input an array which should have already been checked for legnth
            else:
                time_shift_arr = self.time_shift_s
        time_shift_arr = np.array(time_shift_arr)

        # Further change the time shift if we have move out input
        if self.move_out is not None:
            print(f"apply {self.move_out} {self.distance_units}/s move out")
            move_out_arr = self.distances / self.move_out
            time_shift_arr -= move_out_arr

        return time_shift_arr

    def get_srcrcv_dist_az_baz(self):
        """
        Convenience function to wrap _get_srcrcv_info_trace into a loop over the
        whole stream and return lists of distances, azimuths, and backazimuths

        :rtype distances: np.array
        :return distances: source-receiver distances in user-defined units in
            the original order of Stream
        :rtype azimuths: np.array
        :return azimuths: source-receiver azimuths (deg) in the original
            order of Stream
        :rtype backazimuths: np.array
        :return backazimuths: source-receiver azimuths (deg) in the original
            order of Stream
        """
        print("calculating source-receiver distance and (back)azimuths")

        distances, azimuths, backazimuths = [], [], []
        for tr in self.st:
            gcd, az, baz = self._get_srcrcv_dist_az_baz_trace(tr=tr)
            distances.append(gcd)
            azimuths.append(az)
            backazimuths.append(baz)

        distances = np.array(distances)
        azimuths = np.array(azimuths)
        backazimuths = np.array(backazimuths)

        return distances, azimuths, backazimuths

    def _get_srcrcv_dist_az_baz_trace(self, tr=None, idx=0):
        """
        Check the source-receiver characteristics such as src-rcv distance,
        azimuth, backazimuth for a given trace.

        .. note::
            This script ASSUMES that SAC headers have been written to the
            traces. Otherwise we will need more complicated ways to get
            event lat and lon

        :type tr: obspy.core.trace.Trace
        :param tr: trace to get srcrcv information for. If None, will use `idx`
        :type idx: int
        :param idx: if no trace given, give index to check self.st (Stream),
            defaults to index 0
        :rtype gcdist: float
        :return gcdist: great circle distance in units specified by
            `distance_units`, can be 'km' or 'deg'
        :rtype az: float
        :return az: azimuth (degrees) between source and receiver
        :rtype baz: float
        :return baz: azimuth (degrees) between source and receiver
        """
        # Default to the first trace in the Stream
        if tr is None:
            tr = self.st[int(idx)]
        try:
            dist = tr.stats.sac.dist  # units: km
            az = tr.stats.sac.az        # units: deg
            baz = tr.stats.sac.baz      # units: deg
        except AttributeError:
            # If for whatever reason SAC headers dont contain this info already
            # Use ObsPy to get great circle distance, azimuth and backazimuth
            dist, az, baz = gps2dist_azimuth(lat1=tr.stats.sac.evla,
                                             lon1=tr.stats.sac.evlo,
                                             lat2=tr.stats.sac.stla,
                                             lon2=tr.stats.sac.stlo)
            dist *= 1E-3   # units: m -> km

        # Ensure that 0 <= (b)az <= 360
        az = az % 360
        baz = baz % 360

        # Make sure distance is in the correct units, default units 'km'
        if self.distance_units == "deg":
            dist = kilometers2degrees(dist)  # units: km -> deg
        elif self.distance_units == "km_utm":
            # Overwrite `dist`, could probably skip that calc above but
            # leaving for now as I don't think this option will be used heavily.
            dist_deg = np.sqrt(
                ((tr.stats.sac.stlo - tr.stats.sac.evlo) ** 2) /
                ((tr.stats.sac.stla - tr.stats.sac.evla) ** 2)
            )
            dist = kilometers2degrees(dist_deg)  # units: km

        return dist, az, baz

    def get_srcrcv_stats(self):
        """
        Get source receiver information such as min max values, and
        count-related numbers (e.g., num stations) to be used mainly for print
        statements and text information

        Stats Arguments
        ::
            np.array event_names:
                unique event names taken from the SAC header
            int nevents:
                number of unique events in the stream
            np.array unique_sta_ids:
                unique station codes taken from trace stats
            int nstation_ids:
                number of unique station codes
            np.array network_codes:
                unique network codes taken from station ids
            int nnetwork:
                number of unique network codes
            np.array station_codes:
                unique station codes taken from station ids
            int nstation:
                number of unique station codes
            np.array location_codes:
                unique location codes taken from station ids
            int nlocation:
                number of unique location codes
            np.array channel_codes:
                unique channel codes taken from station ids
            int nchannel:
                number of unique channel codes
            bool reverse_sort:
                determine if the user wants to reverse their sort, they do this
                by appending '_r' to the end of the `sort_by` argument
        """
        print("getting source-receiver stats")

        def _unique(list_):
            """return a unique numpy array derived from a list"""
            return np.unique(np.array(list_, dtype=str))

        stats = Dict()
        stats.event_names = _unique([tr.stats.sac.kevnm for tr in self.st])
        stats.nevents = len(stats.event_names)
        stats.unique_sta_ids = _unique([tr.get_id() for tr in self.st])
        stats.longest_id = max([len(_) for _ in stats.unique_sta_ids])
        stats.nstation_ids = len(stats.unique_sta_ids)
        # Get unique network, station, location and channel codes. Also numbers
        for name in ["network", "station", "location", "channel", "component"]:
            stats[f"{name}_codes"] = _unique(
                [getattr(tr.stats, name) for tr in self.st]
            )
            stats[f"n{name}"] = len(stats[f"{name}_codes"])
        # Include the `not` in `reverse_sort` to make the top of the y-axis the
        # starting point, which seems more intuitive for record sections, but
        # is opposite the behavior when you increment from 0
        stats.reverse_sort = not bool("_r" in self.sort_by)

        # Initiate empty lists for _plot_trace to fill with min and max data
        # values which can be used for global plotting parameters like xlims
        stats.xmin, stats.xmax, stats.ymin, stats.ymax = [], [], [], []

        return stats

    def get_amplitude_scaling(self):
        """
        Scale the amplitudes of all the waveforms by producing a Stream
        dependent scale factor based on user choice. It is expected that the
        output array will be DIVIDED by the data arrays

        .. note::
            Needs to be run AFTER preprocessing because filtering etc. will
            affect the final amplitudes of the waveforms

        :rtype: np.array
        :return: an array corresponding to the Stream indexes which provides
            a per-trace scaling coefficient
        """
        print(f"determining amplitude scaling with: {self.scale_amp}")

        # Don't scale by anything
        if self.scale_amp is None:
            amp_scaling = np.ones(len(self.st))
        # Scale by the max amplitude of each trace
        elif self.scale_amp == "normalize":
            amp_scaling = self.max_amplitudes
            # When using absolute distance scale, scale waveforms to minmax dist
            if "abs" in self.sort_by:
                if "distance" in self.sort_by:
                    print("scaling amplitudes for absolute distance")
                    # !!! This scale is from trial and error, it may need to be
                    # !!! made more robust in the future
                    scale = np.mean(self.distances) / 10
                elif "backazimuth" in self.sort_by:
                    print("scaling amplitudes for absolute backazimuth")
                    scale = self.backazimuths.max() - self.backazimuths.min()
                elif "azimuth" in self.sort_by:
                    print("scaling amplitudes for absolute azimuth")
                    scale = self.azimuths.max() - self.azimuths.min()
            # Divide to make waveforms larger
            amp_scaling /= scale
        # Scale by the theoretical geometrical spreading factor
        elif self.scale_amp == "geometric_spreading":
            amp_scaling = self._calculate_geometric_spreading()

        return amp_scaling

    def _calculate_geometric_spreading(self):
        """
        Stations with larger source-receiver distances will have their amplitude
        scaled by a larger value.

        For information on geometric spreading, see Stein and Wysession,
        Section 4.3.4, Eq 20 (for Rayleigh waves, geometrical spreading
        factor = 0.5). For our purposes, this will fold the attenuation factor
        into the same single factor that accounts for geometric spreading.

        .. note::
            This does not take into account the variation in amplitude.

        TODO Plot geometric spreading with best-fitting curve

        :type max_amplitudes: list
        :param max_amplitudes: list of maximum amplitudes for each trace in the
            stream. IF not given, will be calculated on the fly. Optional incase
            this has been calculated already and you want to save time
        :rtype: list
        :return: scale factor per trace in the stream based on theoretical
            geometrical spreading factor. This is meant to be MULTIPLIED by the
            data arrays
        """
        print("calculating geometrical spreading for amplitude normalization")

        # Create a sinusoidal function based on distances in degrees
        sin_del = np.sin(np.array(self.dist) / (180 / np.pi))

        # !!! TODO look at Stein and Wysession and figure out what these vector
        # !!! TODO names are, sort of ambiguous right meow
        if self.geometric_spreading_k_val is not None:
            k_vector = self.max_amplitudes * \
                       (sin_del ** self.geometric_spreading_factor)
            k_val = np.median(k_vector)
        else:
            k_val = self.geometric_spreading_k_val

        w_vector = k_val / (sin_del ** self.geometric_spreading_factor)

        return w_vector

    def get_sorted_idx(self):
        """
        Sort the source-receiver pairs by the chosen parameters. Except for in
        `default` sorting, always maintains component ordering defined by the
        user, i.e., for the same station, components will be ordered by the
        `component` parameter.

        :rtype: np.array
        :return: returns an indexing list which is sorted based on the user
            defined `sort_by` argument.
        """
        print(f"determining sort order with parameter: {self.sort_by}")

        # Retain input ordering, don't change anything, allow reversing
        if "default" in self.sort_by:
            sorted_idx = sorted(self.idx, reverse=self.stats.reverse_sort)
        # Sort alphabetically BUT allow component sorting
        elif "alphabetical" in self.sort_by:
            sorted_idx = self._sort_by_alphabetical()
        else:
            components = self._get_component_order()
            # Azimuthal sorts are allowed to start at a value other than 0
            # Backazimuth needs to come first because 'azimuth' can return both
            if "backazimuth" in self.sort_by:
                sort_list = self.backazimuths - self.sort_by_azimuth_start_deg
            elif "azimuth" in self.sort_by:
                sort_list = self.azimuths - self.sort_by_azimuth_start_deg
            elif "distance" in self.sort_by:
                sort_list = self.distances

            # Sort by the values, AND the components (e.g., Z->R->T or whatever)
            _, _, sorted_idx = \
                zip(*sorted(zip(sort_list, components, self.idx),
                            reverse=self.stats.reverse_sort)
                    )
        sorted_idx = np.array(list(sorted_idx))

        return sorted_idx

    def _sort_by_alphabetical(self):
        """
        Sort by full station name in order. That is, network is sorted, then
        within network station is sorted, and so on. Components are sorted
        last and NOT in alphabetical order. They are ordered based on the
        user-input `components` parameter.

        :rtype: list
        :return: a list of code identifiers to be used to sort the idx list
        """
        networks = [_.split(".")[0] for _ in self.station_ids]
        stations = [_.split(".")[1] for _ in self.station_ids]
        locations = [_.split(".")[2] for _ in self.station_ids]
        components = self._get_component_order()
        zipped = zip(networks, stations, locations, components, self.idx)
        arr_out = np.array(
            list(zip(*sorted(zipped, reverse=self.stats.reverse_sort)))
        )
        # I used the below commented line to check if things were working
        # by converting the component numbers back into componets to make
        # sure the intended 'ZRT' was honored for both alphabetical and reverse
        # arr_out[3] = [complist[int(_)] for _ in arr_out[3]]
        # print(arr_out)

        return arr_out[-1].astype(int)

    def _get_component_order(self):
        """
        When we are sorting, we want components to be sorted by the
        preferred order (i.e, ZRT, Z index is 0, T is 2), which means the
        components have to be reordered to adhere to the sorting algorithm.
        ALSO we need to ensure that we honor component order even if
        everything else is being sorted reverse alphabetically so the
        components entry needs to be the opposite of stats.reverse_sort

        :rtype: list
        :return: components converted to integer values which can be used for
            sorting algorithms.
        """
        channels = [_.split(".")[3] for _ in self.station_ids]
        # ASSUMING component is the last entry in the channel, SEED convention
        # TODO Can we use tr.stats.component here?
        components = [_[-1] for _ in channels]
        if self.stats.reverse_sort:
            comp_list = self.components
        else:
            comp_list = self.components[::-1]

        # Can't list comp here incase `comp_list` does NOT contain one of the
        # available components, which will throw a ValueError. Use a random
        # number to catch these.
        numbered_components = []
        for comp in components:
            try:
                numbered_components.append(comp_list.index(comp))
            except ValueError:
                numbered_components.append(-999)

        return numbered_components

    def get_y_axis_positions(self):
        """
        Determine how seismograms are vertically separated on the Y-Axis when
        plotting. Allows for constant separation, or absolute separation based
        on distance or (back)azimuth separation.

        If `sort_by` has 'absolute' in the argument, then we will weight
        the values defining the Y-Axis. For example, if
        sort_by=='abs_distance', then the Y-Axis will be a scaled
        representation of distances.

        :rtype weights_arr: np.array
        :return: an array of weight values which should be used to place
            seismograms on the Y-Axis
        """
        print(f"determining y-axis positioning for sort: {self.sort_by}")

        # Default weights provides constant `y_axis_spacing` between seismos
        if self.sort_by == "default" or "abs_" not in self.sort_by:
            y_range = np.arange(0, len(self.sorted_idx), 1)
            y_axis = self.y_axis_spacing * y_range
        # Absolute y-axis (i.e., absolute distance or (back)azimuth) will just
        # plot the actual distance or azimuth value, make sure to account for
        # reverse plotting as well
        else:
            if "backazimuth" in self.sort_by:
                y_axis = self.backazimuths
            elif "azimuth" in self.sort_by:
                y_axis = self.azimuths
            elif "distance" in self.sort_by:
                y_axis = self.distances
        return y_axis

    def process_st(self):
        """
        Preprocess the Stream with optional filtering in place.

        .. note::
            Data in memory will be irretrievably altered by running preprocess.

        TODO Add feature to allow list-like periods to individually filter
            seismograms. At the moment we just apply a blanket filter.
        """
        if self.preprocess is None:
            print("no preprocessing applied")
            return
        elif self.preprocess == "st":
            print(f"preprocessing {len(self.st)} `st` waveforms")
            preprocess_list = [self.st]
        elif self.preprocess == "st_syn":
            print(f"preprocessing {len(self.st_syn)} `st_syn` waveforms")
            preprocess_list = [self.st_syn]
        elif self.preprocess == "both":
            print(f"preprocessing {len(self.st) + len(self.st_syn)} "
                  f"`st` and `st_syn` waveforms")
            preprocess_list = [self.st, self.st_syn]

        for st in preprocess_list:
            # Fill any data gaps with mean of the data, do it on a trace by trace
            # basis to get individual mean values
            for tr in st:
                tr.trim(starttime=tr.stats.starttime, endtime=tr.stats.endtime,
                        pad=True, fill_value=tr.data.mean())
            st.detrend("demean")
            st.taper(max_percentage=0.05, type="cosine")

            # Allow multiple filter options based on user input
            # Min period but no max period == low-pass
            if self.max_period_s is not None and self.min_period_s is None:
                print(f"apply lowpass filter w/ cutoff {1/self.max_period_s}")
                st.filter("lowpass", freq=1/self.max_period_s, zerophase=True)
            # Max period but no min period == high-pass
            elif self.min_period_s is not None and self.max_period_s is None:
                print(f"apply highpass filter w/ cutoff {1/self.min_period_s}")
                st.filter("highpass", freq=1/self.min_period_s, zerophase=True)
            # Both min and max period == band-pass
            elif self.min_period_s is not None and \
                    self.max_period_s is not None:
                print(f"applying bandpass filter w/ "
                      f"[{1/self.max_period_s}, {self.min_period_s}]")
                st.filter("bandpass", freqmin=1/self.max_period_s,
                            freqmax=1/self.min_period_s, zerophase=True)
            else:
                print("no filtering applied")

            # Integrate and differentiate N number of times specified by user
            st.detrend("simple")
            if self.integrate != 0:
                if self.integrate < 0:
                    func = "differentiate"
                elif self.integrate > 0:
                    func = "integrate"
            for i in range(np.abs(self.integrate)):
                print("{func} all waveform data")
                getattr(st, func)()

    def plot(self, subset=None, page_num=None, **kwargs):
        """
        Generate record sections based on internal data

        :type subset: list of int
        :param subset: subset of `sorted_idx` if there are too many waveforms to
            plot on one page (set by `max_traces_per_rs`). e.g., to get the
            first 10 entries, subset=[0,10]
        """
        if subset is None:
            start, stop = 0, None  # None will allow full list traversal
            nwav = len(self.sorted_idx)
        else:
            start, stop = subset
            nwav = stop - start

        print(f"plotting record section for {nwav} waveforms")
        print("PLOTTING LINE CHECK STARTING FROM BOTTOM")
        print("\nIDX\tY\t\tID\tDIST\tAZ\tBAZ\tTSHIFT")
        self.f, self.ax = plt.subplots(figsize=self.figsize)

        # Allow choosing observed or synthetic data, defaults to boserved
        for choice in ["st", "st_syn"]:
            if getattr(self, choice) is None:
                continue
            # Main plotting call. Indexes should already be sorted
            # Allow different waveform looks based on observed vs. synthetic
            for y_idx, idx in enumerate(self.sorted_idx[start:stop]):
                # Absolute scaling requires plotting actual dist or az values
                # relative scaling just requires a linear index to stack seismos
                if "abs_" in self.sort_by:
                    y_index = idx
                else:
                    y_index = y_idx + start
                self._plot_trace(idx=idx, y_index=y_index, choice=choice,
                                 **kwargs)

        # Change the aesthetic look of the figure, should be run before other
        # set functions as they may overwrite what is done here
        self._set_plot_aesthetic()

        # if self.sort_by and "azimuth" in self.sort_by:
        #     self._plot_azimuth_bins()
        self._plot_title(nwav=nwav, page_num=page_num)
        if "abs_" in self.sort_by:
            self._set_y_axis_absolute()
            self._set_y_axis_text_labels(start=start, stop=stop, loc="x_min")
        else:
            self._set_y_axis_text_labels(start=start, stop=stop, loc="y_axis")

        # X-axis label is different if we time shift
        if self.time_shift_s.sum() == 0:
            plt.xlabel("Time [s]")
        else:
            plt.xlabel("Relative Time [s]")

        # Allow user defined x-axis limits
        if self.xlim_s is None:
            self.ax.set_xlim([min(self.stats.xmin), max(self.stats.xmax)])
        else:
            self.ax.set_xlim(self.xlim_s)

        # Reverse the y-axis if we are doing absolute y-axis and reversing
        if "abs_" in self.sort_by and "_r" in self.sort_by:
            self.ax.invert_yaxis()

        plt.tight_layout()

        if self.save:
            # Allow appending page numbers to figure names,
            # e.g., default.png -> default_01.png
            if page_num:
                fid, ext = os.path.splitext(self.save)
                save_fid = f"{fid}_{page_num:0>2}{ext}"
            else:
                save_fid = self.save
            print(f"\nsaving figure to {save_fid}")
            plt.savefig(save_fid)
        if self.show:
            plt.show()

    def _plot_trace(self, idx, y_index, choice="st"):
        """
        Plot a single trace on the record section, with amplitude scaling,
        time shifts, etc.s

        Kwargs passed to matplotlib.pyplot.plot()

        :type idx: int
        :param idx: index of the trace to plot and all trace-ordered values like
            amplitude scaling and time shifts
        :type y_index: int
        :param y_index: numerical order which this trace was plotted, as each
            trace is plotted on top of the previous one. y_index should be
            iterated linearly in the loop that calls this function.
        :type choice: str
        :param choice: choice of 'st' or 'st_syn' depending on whether you want
            to plot the observed or synthetic waveforms
        """
        # Used to differentiate the two types of streams for plotting diffs
        choices = ["st", "st_syn"]
        assert (choice in choices)
        c = choices.index(choice)
        tr = getattr(self, choice)[idx]  # i.e., tr = self.st[idx]

        # Plot actual data on with amplitude scaling, time shift, and yoffset
        tshift = self.time_shift_s[idx]
        x = tr.times() + tshift
        y = tr.data / self.amplitude_scaling[idx] + int(self.y_axis[y_index])
        self.ax.plot(x, y, c=["k", "r"][c],  **self.trace_kwargs)

        # Sanity check print station information to check against plot
        print(f"{idx}"
              f"\t{int(self.y_axis[y_index])}"
              f"\t{tr.get_id():<6}"
              f"\t{self.distances[idx]:6.2f}"
              f"\t{self.azimuths[idx]:6.2f}"
              f"\t{self.backazimuths[idx]:6.2f}"
              f"\t{self.time_shift_s[idx]:4.2f}")

        # Retain some stats for global plot args
        self.stats.xmin.append(x.min())
        self.stats.xmax.append(x.max())
        self.stats.ymin.append(y.min())
        self.stats.ymax.append(y.max())

    def _plot_azimuth_bins(self, azimuth_binsize=45):
        """
        If plotting by azimuth, create visual bin separators
        """
        azimuth_bins = np.arange(self.sort_by_azimuth_start_deg,
                                 self.sort_by_azimuth_start_deg + 360,
                                 azimuth_binsize)

        # Make sure that the bins go from 0 <= azimuth_bins <= 360
        azimuth_bins = azimuth_bins % 360
        for azbin in azimuth_bins:
            plt.axhline(y=azbin, c="k", linewidth=2.)

    def _set_y_axis_absolute(self):
        """
        If 'abs_' in sort_by, then the Y-axis should be shown in absolute scale.
        That means we need to format the text labels, add some labelling etc.
        """
        degree_char = u'\N{DEGREE SIGN}'

        # Reset tick label size to be larger to match absolute x-axis size
        ytick_fontsize = self.kwargs.get("ytick_fontsize", 12)
        self.ax.tick_params(axis="y", labelsize=ytick_fontsize)

        if "distance" in self.sort_by:
            ytick_minor = self.kwargs.get("ytick_minor", 25)
            ytick_major = self.kwargs.get("ytick_major", 100)
            ylabel = f"Distance [{self.distance_units}]"
        elif "aziumuth" in self.sort_by:
            ytick_minor = self.kwargs.get("ytick_minor", 45)
            ytick_major = self.kwargs.get("ytick_major", 90)
            ylabel = f"Azimuth [{degree_char}]"

        # Set ytick label major and minor which is either dist or az
        self.ax.yaxis.set_major_locator(MultipleLocator(ytick_major))
        self.ax.yaxis.set_minor_locator(MultipleLocator(ytick_minor))
        self.ax.set_ylabel(ylabel)

    def _set_y_axis_text_labels(self, start=0, stop=-1, loc="y_axis"):
        """
        Plot a text label next to each trace describing the station,
        azimuth and source-receiver distance. We need to do this all at once
        because we have to replace all ticks and tick labels at once.

        .. note::
            if using the 'subset' option in plot, need to also tell the y-axis
            plotter that we are only plotting a subset of data by using the
            `start` and `stop` parameters

        :type start: int
        :param start: starting index for plotting, default to start 0
        :type stop: int
        :param stop: stop index for plotting, default to end -1
        :type loc: str
        :param loc: location to place the y_axis text labels, available:
            - y_axis: Place labels along the y-axis (left side of the figure)
                Will replace the actual y-tick labels so this is probably not
            - x_min: Place labels on the waveforms at the minimum x value
            - x_max: Place labels on the waveforms at the maximum x value
        """
        assert loc in ["y_axis", "x_min", "x_max"]
        degree_char = u'\N{DEGREE SIGN}'

        y_tick_labels = []
        for idx in self.sorted_idx[start:stop]:
            str_id = self.station_ids[idx]
            if self.sort_by is not None and "backazimuth" in self.sort_by:
                # This is named `str_az` but it's actually backazimuths
                str_az = f"{self.backazimuths[idx]:6.2f}{degree_char}"
            else:
                str_az = f"{self.azimuths[idx]:6.2f}{degree_char}"
            str_dist = f"{self.distances[idx]:5.2f}km"

            # Looks like: NN.SSS.LL.CC|30*|250.03km
            label = \
                f"{str_id:>{self.stats.longest_id}}|{str_az:>8}|{str_dist:>8}"
            # Add time shift if we have shifted at all
            if self.time_shift_s[idx] != 0:
                label += f"|{self.time_shift_s[idx]:.2f}s"
            y_tick_labels.append(label)

        if loc == "y_axis":
            # For relative plotting (not abs_), replace y_tick labels with
            # station information
            self.ax.set_yticks(self.y_axis[start:stop])
            self.ax.set_yticklabels(y_tick_labels)
        elif loc == "x_min":
            # Trying to figure out where the minimum X value is on the plot
            if self.xlim_s is not None:
                x_min = min([min(self.xlim_s), min(self.stats.xmin)])
            else:
                x_min = min(self.stats.xmin)

            for idx, s_val in zip(self.sorted_idx[start:stop], y_tick_labels):
                plt.text(x=x_min, y=self.y_axis[idx], s=s_val)

    def _plot_title(self, nwav=None, page_num=None):
        """
        Create the title of the plot based on event and station information
        Allow dynamic creation of title based on user input parameters
        TODO Can we make this two-column to save space?

        :type nwav: int
        :param nwav: if using subset, the title needs to know how many waveforms
            it's showing on the page. self.plot() should tell it
        :type page_num: int
        :param page_num: same as nwav, we need to know what page number were on
        """
        # Defines the number of waveforms plotted on a single page, allowing
        # for subsets per page
        if nwav is None:
            nwav = len(self.sorted_idx)

        # Allow appending page numbers to title
        title_top = "RECORD SECTION"
        if page_num:
            title_top += f" PAGE {page_num:0>2}"

        # The y-label will show baz or az depending on user choice, distinguish
        if self.sort_by is not None and "backazimuth" in self.sort_by:
            az_str = "BAZ"
        else:
            az_str = "AZ"

        # Get the unique components that have been plotted, only
        cmp = "".join(np.unique([self.st[i].stats.component
                                 for i in self.sorted_idx]))

        title = "\n".join([
            title_top,
            f"{'/' * len(title_top*2)}",
            f"ORIGINTIME: {min([tr.stats.starttime for tr in self.st])}",
            f"Y_FMT: NET.STA.LOC.CHA|{az_str}|DIST",
            f"NWAV: {nwav}; NEVT: {self.stats.nevents}; "
            f"NSTA: {self.stats.nstation}; COMP: {cmp}",
            f"SORT_BY: {self.sort_by}; SCALE_AMP: {self.scale_amp}",
            f"FILT: [{self.min_period_s}, {self.max_period_s}]s; "
            f"MOVE_OUT: {self.move_out or 0}{self.distance_units}/s",
        ])
        self.ax.set_title(title)

    def _set_plot_aesthetic(self):
        """
        Give a nice look to the output figure by creating thicc borders on the
        axis, adjusting fontsize etc. All plot aesthetics should be placed here
        so it's easiest to find.

        .. note::
            This was copy-pasted from Pyatoa.visuals.insp_plot.default_axes()
        """
        ytick_fontsize = self.kwargs.get("ytick_fontsize", 8)
        xtick_fontsize = self.kwargs.get("xtick_fontsize", 12)
        tick_linewidth = self.kwargs.get("tick_linewidth", 1.5)
        tick_length = self.kwargs.get("tick_length", 5)
        tick_direction = self.kwargs.get("tick_direction", "in")
        label_fontsize = self.kwargs.get("label_fontsize", 10)
        axis_linewidth = self.kwargs.get("axis_linewidth", 2.)
        title_fontsize = self.kwargs.get("title_fontsize", 10)
        xtick_minor = self.kwargs.get("xtick_minor", 25)
        xtick_major = self.kwargs.get("xtick_major", 100)

        # Re-set font sizes for labels already created
        self.ax.title.set_fontsize(title_fontsize)
        self.ax.tick_params(axis="both", which="both", width=tick_linewidth,
                            direction=tick_direction, length=tick_length)
        self.ax.tick_params(axis="x", labelsize=xtick_fontsize)
        self.ax.tick_params(axis="y", labelsize=ytick_fontsize)
        self.ax.xaxis.label.set_size(label_fontsize)

        # Thicken up the bounding axis lines
        for axis in ["top", "bottom", "left", "right"]:
            self.ax.spines[axis].set_linewidth(axis_linewidth)

        # Set xtick label major and minor which is assumed to be a time series
        self.ax.xaxis.set_major_locator(MultipleLocator(xtick_major))
        self.ax.xaxis.set_minor_locator(MultipleLocator(xtick_minor))

        plt.grid(visible=True, which="major", axis="x", alpha=0.5, linewidth=1)
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)


def plotw_rs(*args, **kwargs):
    """
    Main. Run the record section plotting functions in order.
    """
    _start = datetime.now()
    print(f"STARTING RECORD SECTION PLOTTER")

    rs = RecordSection(*args, **kwargs)
    rs.process_st()
    rs.get_parameters()
    # Simple case where all waveforms will fit on one page
    if len(rs.st) <= rs.max_traces_per_rs:
        rs.plot()
    # More complicated case where we need to split onto multiple pages
    else:
        for i, start in enumerate(np.arange(0, len(rs.st),
                                            rs.max_traces_per_rs)):
            stop = start + rs.max_traces_per_rs
            # Case where the num waveforms is less than max_traces_per_rs
            if stop < rs.max_traces_per_rs:
                stop = len(rs.st)
            rs.plot(subset=[start, stop], page_num=i+1)

    _end = datetime.now()
    print(f"FINISHED RECORD SECTION (t={_end - _start}s)")


def testw_rs():
    """
    Testing function for dev purposes, to be deleted when published
    """
    from glob import glob
    from obspy import read, Stream

    st = Stream()
    path = "/Users/Chow/Repositories/pysep/20090407201255351"
    for fid in glob(os.path.join(path, "20090407201255351.??.*.*.*.?"))[:50]:
        st += read(fid)

    plotw_rs(
        st,
        st_syn=None,
        # sort_by="distance",
        sort_by="abs_distance_r",
        scale_amp="normalize",
        # scale_amp=None,
        time_shift_s=None,
        time_marker=None,
        min_period_s=2,
        max_period_s=50,
        # move_out=4,
        preprocess="st",
        max_traces_per_rs=50,
        integrate=0,
        components="Z",
        # components="ZRTNE12",
        #xlim_s=[50, 250],
        y_axis_spacing=1,
        sort_by_azimuth_start_deg=0,
        distance_units="km",
        geometric_spreading_factor=0.5,
        geometric_spreading_k_val=None,
        figsize=(9, 11),
        show=False,
        save="./output/record_section.png",
        trace_kwargs={
            "linewidth": 1.
        },
        overwrite=True,
        # Kwargs below
        xtick_minor=25,
        xtick_major=50,
    )


def create_examples():
    """
    Create example figures with a variety of input parameters to test out the
    functionality of the plotting script

    Tests:
        sort_by, xlim, filter, comp, multipage, distance_units
    """
    from glob import glob
    from obspy import read, Stream

    st = Stream()
    path = "/Users/Chow/Repositories/pysep/20090407201255351"
    for fid in glob(os.path.join(path, "20090407201255351.??.*.*.*.?"))[:50]:
        st += read(fid)

    for sort_by in ["default", "default_r", "alphabetical", "alphabetical_r",
                    "distance", "distance_r", "azimuth", "azimuth_r",
                    "backazimuth", "backazimuth_r"]:
        plotw_rs(st, sort_by=sort_by, scale_amp="normalize", min_period_s=2,
                 max_period_s=50, preprocess="st", max_traces_per_rs=50,
                 components="ZRTNE12", y_axis_spacing=1, distance_units="km",
                 figsize=(9, 11), show=False, save=f"./output/rs_{sort_by}.png",
                 trace_kwargs={"linewidth": 1.}, overwrite=True,)

    plotw_rs(st, sort_by="distance", scale_amp="normalize", min_period_s=2,
             max_period_s=50, preprocess="st", max_traces_per_rs=50,
             move_out=4, y_axis_spacing=1, distance_units="km",
             figsize=(9, 11), show=False, save=f"./output/rs_move_out.png",
             trace_kwargs={"linewidth": 1.}, overwrite=True, )


if __name__ == "__main__":
    # create_examples()
    testw_rs()
    # plotw_rs()

