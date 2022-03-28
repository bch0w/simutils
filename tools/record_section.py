#!/usr/bin/evn python3
"""
Record Section plotting tool for seismic waveforms (observed and synthetic)

This is a refactor of Pysep's Python utility `plotw_rs`, a record section
plotting script. The intent of this script is to plot multiple time series'
based on source-receiver characteristics (i.e., src-rcv distance, backazimuth).

.. note:: Code History
    Written by Carl Tape (11/21/2011) and Yun Wang (11/2011) as a Matlab script
    Translated to Python by Nealy Sims (1/2021)
    Upgraded by Aakash Gupta (9/2021)
    Refactored by Bryant Chow (3/2022)
"""
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth,


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
    def __init__(self, st, st_syn=None, sort_by=None, scale_amp=None,
                 time_shift_s=None, time_marker=None, min_period_s=None,
                 max_period_s=None, max_traces_per_rs=50,
                 integrate=0, differentiate=0, xlim_s=None,
                 sort_by_weight_factor=None, sort_by_azimuth_start_deg=None,
                 distance_units="km", geometric_spreading_factor=0.5,
                 geometric_spreading_k_val=None, plot_map=False, show=True):
        """
        Set the default record section plotting parameters and enforce some types.

        :type st: obspy.core.stream.Stream
        :param st: Stream objects containing observed time series to be plotted on
            the record section. Can contain any number of traces
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream objects containing synthetic time series to be plotted
            on the record section. Must contain the same number of traces as `st`
        :type sort_by: str
        :param sort_by: How to sort the Y-axis of the record section, available:
            None: Not set, don't sort, just iterate directly through the Stream
            'azimuth': sort by source-receiver azimuth (deg) with constant
                vertical spacing
            'backazimuth': sort by source-receiver backazimuth (deg) with constant
                vertical spacing. Requires `sort_by_azimuth_start_deg`
            'distance': sort by source-receiver distance (km) with constant
                vertical spacing. Requires `sort_by_azimuth_start_deg` AND
                `distance_units`
            'weighted_distance': weight vertical spacing of waveforms by the
                source-receiver distance (km). Requires `distance_units`
            'weighted_azimuth': weight vertical spacing of waveforms by
                source-receiver azimuth (deg). Requires `sort_by_weight_factor` AND
                `sort_by_azimuth_start_deg`
            'weighted_backazimuth': weight vertical spacing of waveforms by
                source-receiver backazimuth (deg). Requires `sort_by_weight_factor`
        :type scale_amp: list
        :param scale_amp: scale amplitude of waveforms by available:
            None: Not set, no amplitude scaling, waveforms shown raw
            'normalize': scale each trace by the maximum amplitude,
                i.e., > a /= max(abs(a))  # where a is the time series amplitudes
            'geometric_spreading': scale amplitudes globally by predicting the
                expected geometric spreading amplitude reduction and correcting
                for this factor. Requires `geometric_spreading_factor`, optional
                `geometric_spreading_k_val`
        :type time_shift_s: float OR list of float
        :param time_shift_s: apply a static time shift to waveforms, two options:
            1. float (e.g., -10.2), will shift ALL waveforms by
                that number (i.e., -10.2 second time shift applied)
            2. list (e.g., [5., -2., ... 11.2]), will apply individual time shifts
                to EACH trace in the stream. The length of this list MUST match
                the number of traces in your input stream.
        :type time_marker: list of float
        :param time_marker: absolute time markers. If this parameter is set,
            `time_shift_s` MUST be set as Option 1 (single float)
        :type min_period_s: float
        :param min_period_s: minimum filter period in seconds
        :type max_period_s: float
        :param max_period_s: maximum filter period in seconds
        :type max_traces_per_rs: int
        :param max_traces_per_rs: maximum number of traces to show on a single
            record section plot
        :type integrate: int
        :param integrate: apply integration `integrate` times on all traces.
            i.e., if integrate == 2, will integrate each trace twice.
        :type differentiate: int
        :param differentiate: apply differentiation `differentiate` times on all
            traces. i.e., if differentiate == 1, will differentiate each trace
            one time.
        :type sort_by_weight_factor: float
        :param sort_by_weight_factor: a trial-and-error weight factor applied to
            Y-axis sorting which have 'weighted' in their name. Use smaller values
            for more prominent peaks (TODO Check this statement)
            Set to 1 for default behavior
        :type sort_by_azimuth_start_deg: float
        :param sort_by_azimuth_start_deg: azimuthal angle for the top record,
            applied to any Y-axis sorting which has 'azimuth' in its name.
            Set to 0 for default behavior
        :type distance_units: str
        :param distance_units: Y-axis units for epicentral distance record sections
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
        :type plot_map: bool
        :param plot_map: plot a source-receiver map alongside the record section
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
        self.sort_by_weight_factor = float(sort_by_weight_factor)
        self.sort_by_azimuth_start_deg = float(sort_by_azimuth_start_deg)
        # Amplitude scaling parameters
        self.scale_amp = scale_amp
        self.geometric_spreading_factor = float(geometric_spreading_factor)
        self.geometric_spreading_k_val = float(geometric_spreading_k_val)
        # Time shift parameters
        self.time_shift_s = time_shift_s
        self.time_marker = time_marker
        # Filtering parameters
        self.min_period_s = float(min_period_s)
        self.max_period_s = float(max_period_s)
        self.max_traces_per_rs = int(max_traces_per_rs)
        self.integrate = int(integrate)
        self.differentiate = int(differentiate)
        # Plotting parameters
        self.xlim_s = xlim_s
        self.distance_units = distance_units.lower()
        self.plot_map = bool(plot_map)
        self.show = bool(show)

        # Run checks to ensure that all the parameters are set properly
        self.check_parameters()

        # Internally derive parameters used to keep track of srcrcv info
        self.idx = range(0, len(self.st), 1)
        self.sorted_idx = self.get_sorted_idx()
        self.amplitude_scaling = self.get_amplitude_scaling()

        self.distances, self.azimuths, self.backazimuths = [], [], []
        for tr in self.st:
            _gcd, _az, _baz = self.get_srcrcv_info(tr=tr)
            self.distances.append(_gcd)
            self.azimuths.append(_az)
            self.backazimuths.append(_baz)

    def check_parameters(self):
        """
        Check that parameters are set properly and in line with how they
        are expected by the program

        .. note::
            Not using assertions here because we want all incorrect parameters to be
            evaluated together and displayed at once, that way the user doesn't
            have to run this function multiple times to figure out how to set their
            parameters correct
        """
        err = Dict()

        # Pysep should have created SAC headers if we want to sort by dist or az
        if self.sort_by is not None and \
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

        if self.sort_by is not None:
            acceptable_sort_by = ["azimuth", "backazimuth", "distance",
                                  "alphabetical", "weighted_azimuth",
                                  "weighted_distance"]
            acceptable_distance_units = ["km", "km_utm", "deg"]
            if self.sort_by not in acceptable_sort_by:
                err.sort_by = f"must be in {acceptable_sort_by}"

            if ("weighted" in self.sort_by) and \
                    self.sort_by_weight_factor is None:
                err.sort_by_weight_factor = f"must be set for weighted sort"

            if ("azimuth" in self.sort_by) and \
                    (self.sort_by_azimuth_start_deg is None):
                err.sort_by_azimuth_start_deg = f"must be set for azimuth sort"
            else:
                if not 0 <= self.sort_by_azimuth_start_deg <= 360:
                    err.sort_by_azimuth_start_deg = f"0 < azi < 360"

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


        if self.xlim_s is not None:
            if len(self.xlim_s) != 2:
                err.xlim_s = f"must be of length 2, [start, stop]"
            elif self.xlim_s[0] > self.xlim_s[1]:
                err.xlim_s = f"start time must be less than stop time"

        if err:
            out = "Parameter errors found, please make the following changes:"
            out += "\n".join([f"{key}: {val}" for key, val in err.items()])
            raise AssertionError(out)

    def get_srcrcv_info(self, tr=None, idx=0):
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
            gcdist = tr.stats.sac.dist
            az = tr.stats.sac.az
            baz = tr.stats.sac.baz
        except AttributeError:
            # If for whatever reason SAC headers dont contain this info already
            # Use ObsPy to get great circle distance, azimuth and backazimuth
            gcdist, az, baz = gps2dist_azimuth(lat1=tr.stats.sac.evla,
                                               lon1=tr.stats.sac.evlo,
                                               lat2=tr.stats.sac.stla,
                                               lon2=tr.stats.sac.stlo)
        # Don't allow negative degree values
        if az < 0:
            az += 360
        if baz < 0:
            baz += 360

        # Make sure distance is in the correct units, default units 'km'
        gcdist *= 1E-3  # units: m -> km
        if self.distance_units == "deg":
            gcdist = kilometers2degrees(gcdist)  # units: km -> deg
        elif self.distance_units == "km_utm":
            # Overwrite `gcdist`, could probably skip that calc above but
            # leaving for now as I don't think this option will be used heavily.
            gcdist = np.sqrt(
                ((tr.stats.sac.stlo - tr.stats.sac.evlo) ** 2) /
                ((tr.stats.sac.stla - tr.stats.sac.evla) ** 2)
            )
            gcdist *= 1E-3  # units: m -> km

        return gcdist, az, baz

    def get_sorted_idx(self, reverse=False):
        """
        Sort the source-receiver pairs by the chosen parameters.
        Overwrites the internal `idx` parameter, which defines the order
        in which the Stream traces should be organized

        TODO allow the user to set sort arguments of alphabetical
        """
        if self.sort_by is None:
            sorted_idx = self.idx
        else:
            if self.sort_by == "alphabetical":
                sort_list = self._sort_by_alphabetical(net_or_sta="net",
                                                       reverse=False)
            elif self.sort_by == "azimuth":
                sort_list = self.azimiths
            elif self.sort_by == "backazimuth":
                sort_list = self.backazimuths
            elif self.sort_by == "distance":
                sort_list = self.distances
            sorted_tup = zip(*sorted(zip(sort_list, self.idx), reverse=reverse))
            sorted_idx = list(sorted_tup[1])

        return sorted_idx

    def _sort_by_alphabetical(self, net_or_sta="net"):
        """
        Sort by full station name, allow choosing sort based on network or
        station. The expected station code looks like NN.SSS.LL.CCC
        where N=network, S=station, L=location, C=channel

        :type net_or_sta: str
        :param net_or_sta: sort by network 'net', or station 'sta'
        :rtype: list
        :return: a list of code identifiers to be used to sort the idx list
        """
        names = [tr.get_id() for tr in self.st]
        if net_or_sta == "net":
            sort_by = [_.split(".")[0] for _ in names]  # part 0 is network
        elif net_or_sta == "sta":
            sort_by = [_.split(".")[1] for _ in names]  # part 1 is station
        return sort_by

    def preprocess(self, st_input=None):
        """
        Preprocess the Stream with optional filtering

        TODO Add feature to allow list-like periods to individually filter
            seismograms. At the moment we just apply a blanket filter.
        """
        if st_input is None:
            st = self.st.copy()

        # Fill any data gaps with mean of the data, do it on a trace by trace
        # basis to get individual mean values
        for tr in st:
            tr.trim(starttime=tr.stats.starttime, endtime=tr.stats.endtime,
                    pad=True, fill_vale=tr.data.mean())
        st.detrend("demean")
        st.taper(max_percentage=0.05, type="cosine")

        # Allow multiple filter options based on user input
        # Min period but no max period == low-pass
        if self.max_period_s is not None and self.min_period_s is None:
            print(f"applying lowpass filter w/ cutoff {1/self.max_period_s}")
            st.filter("lowpass", freq=1/self.max_period_s, zerophase=True)
        # Max period but no min period == high-pass
        elif self.min_period_s is not None and self.max_period_s is None:
            print(f"applying highpass filter w/ cutoff {1/self.min_period_s}")
            st.filter("highpass", freq=1/self.min_period_s, zerophase=True)
        # Both min and max period == band-pass
        elif self.min_period_s is not None and self.max_period_s is not None:
            print(f"applying bandpass filter w/ "
                  f"[{1/self.max_period_s}, {self.min_period_s}]")
            st.filter("bandpass", freqmin=1/self.max_period_s,
                        freqmax=1/self.min_period_s, zerophase=True)
        else:
            print("no filtering applied")

        # Integrate and differentiate N number of times specified by user
        # this is a "dumb" approach so you can screw up your data by putting
        # values for both integrate and differentiate.
        st.detrend("simple")
        if self.integrate:
            for i in range(self.integrate):
                print("integrating all waveform data")
                st.integrate()
        if self.differentiate:
            for i in range(self.differentiate):
                print("differentiating all waveform data")
                st.differentiate()

        return st

    def get_amplitude_scaling(self):
        """
        Scale the amplitudes of all the waveforms by producing a Stream
        dependent scale factor based on user choice. It is expected that the
        output array will be MULTIPLIED by the data arrays

        .. note::
            Needs to be run AFTER preprocessing because filtering etc. will
            affect the final amplitudes of the waveforms
        """
        max_amplitudes = [1/max(abs(tr.data)) for tr in self.st]
        if self.scale_amp is None:
            amp_scaling = np.ones(len(self.st))
        elif self.scale_amp == "normalize":
            amp_scaling = max_amplitudes
        elif self.scale_amp == "geometric_spreading":
            amp_scaling = self._calculate_geometric_spreading(max_amplitudes)

        return amp_scaling

    def _calculate_geometric_spreading(self, max_amplitudes=None):
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
        if max_amplitudes is None:
            max_amplitudes = [1 / max(abs(tr.data)) for tr in self.st]

        # Create a sinusoidal function based on distances in degrees
        # !!! TODO what does "sindel" mean?
        sindel = np.sin(np.array(self.dist) / (180 / np.pi))

        # !!! TODO look at Stein and Wysession and figure out what these vector
        # !!! TODO names are, sort of ambiguous right meow
        if self.geometric_spreading_k_val is not None:
            k_vector = max_amplitudes * \
                       (sindel ** self.geometric_spreading_factor)
            k_val = np.median(k_vector)
        else:
            k_val = self.geometric_spreading_k_val

        w_vector = k_val / (sindel ** self.geometric_spreading_factor)
        inverse_w = 1/ w_vector

        return inverse_w

    def plot(self, azimuth_binsize=45):
        """
        Generate record sections based on internal data
        """
        azimuth_bins = np.arange(self.sort_by_azimuth_start_deg,
                                 self.sort_by_azimuth_start_deg + 360,
                                 azimuth_binsize)
        # Make sure that the bins go from 0 <= azimuth_bins <= 360
        azimuth_bins = azimuth_bins % 360

        


def plotw_rs(*args, **kwargs):
    """
    Main. Run the record section plotting functions in order.
    """
    _start = datetime.now()

    rs = RecordSection(*args, **kwargs)

    print(f"starting plotw_rs to create record sections")

    _end = datetime.now()
    print(f"completed record section plotting in {_start - _end}s")

if __name__ == "__main__":
    plotw_rs()

