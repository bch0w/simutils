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
import sys
import numpy as np
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
                 plot_map=False, show=True):
        """
        Set the default record section plotting parameters and do some simple checks
        to make sure that the input parameters are good to go.

        .. note::
            Not using assertions here because we want all incorrect parameters to be
            evaluated together and displayed at once, that way the user doesn't
            have to run this function multiple times to figure out how to set their
            parameters correct

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
                for this factor. Requires `geometric_spreading_factor`
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
            'km_sphere': kilometers on the sphere
            'deg_sphere': degrees on the sphere
            'km_cartesian': kilometers on flat plane
        :type geometric_spreading_factor: float
        :param geometric_spreading_factor: geometrical spreading factor when
            using the `scale_amp` parameter. Defaults to 0.5 for surface waves.
            Use values of 0.5 to 1.0 for regional surface waves
        :type plot_map: bool
        :param plot_map: plot a source-receiver map alongside the record section
        :rtype: Dict
        :return: a Dictionary object containing all the user set parameters
        :raises AssertionError: if any parameters are set incorrectly
        """
        # To be used to state how input parameters are set incorrectly
        err_msgs = Dict()

        # No type checking, assuming that the User knows these are streams
        assert st, f"Stream input must contain at least 1 (one) trace!"
        self.st = st.copy()
        try:
            self.st_syn = st_syn.copy()
            if len(st) != len(st_syn):
                err_msgs.st_syn = f"length must match `st` (which is {len(st)})"
        except AttributeError:
            self.st_syn = None

        # Set and check all the sorting parameters and their required sub arguments
        self.sort_by = sort_by
        self.sort_by_weight_factor = sort_by_weight_factor
        self.sort_by_azimuth_start_deg = sort_by_azimuth_start_deg
        self.distance_units = distance_units
        if sort_by is not None:
            acceptable_sort_by = ["azimuth", "backazimuth", "distance",
                                  "alphabetical", "weighted_azimuth",
                                  "weighted_distance"]
            acceptable_distance_units = ["km_sphere", "km_flat", "deg_sphere"]
            if sort_by.lower() not in acceptable_sort_by:
                err_msgs.sort_by = f"must be in {acceptable_sort_by}"

            if ("weighted" in sort_by.lower()) and \
                    sort_by_weight_factor is None:
                err_msgs.sort_by_weight_factor = \
                    f"must be set for weighted sort"

            if ("azimuth" in sort_by.lower()) and \
                    (sort_by_azimuth_start_deg is None):
                err_msgs.sort_by_azimuth_start_deg = \
                    f"must be set for azimuth sort"
            else:
                if not 0 <= sort_by_azimuth_start_deg <= 360:
                    err_msgs.sort_by_azimuth_start_deg = f"0 < azi < 360"

            if ("distance" in sort_by.lower()) and \
                    (distance_units not in acceptable_distance_units):
                err_msgs.sort_by_azimuth_start_deg = \
                    f"must be in {acceptable_distance_units}"

        self.scale_amp = scale_amp
        self.geometric_spreading_factor = geometric_spreading_factor
        if scale_amp is not None:
            acceptable_scale_amp = ["normalize", "geometric_spreading"]
            if scale_amp.lower() not in acceptable_scale_amp:
                err_msgs.scale_amp = f"must be in {acceptable_scale_amp}"

        self.time_shift_s = time_shift_s
        if time_shift_s is not None:
            acceptable_time_shift_s = [1, len(st)]
            if len(time_shift_s) not in acceptable_time_shift_s:
                err_msgs.time_shift_s = f"must be in {acceptable_time_shift_s}"

        self.time_marker = time_marker
        if time_marker is not None:
            if (time_shift_s is not None) and (time_shift_s != 1):
                err_msgs.time_marker = ("time_shift_s cannot be variable when "
                                        "using this option")

        self.min_period_s = min_period_s
        self.max_period_s = max_period_s
        if min_period_s is not None and max_period_s is not None:
            if min_period_s >= max_period_s:
                err_msgs.min_period_s = "must be less than `max_period_s`"

        self.max_traces_per_rs = max_traces_per_rs
        if not isinstance(max_traces_per_rs, int):
            err_msgs.max_traces_per_rs = f"must be of type integer"

        self.integrate = integrate
        if not isinstance(integrate, int):
            err_msgs.integrate = f"must be of type integer"

        self.differentiate = differentiate
        if not isinstance(differentiate, int):
            err_msgs.differentiate = f"must be of type integer"

        self.xlim_s = xlim_s
        if xlim_s is not None:
            if len(xlim_s) != 2:
                err_msgs.xlim_s = f"must be of length 2, [start, stop]"
            elif xlim_s[0] > xlim_s[1]:
                err_msgs.xlim_s = f"start time must be less than stop time"

        self.plot_map = plot_map
        self.show = show

        # Check if any error messages were created, if so crash out gracefully
        if err_msgs:
            out = "Parameter errors found, please make the following changes:"
            out += "\n".join([f"{key}: {val}" for key, val in err_msgs.items()])
            raise AssertionError(out)

    def get_srcrcv_characteristics(self, tr=None, event=None, idx=0):
        """
        Check the source-receiver characteristics such as src-rcv distance,
        azimuth, backazimuth for a given trace.
        """
        # Default to the first trace in the Stream
        if tr is None:
            tr = self.st[0]

        # Pysep should have created SAC headers
        if hasattr(tr.stats, "sac"):
            sta_lat = tr.stats.sac.stla
            sta_lon = tr.stats.sca.stlo
        elif hasattr(tr.stats, "coordinates"):
            sta_lat = tr.stats.coordinates[0]
            sta_lon = tr.stats.coordinates[1]
        else:
            raise AttributeError(f"Trace {tr.get_id}) (idx = {idx} has no "
                                 f"'sac' or 'coordinates' attributes, cannot "
                                 f"sort by distance or azimuth")

        # Use ObsPy to get the great circle distance, azimuth and backazimuth
        gcdist, az, baz = gps2dist_azimuth(
            lat1=event.preferred_origin().latitude,
            lon1=event.preferred_origin().longitude,
            lat2=sta_lat, lon2=sta_lon
        )

        return gcdist, az, baz






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

