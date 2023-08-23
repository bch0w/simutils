"""
Utilities for reading various file types, mostly from Specfem3D to ObsPy classes
These are meant to be standalone functions so they may repeat some functionality
found elsewhere in the package.
"""
import os
import numpy as np
from glob import glob

from obspy import Stream, Catalog, read, read_events
from pysep.utils.io import read_specfem2d_source, read_forcesolution
from pyatoa import logger
from pyatoa.utils.calculate import overlapping_days


def read_waveforms_from_seed_directory(
        code, origin_time, base_path="./",
        obs_dir_template="{year}/{net}/{sta}/{cha}",
        obs_fid_template="{net}.{sta}.{loc}.{cha}.{year}.{jday:0>3}",
        start_pad=3600, end_pad=3600):
    """
    Fetch seismic data (waveforms) via a very specific directory structure that
    is typically used in SEED datacenters, where data are stored as 24-hour
    MSEED files, and organized by their network, station, channel and day.

    .. note::

        Default waveform directory structure assumed to follow SEED
        convention. That is:
        base_path/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
        e.g. base_path/2017/NZ/OPRZ/HHZ.D/NZ.OPRZ.10.HHZ.D

    :type code: str
    :param code: Station code following SEED naming convention.
        This must be in the form NN.SSSS.LL.CCC (N=network, S=station,
        L=location, C=channel). Allows for wildcard naming. By default
        the pyatoa workflow wants three orthogonal components in the N/E/Z
        coordinate system. Example station code: NZ.OPRZ.10.HH?
    :type origin_time: UTCDateTime
    :param origin_time: the origin time of the event or waveform used to
        determine which files to read from. Parameters `start_pad` and `end_pad`
        are used to set a buffer time region around the `origin_time` incase
        waveforms are requested across multiple days.
    :type base_path: str
    :param base_path: the base path where the MSEED directories are presumed
        to start, and where the sub-directories will be built from (see note
        above). Defaults to CWD
    :type obs_dir_template: str
    :param obs_dir_template: directory structure to search for observation
        data. Follows the SEED convention:
        'path/to/obs_data/{year}/{net}/{sta}/{cha}'
    :type obs_fid_template: str
    :param obs_fid_template: File naming template to search for observation
        data. Follows the SEED convention:
        '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
    :type start_pad: int
    :param start_pad: buffer time BEFORE `origin_time` in units of s. Defaults
        to 3600s (1h)
    :type end_pad: int
    :param end_pad: buffer time AFTER `origin_time` in units of s. Defaults to
        3600s (1h)
    :rtype stream: obspy.core.stream.Stream or None
    :return stream: stream object containing relevant waveforms, else None
    """
    net, sta, loc, cha = code.split('.')

    # If waveforms contain midnight, multiple files need to be read.
    # Checks `start_pad` before and `end_pad` after origintime
    jdays = overlapping_days(origin_time=origin_time, start_pad=start_pad,
                             end_pad=end_pad)

    st = Stream()
    pathlist = []
    full_path = os.path.join(base_path, obs_dir_template, obs_fid_template)
    for jday in jdays:
        pathlist.append(full_path.format(net=net, sta=sta, cha=cha,
                                         loc=loc, jday=jday,
                                         year=origin_time.year)
                        )
    for fid in pathlist:
        logger.debug(f"searching for observations: {fid}")
        for filepath in glob(fid):
            st += read(filepath)
            logger.info(f"retrieved observations locally: {filepath}")

    # Take care of gaps in data by converting to masked data
    if len(st) > 0:
        st.merge()
    else:
        logger.warning("No waveform data found for the given SEED "
                       f"configurations: {pathlist}")

    return st


