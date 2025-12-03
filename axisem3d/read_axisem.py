"""
Read the ASCII outputs from AxiSEM as an ObsPy Stream
"""
import os
import numpy as np
from obspy import Stream, Trace, UTCDateTime


def channel_code(dt):
    """
    Specfem outputs seismograms with channel band codes set by IRIS. Instrument
    codes are always X for synthetics, but band code will vary with the sampling
    rate of the data, return the correct code given a sampling rate.
    Taken from Appenix B of the Specfem3D cartesian manual (June 15, 2018)

    :type dt: float
    :param dt: sampling rate of the data in seconds
    :rtype: str
    :return: band code as specified by Iris
    :raises KeyError: when dt is specified incorrectly
    """
    if dt >= 1:
        return "L"  # long period
    elif 0.1 < dt < 1:
        return "M"  # mid period
    elif 0.0125 < dt <= 0.1:
        return "B"  # broad band
    elif 0.001 <= dt <= 0.0125:
        return "H"  # high broad band
    elif 0.004 <= dt < 0.001:
        return "C"
    elif 0.001 <= dt < 0.0002:
        return "F"
    else:
        raise KeyError("Channel code does not exist for this value of 'dt'")


def read_axisem(fid, time_fid="data_time.ascii", 
                origintime="1970-01-0100:00:00", location="", components="RTZ"):
    """
    Read AxiSEM synthetics as ObsPy Stream
    """
    data = np.loadtxt(fid).T
    times = np.loadtxt(time_fid)

    origintime = UTCDateTime(origintime)

    # gather metadata for the stats header
    delta = times[1] - times[0]  # s
    net, sta, ext = os.path.basename(fid).split(".")

    # data are stored as multi-component
    data_dict = {}
    traces = []
    for i, comp in enumerate(components):
        tr_data = data[i]

        # build channel code
        cha = f"{channel_code(delta)}X{comp}"
        stats = {"network": net, "station": sta, "location": location,
                 "channel": cha, "starttime": origintime, "npts": len(data[0]),
                 "delta": delta, "mseed": {"dataquality": 'D'}, 
                 }
        traces.append(Trace(data=tr_data, header=stats))

    return(Stream(traces))


if __name__ == "__main__":
    import sys
    st = read_axisem(f"VN.S{sys.argv[1]}.ascii")
    st.filter("lowpass", freq=1/60)
    st.plot()
