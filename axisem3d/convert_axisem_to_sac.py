"""
Convert AxiSEM synthetics to SAC files with proper metadata so that they can
be plotted with RecSec
"""
import os
import sys
import numpy as np
import re
import yaml
from glob import glob
from obspy import UTCDateTime
from pysep import RecordSection
from obspy import Stream, Catalog, Trace
from obspy.core.event import Event, Origin, Magnitude
from pysep import read_stations
from pysep.utils.cap_sac import append_sac_headers
from pysep.utils.fmt import channel_code
from pysep.utils.mt import seismic_moment, moment_magnitude

# Dummy variable
ORIGINTIME = "2000-01-01T00:00:00"


def _get_resource_id(name, res_type, tag=None):
    """
    Helper function to create consistent resource ids, from ObsPy. Used to
    create resource ID's when generating Event objects
    """
    res_id = f"smi:local/source/{name:s}/{res_type:s}"
    if tag is not None:
        res_id += "#" + tag
    return res_id


def read_axisem(fid, time_fid="data_time.ascii",
                origintime=ORIGINTIME, location="", components="RTZ"):
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


def read_axisource(fid="./inparam.source.yaml"):
    """
    Read source information from the `inparam.source.yaml` file. This is
    slightly hardcoded for a specific parameter set so it is not generalizable,
    use at your own risk.
    """
    yaml.SafeLoader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.')
        )

    with open(fid, "r") as f:
        srcdict = yaml.safe_load(f)

    cat = Catalog()
    
    # Multiple sources allowed in the input source file
    for srcs in srcdict["list_of_sources"]:
        for name, src in srcs.items():
            print(name)

            # Get geographic information
            loc = src["location"]
            latitude, longitude = loc["latitude_longitude"]
            origin_time = UTCDateTime(ORIGINTIME)
            assert(loc["depth_below_solid_surface"]), "required"
            depth = float(loc["depth"]) * 1E3  # units: m

            origin = Origin(
                resource_id=_get_resource_id(name, "origin", tag="source"),
                time=origin_time, longitude=longitude, latitude=latitude,
                depth=depth  
            )

            # Get magnitude information
            mec = src["mechanism"]
            mt = [float(_) for _ in mec["data"]]  # dyn*cm
            m0 = seismic_moment(mt)  
            mw = moment_magnitude(m0)

            
            magnitude = Magnitude(
                resource_id=_get_resource_id(name, "magnitude"),
                mag=mw, magnitude_type="Mw", origin_id=origin.resource_id.id
            )

            event = Event(resource_id=_get_resource_id(name=name,
                                                       res_type="event"))
            event.origins.append(origin)
            event.magnitudes.append(magnitude)

            event.preferred_origin_id = origin.resource_id.id
            event.preferred_magnitude_id = magnitude.resource_id.id

            cat.append(event)

    return cat


if __name__ == "__main__":
    path = sys.argv[1]
    AXISEM = glob(f"{path}/ascii/VN.*.ascii")
    TIME   = f"{path}/ascii/data_time.ascii"
    SRCFID = "inparam.source.yaml"
    STAFID = "STATIONS_GRID"
    OUT    = path

    cat = read_axisource(SRCFID)
    assert(len(cat) == 1)  # currently only 1 src allowed
    event = cat[0]
    inv = read_stations(STAFID)

    for fid in AXISEM:
        print(fid)
        st = read_axisem(fid, time_fid=TIME)  # 3-C
        st = append_sac_headers(st, event, inv)
        for tr in st:
            tr.write(f"{OUT}/{tr.get_id()}.SAC", format="SAC")

