"""
Simple script to run the Decluster class for declustering ObsPy catalogs
"""
from pysep import Declust, read_events_plus, read_stations

# User select parameters 
nx = 10  # number of bins in the X-direction, larger numbers = more events
nx = 10  # number of bins in the Y-direction


cat = read_events_plus("events.xml")
dc = Declust(cat=cat, use_magnitudes=True, use_depths=True)

# Remove events that fall below a given magnitud range
dc.threshold_catalog(zedges=[0, 10, 100],  # depth boundaries for vertical bins
                     min_mags=[4, 5],      # min magnitude in each depth bin
                     )

# Decluster the Catalog, larger events take priority
cat = dc.decluster_events(nx=nx, ny=ny, select_by="magnitude", 
                          plot=True, plot_dir="./figs", dpi=200)

# Get event labels 
event_ids = []
for event in cat:
    eid = event.resource_id.id.replace(".", "_")
    eid = f"CMTSOLUTION_{eid}"
    event_ids.append(eid)

# Write out list of IDs for chosen events
event_ids = "\n".join(sorted(event_ids))
with open("declustered_events.txt", "w") as f:
    f.write(event_ids)
