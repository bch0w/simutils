# Getting Started with ObsPy
Notes for Il-Sang and Jay started 7/14

## Links
- [Conda: Package manager for Python](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)
- [ObsPy: Python library for seismology](https://docs.obspy.org/)

## Install Instructions (run this once)
1. Open your `Terminal` application
2. Install [**Conda**](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)
  - *Follow the instructions for your given operating system*
  - *Note: in data science traditionally use Anaconda/Miniconda/Conda, a package manager
    used to install common data science libraries. Miniconda is the lightweight version*
3. Check that you have successfully installed **Conda** by running the following command 
   in your `Terminal`
```bash
conda
```
  You should see a help message displayed. If not, something has gone wrong.
4. In your `Terminal` copy-paste the following commands:
```bash
conda create -n obspy
conda activate obspy
conda install obspy ipython
```
5. You can test that the installation was successful by running Python and trying
  to import ObsPy. ObsPy is a Python-based toolkit for seismology that has a lot of
  built-in functionality for analyzing and processing seismic data. IPython is
  a Python interpreter app that has some nicer features compared to the standard
  Python interpreter
```bash
ipython  
```
This will start Python in your terminal, in the Python instance type
```python
import obspy
obspy.__version__
exit
```
If the above works then you have successfully installed ObsPy. If it does not 
then something has gone wrong

## Instructions

Now that you have Python, you can run the following codes that will let you
process the seismic data


### To Convert Data Counts to Velocity

The following script will take a MiniSEED file (standard seismic format) and 
convert the amplitudes from counts to velocity in meters per second.

1. Ask Joey/Sylvia what the sensor parameters were (sampling rate, filter phase,
  DC filter, preamp Db)
2. Download the corresponding `StationXML` file by saving it as a text file on
  your computer
3. Modify the script below where it says MODIFY
4. Start `IPython` in your `Terminal` and copy-paste the modified code
5. Converted files will be stored with the same file name but with a suffix 
  to indicate the new units (e.g., <FID>\_VEL.mseed)


```python
import sys
import os
import obspy
from obspy import read, read_inventory, UTCDateTime, Inventory
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.core.inventory.channel import Channel
from obspy.clients.nrl import NRL
from obspy.io.xseed.core import _read_resp


# vvv MODIFY PARAMETERS HERE vvv
data_fid = "590004471.0016.2025.04.08.00.00.00.000.Z.miniseed"   # path to the miniSeed data file
response_fid = "RESP.txt"  # path to RESP file downloaded from IRIS
output = "VEL"  # VEL: velocity, DISP: displacement, ACC: acceleration
# ^^^ MODIFY PARAMETERS HERE ^^^

# Code below will read and convert the SmartSolo data
st = read(data_fid)

# Create a corresponding response and inventory objects
resp = _read_resp(response_fid)[0][0][0].response


# Dummy values to be used for Inventory, these are not actually important 
lat = 0.
lon = 0.
elevation = 0.
depth = 0.
start_date = UTCDateTime("1990-01-01")
end_date = UTCDateTime()

channels = []
for tr in st:
    channel = Channel(code=tr.stats.channel,
                      location_code=tr.stats.location, latitude=lat,
                      longitude=lon, elevation=elevation, depth=depth,
                      response=resp)
    channels.append(channel)

station = Station(code=tr.stats.station, latitude=lat, longitude=lon,
                  elevation=elevation, channels=channels,
                  start_date=start_date, end_date=end_date
                  )

network = Network(code=tr.stats.network, stations=[station])

inv = Inventory(networks=[network])

# Remove the instrument response
st.remove_response(inv, output=output)

# Create a new filename
fid_new, ext = os.path.splitext(os.path.basename(data_fid))
fid_new = f"{fid_new}_{output}.{ext}"

st.write(fid_new, format="MSEED")
print(f"{os.path.basename(data_fid)} -> {os.path.basename(fid_new)} "
      f"with units {output}")
```

