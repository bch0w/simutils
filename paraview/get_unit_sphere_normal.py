"""
Written by ChatGPT-4 for ParaView
"""
import math

def spherical_to_cartesian(lat, lon):
    # Convert latitude and longitude to radians
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    # Calculate x, y, z vector directions on unit sphere
    x = math.cos(lat_rad) * math.cos(lon_rad)
    y = math.cos(lat_rad) * math.sin(lon_rad)
    z = math.sin(lat_rad)

    return (x, y, z)

# Put these directly into the slice to get unit normal
x, y, z = spherical_to_cartesian(lat, lon, depth)

# Calculate normal vector
normal_vector = (x, y, z)

# Calculate radius at given depth
radius = 1 - depth / 6371
