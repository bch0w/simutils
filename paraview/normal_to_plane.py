"""
To cut an arbitrary plane through a rectangular volume, Paraview requires an
origin point, and the plane normal. Simple function returns the values required
to cut the plane which consists of a dot product
"""
import numpy as np

origin = input("Origin point (x,y,z): ")
point2 = input("Second point (x,y,z): ")
point3 = input("Third point (if None, set to origin point with z=-100): " )

origin = [float(_) for _ in origin.split(",")]
point2 = [float(_) for _ in point2.split(",")]
if point3:
    point3 = [float(_) for _ in point3.split(",")]
else:
    point3 = [origin[0], origin[1], -100]

qr = np.subtract(point2, origin)
qs = np.subtract(point3, origin)
normal = np.cross(qr, qs)

# Normal values dont need magnitude
oom = np.floor(np.log10(np.abs(normal[0]))) + 1
normal = [_/10**oom for _ in normal]

print(f"\nOrigin: {origin}")
print(f"Normal: {normal}\n")
