#SS-Homework 3
import numpy as np
import matplotlib.pyplot as plt
from sectionproperties.pre.geometry import Geometry
from sectionproperties.analysis.section import Section

# Inputs (inches)
t = 0.0747
r_i = 0.125
w_i = 1.71875
H_out = 1.171875
mesh_size = 0.06

def arc(cx, cy, r, a0, a1, n=25):
    ang = np.linspace(a0, a1, n)
    pts = [(cx + r*np.cos(a), cy + r*np.sin(a)) for a in ang]
    return pts[1:-1]

# Geometry values
r_o = r_i + t
xL_o, yB_o = 0.0, 0.0
xR_o = w_i + 2*r_i + 2*t
yT = H_out

xL_i = xL_o + t
xR_i = xR_o - t
yB_i = yB_o + t

yTan_o = yB_o + r_o
yTan_i = yB_i + r_i

# Build shape
pts = []
pts += [(xL_o, yT), (xL_o, yTan_o)]
pts += arc(xL_o + r_o, yB_o + r_o, r_o, np.pi, 1.5*np.pi)
pts += [(xL_o + r_o, yB_o), (xR_o - r_o, yB_o)]
pts += arc(xR_o - r_o, yB_o + r_o, r_o, 1.5*np.pi, 2*np.pi)
pts += [(xR_o, yTan_o), (xR_o, yT)]

pts += [(xR_i, yT), (xR_i, yTan_i)]
pts += arc(xR_i - r_i, yB_i + r_i, r_i, 0.0, -0.5*np.pi)
pts += [(xR_i - r_i, yB_i), (xL_i + r_i, yB_i)]
pts += arc(xL_i + r_i, yB_i + r_i, r_i, -0.5*np.pi, -np.pi)
pts += [(xL_i, yTan_i), (xL_i, yT)]
pts += [(xL_o, yT)]  # close

facets = [[i, i + 1] for i in range(len(pts) - 1)]
control_points = [(xL_o + 0.5*t, 0.5*(yB_o + yT))]

# Create and analyze
geom = Geometry.from_points(points=pts, facets=facets, control_points=control_points)
geom.create_mesh(mesh_sizes=[mesh_size])

sec = Section(geom)
sec.calculate_geometric_properties()
sec.calculate_warping_properties()

# Plot 
fig, ax = plt.subplots()
geom.plot_geometry(ax=ax)
ax.set_aspect("equal", "box")
ax.set_xlabel("x [in]")
ax.set_ylabel("y [in]")
plt.show()

# Results
area = sec.get_area()
cx, cy = sec.get_c()
J = sec.get_j()

# Ixx, Iyy 
try:
    Ixx_c, Iyy_c, _ = sec.get_ic()
except Exception:
    Ixx_c = sec.get_ixx_c()
    Iyy_c = sec.get_iyy_c()

print(f"Area = {area:.6f} in^2")
print(f"Centroid (cx, cy) = ({cx:.6f}, {cy:.6f}) in")
print(f"Ixx (centroid) = {Ixx_c:.6f} in^4")
print(f"Iyy (centroid) = {Iyy_c:.6f} in^4")
print(f"J = {J:.6f} in^4")
print(f"Outside width = {xR_o - xL_o:.6f} in")
print(f"Total height  = {yT - yB_o:.6f} in")
