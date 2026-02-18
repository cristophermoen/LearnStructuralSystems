
using CrossSectionGeometry, SectionProperties, AISIS100, CairoMakie

#define cross section 
t = 5/64
L = [1 + 11/64, 1 + 23/32, 1 + 11/64]
θ = [-π/2, 0.0, π/2]
n = [3, 3, 3]
r = [3t, 3t]
n_r = [5, 5]

section = CrossSectionGeometry.create_thin_walled_cross_section_geometry(L, θ, n, r, n_r, t, centerline="to left", offset=(0.0, 0.0))

#get centerline coordinates 
X = [section.centerline_node_XY[i][1] for i in eachindex(section.centerline_node_XY)]
Y = [section.centerline_node_XY[i][2] for i in eachindex(section.centerline_node_XY)]

#zero out cross section on X-Y axis 
X .-= minimum(X)
Y .-= minimum(Y)

X .+= t/2
Y .+= t/2


#show cross section 

aspect_ratio = (1 + 23/32)/(1 + 11/64)
f = Figure(size = (800, 800 / aspect_ratio))
ax = Axis(f[1, 1])
scatterlines!(ax, X, Y)
f


#calculate section properties 
coord = [X Y]
ends = [1:length(X)-1 2:length(X) fill(t, length(X)-1)]
properties = SectionProperties.cutwp_prop2(coord, ends)


#calculate critical elastic buckling loads
E = 29000000.0 #psi 
L = 51.0 / 2 # inches 
Ixx = properties.Ixx
Iyy = properties.Iyy
A = properties.A

Fcre(E, I, L, A) = π^2 * E * I / (L^2 * A)

Fcre_diagonal = minimum([Fcre(E, Ixx, L, A), Fcre(E, Iyy, L, A)])

#calculate web diagonal compressive strength 
Fy = 50000.0 #psi 
strength = AISIS100.v2024.e2(Fcre_diagonal, Fy, A, "AISI S100-2024 LRFD")
Pne = strength.Rn