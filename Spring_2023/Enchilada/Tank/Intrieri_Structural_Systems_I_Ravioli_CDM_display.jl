#EN.560.301 Structural Systems I 
#Final Ravioli
#April 2023

using XLSX, InstantFrame, GLMakie 


xf = XLSX.readxlsx("/Users/crismoen/Documents/dev/LearnStructuralSystems/Spring_2023/Enchilada/Tank/Intrieri Structural Systems I Ravioli Nodes.xlsx")
sh = xf["Sheet1"]

#Calculate the deformation of a 2D joist truss cell.

#Define section properties.
A = [5.25] #in^2

I = [5.359375]  #in^4

h = [5.25/2] #in

#Define material properties.
E = [500000] #psi

#Define nodal coordinates.
node_coordinates = Vector{Vector{Float64}}(undef, 69)
for i in range(1,69)
    node_coordinates[i] = [sh[i+1, 2], sh[i+1, 3]]
end

#Define element connectivity.
element_connectivity = Vector{Tuple{Int64, Int64}}(undef, 137)
for i in range(1,137)
    element_connectivity[i] = (sh[i+1, 9], sh[i+1, 10])
end


#Define support degrees of freedom.
supports = Vector{Int64}(undef, 37)
for i in range(1, 37)
    supports[i] = sh[i+1, 14]
end

#Define Live Load [Case 1: Top of hill, Case 2: Middle of hill, Case 3: Bottom of Hill]
#Un-comment load conditions for preferred case ONLY
num_dof = size(node_coordinates, 1) * 3  #3 dof per node
F = zeros(Float64, num_dof)

#CASE 1
#node 13
# F[38] = -5291.094 #lbf


#CASE 2
#node 43
# F[127] = -4562.63176163
# F[128] = -4562.63176163


#CASE 3
#node 69
F[206] = -6452.53571741


#Define imposed support displacements.
us = Vector{Float64}(undef, 37)
for i in range(1, 37)
    us[i] = 0.0
end


#####InstantFrame mapping 


material = InstantFrame.Material(names=["wood"], E=[500000.0], ν=[0.40], ρ=[23.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["framing"], A=[5.25], Iy=[5.359375], Iz=[5.359375], J=[0.001])   #Iy and J not really important here

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

numbers = 1:size(node_coordinates, 1)
coordinates = [(node_coordinates[i][1], node_coordinates[i][2], 0.0) for i in eachindex(node_coordinates)]
node = InstantFrame.Node(numbers= numbers, coordinates= coordinates)

numbers = 1:size(element_connectivity, 1)
types = fill("frame", length(numbers))
nodes = element_connectivity
orientation = fill(0.0, length(numbers))
connections = fill(("rigid", "rigid"), length(numbers))
cross_section_type = fill("framing", length(numbers))
material_type = fill("wood", length(numbers))
element = InstantFrame.Element(types=types, numbers=numbers, nodes=nodes, orientation=orientation, connections=connections, cross_section=cross_section_type, material=material_type)


nodes = 2:13
uX = fill(Inf, length(nodes))
uY = fill(0.0, length(nodes))
uZ = fill(Inf, length(nodes))
rX = fill(0.0, length(nodes))
rY = fill(0.0, length(nodes))
rZ = fill(0.0, length(nodes))

nodes = [nodes; [1
14
26
36
44
50
55
58
61
64
66
68
14
26
36
44
50
55
58
61
64
66]]
uX = [uX; fill(Inf, 22)]
uY = [uY; fill(Inf, 22)]
uZ = [uZ; fill(Inf, 22)]
rX = [rX; fill(0.0, 22)]
rY = [rY; fill(0.0, 22)]
rZ = [rZ; fill(0.0, 22)]

nodes = [nodes; 68:69]
uX = [uX; fill(Inf, 2)]
uY = [uY; fill(0.0, 2)]
uZ = [uZ; fill(Inf, 2)]
rX = [rX; fill(0.0, 2)]
rY = [rY; fill(0.0, 2)]
rZ = [rZ; fill(0.0, 2)]


support = InstantFrame.Support(nodes=nodes, stiffness=(uX=uX, uY=uY, uZ=uZ, rX=rX, rY=rY, rZ=rZ))

uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["case 2"], nodes=[43], magnitudes=(FX=[-4562.63176163], FY=[-4562.63176163], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, "first order")


####show Wildcat 

element_nodal_coords = InstantFrame.Show.define_element_nodal_start_end_coordinates(element, node)

X, Y, Z = InstantFrame.Show.get_node_XYZ(node)

X_range = abs(maximum(X) - minimum(X))
Y_range = abs(maximum(Y) - minimum(Y))
Z_range = abs(maximum(Z) - minimum(Z))



figure = Figure()
ax = Axis3(figure[1,1])
# ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)
# ax.yticks = WilkinsonTicks(2)
# ylims!(ax, 0.0, 50.0)

color = :brown
InstantFrame.Show.elements!(ax, element_nodal_coords, color)
figure

# markersize = 10
# color = :blue
# InstantFrame.Show.nodes!(ax, X, Y, Z, markersize, color)
# figure


# unit_arrow_head_size = [1.0, 1.0, 1.0]
# arrow_head_scale = 3.0
# arrow_scale = 8.0
# arrowcolor = :orange 
# linecolor = :orange 
# linewidth = 1

# InstantFrame.Show.element_local_axes!(ax, element, node, model, unit_arrow_head_size, arrow_head_scale, arrow_scale, arrowcolor, linecolor, linewidth)
# figure



# textsize = 10
# color = :darkred
# InstantFrame.Show.element_numbers!(ax, element, node, textsize, color)
# figure

textsize = 10
color = :green
InstantFrame.Show.node_numbers!(ax, node, textsize, color)
figure


unit_arrow_head_size = [1.0, 1.0, 1.0]
arrow_head_scale = 1.0
arrow_scale = 1.0
arrowcolor = :green 
linecolor = :green 
linewidth = 1.0
InstantFrame.Show.point_loads!(ax, point_load, node, arrow_scale, arrow_head_scale, unit_arrow_head_size, arrowcolor, linecolor, linewidth)
figure

n = fill(5, length(element.numbers))
scale = (1.0, 1.0, 1.0)
linecolor = :blue
InstantFrame.Show.deformed_shape!(ax, model.solution.displacements, model.properties.global_dof, element, node, model.properties, model.solution.connections, n, scale, linecolor)