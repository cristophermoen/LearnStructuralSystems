
using LinearAlgebra
include("/Users/crismoen/Documents/dev/LearnStructuralSystems/Spring_2025/Lecture_7_stiffness/truss_helpers.jl")
#Analyze a truss, calculate displacement field when some loads are applied 


nodes = Vector{Vector{Float64}}(undef, 0)

#[x, y]
push!(nodes, [0.0, 0.0])
push!(nodes, [-10.0, -10.0])
push!(nodes, [10.0, -10.0])


#material property 
#                    E  steel        wood    carbon fiber
material_properties = [29000000.0, 10000.0, 15000000.0]

#define section properties

         #cross-sectional area              
section_properties = [1.0, 10.0, 3.6]

#[node i  to node j  material property, section properties] 
elements = [[1, 2, 1, 2], [1, 3, 3, 1], [2, 3, 2, 3]]

e = 1 #set element number!

E = material_properties[elements[e][3]]

node_i_coordinates = nodes[elements[e][1]]
node_j_coordinates = nodes[elements[e][2]]
L = calculate_element_length(node_i_coordinates, node_j_coordinates)

A = section_properties[elements[e][4]]

ke_local = calculate_element_local_stiffness_matrix(E, A, L)