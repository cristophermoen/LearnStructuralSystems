
using LinearAlgebra
include("/Users/crismoen/Documents/dev/LearnStructuralSystems/Spring_2025/Lecture_8_truss_system_stiffness/truss_helpers.jl")
#Analyze a truss, calculate displacement field when some loads are applied 


nodes = Vector{Vector{Float64}}(undef, 0)

#[x, y]   #inches
push!(nodes, [0.0, 0.0])
push!(nodes, [3, 4])
push!(nodes, [10.0, -10.0])


#material properties    #inches lbf 
# library 
#                    E  steel        wood    carbon fiber
material_properties = [29500000.0, 10000.0, 15000000.0]

#define section properties
#library 

         #cross-sectional area A   inches^2            
section_properties = [0.1, 10.0, 3.6]

#[node i  to node j  material property, section properties] 
elements = [[1, 2, 1, 1], [1, 3, 3, 1], [2, 3, 2, 3]]


###############

A = [section_properties[elements[i][4]] for i in eachindex(elements)]
E = [material_properties[elements[i][3]] for i in eachindex(elements)]
L = [calculate_element_length(nodes[elements[i][1]], nodes[elements[i][2]]) for i in eachindex(elements)]

θ = [calculate_element_orientation(nodes[elements[i][1]], nodes[elements[i][2]]) for i in eachindex(elements)]

ke_local = [calculate_element_local_stiffness_matrix(E[i], A[i], L[i]) for i in eachindex(L)]

β = [define_element_rotation_matrix(θ[i]) for i in eachindex(θ)]

ke_global = [define_global_element_stiffness_matrix(ke_local[i], β[i]) for i in eachindex(ke_local)]

element_global_dof = [define_global_dof_for_element(elements[i][1], elements[i][2]) for i in eachindex(elements)]

num_dof = size(nodes)[1] * 2

#Big Mac
K = zeros(Float64, num_dof, num_dof)
K = assemble_global_stiffness_matrix(num_dof, ke_global, element_global_dof)


# for i in eachindex(ke_global)

#     i = 2
#     K[element_global_dof[i], element_global_dof[i]] += ke_global[i]
    
# end


