
using LinearAlgebra
include("/Users/crismoen/Documents/dev/LearnStructuralSystems/Spring_2025/Lecture_8_truss_system_stiffness/truss_helpers.jl")
#Analyze a truss, calculate displacement field when some loads are applied 


nodes = Vector{Vector{Float64}}(undef, 0)

#[x, y]   #inches
push!(nodes, [0.0, 0.0])
push!(nodes, [-100.0, -100.0])
push!(nodes, [100.0, -100.0])


#material properties    #inches lbf 
# library 
#                    E  steel        wood    carbon fiber
material_properties = [29000000.0, 10000.0, 15000000.0]

#define section properties
#library 

         #cross-sectional area A   inches^2            
section_properties = [1.0, 2.0, 3.6]

#[node i  to node j  material property, section properties] 
elements = [[1, 2, 1, 2], [1, 3, 3, 1], [2, 3, 2, 3]]

#supports 
#Define fixed degrees of freedom
s = [3, 4, 5, 6]

#Define external forces
F = zeros(Float64, size(nodes)[1]*2)
F[1] = 100000.0  #lbs 

#Define any imposed displacements 
u = zeros(Float64, size(nodes)[1]*2)


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


#partition Big Mac 

p = setdiff(1:num_dof, s)  #free dof

Kpp = K[p, p]
Kss = K[s, s]
Kps = K[p, s]
Ksp = K[s, p]

Fp = F[p]


####SOLVE FOR GLOBAL DISPLACEMENTS
us = u[s] #if there are any imposed displacements, get them here
up = Kpp \ (Fp - Kps * us)
# up = Kpp^-1 * (Fp - Kps * us)

#fill displacement vector with solution
u[p] = up

#SOLVE FOR REACTIONS
Fs = Ksp * up + Kss * us

####SOLVE FOR ELEMENT INTERNAL FORCES

#find global displacements at node i and node j of each element
u_element_global = [u[element_global_dof[i]] for i in eachindex(element_global_dof)]

#find element displacements in element local coordinate system 
δ = [β[i] * u_element_global[i] for i in eachindex(u_element_global)]

#solve for internal forces in each element 
P = [ke_local[i] * δ[i] for i in eachindex(ke_local)] 




