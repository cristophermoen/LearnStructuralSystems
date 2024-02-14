#EN.560.301 Structural Systems I 
#Lecture 8
#February 2024

using LinearAlgebra

#Sennett Example 2.1

###DEFINITIONS

#Define section properties.
A = [1.0, 1.0, 1.0]  #in^2

#Define material properties.
E = [1.0, 1.0, 1.0] #ksi

#Define nodal coordinates.
node_coordinates = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]  #in

#Define element connectivity.  i to j
element_connectivity = [(1, 2), (1, 3), (2, 3)]

#Define fixed degrees of freedom
s = [1, 2, 5, 6]

#Define external forces
F = zeros(Float64, 6)
F[4] = -1.0

#Define imposed displacements
u = zeros(Float64, 6)

####FUNCTIONS

function calculate_element_length(node_i_coordinates, node_j_coordinates)

    L = norm(node_j_coordinates - node_i_coordinates)

    return L

end

function calculate_element_orientation(node_i_coordinates, node_j_coordinates)

    Δ = node_j_coordinates - node_i_coordinates
    
    θ = atan(Δ[2], Δ[1])

    return θ

end

function define_local_element_stiffness_matrix(E, A, L)
    
    ke_local = zeros(Float64, 4, 4)
    ke_local[1, 1] = E*A/L
    ke_local[1, 3] = -E*A/L
    ke_local[3, 1] = -E*A/L
    ke_local[3, 3] = E*A/L

    return ke_local

end

function define_element_rotation_matrix(θ)

    β = zeros(Float64, 4, 4)

    unit_rotation = [cos(θ) sin(θ)
                     -sin(θ) cos(θ)]

    β[1:2, 1:2] = unit_rotation
    β[3:4, 3:4] = unit_rotation

    return β

end

function define_global_element_stiffness_matrix(ke_local, β)

    ke_global = β' * ke_local * β

    return ke_global

end

function define_global_dof_for_element(node_number_i, node_number_j)

    dof_i = [node_number_i * 2 - 1, node_number_i * 2]
    dof_j = [node_number_j * 2 - 1, node_number_j * 2]

    dof = [dof_i; dof_j]   # [1, 2],  [3, 4]  => [1, 2, 3, 4]   

    return dof

end

function assemble_global_stiffness_matrix(num_dof, ke_global_all_elements, element_global_dof_all_elements)

    K = zeros(Float64, num_dof, num_dof)

    for i in eachindex(ke_global_all_elements)

        K[element_global_dof_all_elements[i], element_global_dof_all_elements[i]] += ke_global_all_elements[i]
        
    end

    return K

end

####GLOBAL STIFFNESS MATRIX (BIG MAC)

L = [calculate_element_length(node_coordinates[element_connectivity[i][1]], node_coordinates[element_connectivity[i][2]]) for i in eachindex(element_connectivity)]

θ = [calculate_element_orientation(node_coordinates[element_connectivity[i][1]], node_coordinates[element_connectivity[i][2]]) for i in eachindex(element_connectivity)]

ke_local = [define_local_element_stiffness_matrix(E[i], A[i], L[i]) for i in eachindex(L)]

β = [define_element_rotation_matrix(θ[i]) for i in eachindex(θ)]

ke_global = [define_global_element_stiffness_matrix(ke_local[i], β[i]) for i in eachindex(ke_local)]

global_element_dof = [define_global_dof_for_element(element_connectivity[i][1], element_connectivity[i][2]) for i in eachindex(element_connectivity)]

num_dof = size(node_coordinates)[1] * 2

K = assemble_global_stiffness_matrix(num_dof, ke_global, global_element_dof)

####PARTITION K matrix and F vector
p = setdiff(1:num_dof, s)  #free dof

Kpp = K[p, p]
Kss = K[s, s]
Kps = K[p, s]
Ksp = K[s, p]

Fp = F[p]


####SOLVE FOR GLOBAL DISPLACEMENTS
us = u[s] #if there are any imposed displacements, get them here
up = Kpp \ (Fp - Kps * us)

#fill displacement vector with solution
u[p] = up

#SOLVE FOR REACTIONS
Fs = Ksp * up + Kss * us

####SOLVE FOR ELEMENT INTERNAL FORCES

#find global displacements at node i and node j of each element
u_element_global = [u[global_element_dof[i]] for i in eachindex(global_element_dof)]

#find element displacements in element local coordinate system 
δ = [β[i] * u_element_global[i] for i in eachindex(u_element_global)]

#solve for internal forces in each element 
P = [ke_local[i] * δ[i] for i in eachindex(ke_local)] 

