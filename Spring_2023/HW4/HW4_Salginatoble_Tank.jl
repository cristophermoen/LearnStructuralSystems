#EN.560.301 Structural Systems I 
#HW 4 solution
#March 2023

using LinearAlgebra

#Calculate displacements and internal forces for the Salginatoble bridge.

#Define section properties.
A = [0.32, 0.877, 0.877, 0.877, 0.877, 0.877, 0.877]  #mm^2

I = [1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]  #mm^4

#Define material properties.
E = [30.0] #kN/mm^2 

#Define nodal coordinates.
node_coordinates = [[0.0, -12000.0], [6000.0, -9000.0], [12000.0, -6000.0], [18000.0, -4000.0], [24000.0, -2500.0], [30000.0, -1000.0], [36000.0, -750.0], [40000.0, 0.0]]

#Define element connectivity.
element_connectivity = [(1, 2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8)]

#Define element material properties
element_mat_props = [1, 1, 1, 1, 1, 1, 1]

#Define element section properties
element_sect_props = [1, 2, 3, 4, 5, 6, 7]

#Define support degrees of freedom.
supports = [1, 2, 22]

#Define external loads.
num_dof = size(node_coordinates, 1) * 3  #3 dof per node   
F = zeros(Float64, num_dof)
F[2] = -1.875 #kN
F[5] = -1.875 
F[8] = -1.875 
F[11] = -1.875 
F[14] = -1.875 
F[17] = -1.875 
F[20] = -1.875 
F[23] = -1.875 

#Define imposed support displacements.
us = [0.0, 0.0, 0.0]


#We could also define imposed temperature changes here if we wanted to, at some point...


####Calculations#####


#element length
function calculate_element_length(node_i, node_j)
    L = norm(node_j - node_i)
end

L = [calculate_element_length(node_coordinates[element_connectivity[i][1]], node_coordinates[element_connectivity[i][2]]) for i in eachindex(element_connectivity)] 



#element orientations
function calculate_element_orientation(node_i, node_j)

    Δy = node_j[2] - node_i[2]
    Δx = node_j[1] - node_i[1]
    θ = atan(Δy, Δx)

    return θ

end

θ = [calculate_element_orientation(node_coordinates[element_connectivity[i][1]], node_coordinates[element_connectivity[i][2]]) for i in eachindex(element_connectivity)]



#local element stiffness matrix
function define_k_local(E, I, A, L)
    
    k_local = zeros(Float64, 6, 6)

    k_local[1, 1] = E*A/L
    k_local[1, 4] = -E*A/L
    k_local[2, 2] = 12*E*I/L^3
    k_local[2, 3] = 6*E*I/L^2
    k_local[2, 5] = -12*E*I/L^3
    k_local[2, 6] = 6*E*I/L^2
    k_local[3, 3] = 4*E*I/L
    k_local[3, 5] =  -6*E*I/L^2
    k_local[3, 6] = 2*E*I/L
    k_local[4, 4] = E*A/L
    k_local[5, 5] = 12*E*I/L^3
    k_local[5, 6] = -6*E*I/L^2
    k_local[6, 6] = 4*E*I/L

    #reflect symmetric values that aren't on the diagonal
    i = [1, 2, 2, 2, 3, 3, 5]
    j = [4, 3, 5, 6, 5, 6, 6]
    
    for z in eachindex(i)
        k_local[j[z], i[z]] = k_local[i[z], j[z]]
    end


    return k_local

end

k_local = [define_k_local(E[element_mat_props[i]], I[element_sect_props[i]], A[element_sect_props[i]], L[i]) for i in eachindex(element_connectivity)]


#rotation matrix 
function define_rotation_matrix(θ)

    β = zeros(Float64, 6, 6)

    unit_rotation = [cos(θ) sin(θ)
                     -sin(θ) cos(θ)]

    β[3, 3] = 1.0
    β[6, 6] = 1.0

    β[1:2, 1:2] = unit_rotation
    β[4:5, 4:5] = unit_rotation

    return β

end

β = [define_rotation_matrix(θ[i]) for i in eachindex(θ)]


function transform_local_to_global_element_stiffness(β, k_local)

    k_global = β' * k_local * β

end


#calculate element stiffness matrices in global coordinates
k_global = [transform_local_to_global_element_stiffness(β[i], k_local[i]) for i in eachindex(element_connectivity)]


#define global dof associated with each element 
function define_global_dof_for_element(node_number_i, node_number_j)

    
    dof_i = [node_number_i * 3 - 2, node_number_i * 3 - 1, node_number_i * 3]
    dof_j = [node_number_j * 3 - 2, node_number_j * 3 - 1, node_number_j * 3]

    dof = [dof_i; dof_j]   # [1, 2, 3],  [4, 5, 6]  => [1, 2, 3, 4, 5, 6]   

    return dof

end

element_global_dof = [define_global_dof_for_element(element_connectivity[i][1], element_connectivity[i][2]) for i in eachindex(element_connectivity)]



#assemble Big Mac K, i.e. the global K matrix
function assemble_global_stiffness_matrix(k_global, num_dof)

    K = zeros(Float64, num_dof, num_dof)

    for i in eachindex(k_global)

        K[element_global_dof[i], element_global_dof[i]] += k_global[i]
        
    end

    return K

end

K = assemble_global_stiffness_matrix(k_global, num_dof)


#partition global K matrix 
p = setdiff(1:num_dof, supports)  #generalize 
s = supports

Kpp = K[p, p]
Kss = K[s, s]
Kps = K[p, s]
Ksp = K[s, p]

#partition global external force vector
Fp = F[p]

#solve for displacements at free degrees of freedom
up = Kpp \ (Fp - Kps * us)

#fill out complete displacement vector
u = zeros(Float64, num_dof)
u[s] = us
u[p] = up

#solve for reactions 
Fs = Ksp * up + Kss * us

#find global displacements at node i and node j of each element
u_element_global = [u[element_global_dof[i]] for i in eachindex(element_global_dof)]

#find element displacements in element local coordinate system 
δ = [β[i] * u_element_global[i] for i in eachindex(u_element_global)]

#solve for internal forces in each element 
P = [k_local[i] * δ[i] for i in eachindex(k_local)] 