#EN.560.301 Structural Systems I 
#HW 5 solution
#March 2023

using LinearAlgebra

#Calculate the deformation of a 2D joist truss cell.

#Define section properties.
A = [0.5244, 0.7866, 1.8354, 3.1464, 2.622, 1.311, 0.6992] #m^2

I = [0.107828753334, 0.261918690834, 3.14535969083, 8.17548125334, 6.13346075334, 1.10333919084, 0.211154336667]  #m^4

h = [0.75/2, 1.125/2, 2.625/2, 4.5/2, 3.75/2, 1.875/2, 1/2]

#Define material properties.
E = [27579029.17] #kPa #(Converted from ~4000 ksi)

#Define nodal coordinates.
node_coordinates = [[0.0, -12.0], [6.0, -9.0], [12.0, -6.0], [18.0, -4.0], [24.0, -2.5], [30.0, -1.0], [36.0, -0.3], [42.0, 0.0]] #m

#Define element connectivity.
element_connectivity = [(1, 2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8)]

#Define element material properties
element_mat_props = [1, 1, 1, 1, 1, 1, 1]

#Define element section properties
element_sect_props = [1, 2, 3, 4, 5, 6, 7]

#Define support degrees of freedom.
supports = [1, 2, 22]

#Define external loads (Live + Dead)
num_dof = size(node_coordinates, 1) * 3  #3 dof per node
F = zeros(Float64, num_dof)
#calculate weight below, then finalize F vector.



#Define imposed support displacements.
us = [0.0, 0.0, 0.0]


#We could also define imposed temperature changes here if we wanted to, at some point...


####Calculations#####


#element length
function calculate_element_length(node_i, node_j)
    L = norm(node_j - node_i)
end

L = [calculate_element_length(node_coordinates[element_connectivity[i][1]], node_coordinates[element_connectivity[i][2]]) for i in eachindex(element_connectivity)] 

#element weight
function calculate_element_weight(A, L, ρ)
    W = A * L * ρ
end

W = [calculate_element_weight(A[i], L[i], 19.6) for i in eachindex(A)]

sum(W)


#element orientations
function calculate_element_orientation(node_i, node_j)

    Δy = node_j[2] - node_i[2]
    Δx = node_j[1] - node_i[1]
    θ = atan(Δy, Δx)

    return θ

end

θ = [calculate_element_orientation(node_coordinates[element_connectivity[i][1]], node_coordinates[element_connectivity[i][2]]) for i in eachindex(element_connectivity)]


#local element stiffness matrix
function define_k_local(E, A, I, L)
    
    k_local = zeros(Float64, 6, 6)
    k_local[1, 1] = E*A/L
    k_local[1, 4] = -E*A/L
    k_local[2, 2] = 12*E*I/L^3
    k_local[2, 3] = 6*E*I/L^2
    k_local[2, 5] = -12*E*I/L^3
    k_local[2, 6] = 6*E*I/L^2
    k_local[3, 2] = 6*E*I/L^2
    k_local[3, 3] = 4*E*I/L
    k_local[3, 5] = -6*E*I/L^2
    k_local[3, 6] = 2*E*I/L
    k_local[4, 1] = -E*A/L
    k_local[4, 4] = E*A/L
    k_local[5, 2] = -12*E*I/L^3
    k_local[5, 3] = -6*E*I/L^2
    k_local[5, 5] = 12*E*I/L^3
    k_local[5, 6] = -6*E*I/L^2
    k_local[6, 2] = 6*E*I/L^2
    k_local[6, 3] = 2*E*I/L
    k_local[6, 5] = -6*E*I/L^2
    k_local[6, 6] = 4*E*I/L
    return k_local

end

k_local = [define_k_local(E[element_mat_props[i]], A[element_sect_props[i]], I[element_sect_props[i]], L[i]) for i in eachindex(element_connectivity)]

#update this... Tank 
#rotation matrix 
function define_rotation_matrix(θ)

    β = zeros(Float64, 6, 6)

    unit_rotation = [cos(θ) sin(θ)
                     -sin(θ) cos(θ)]

    β[1:2, 1:2] = unit_rotation
    β[3, 3] = 1.0
    β[4:5, 4:5] = unit_rotation
    β[6, 6] = 1.0

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


#define dead load 

num_dof = size(node_coordinates, 1) * 3  #3 dof per node
F = zeros(Float64, num_dof)
F[2] = -W[1]/2
F[5] = -(W[1] + W[2])/2
F[8] = -(W[2] + W[3])/2
F[11] = -(W[3] + W[4])/2
F[14] = -(W[4] + W[5])/2
F[17] = -(W[5] + W[6])/2
F[20] = -(W[6] + W[7])/2
F[23] = -W[7]/2



#partition global K matrix 
p = setdiff(1:num_dof, supports)
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


using Plots 

x = [node_coordinates[i][1] for i in eachindex(node_coordinates)]
y = [node_coordinates[i][2] for i in eachindex(node_coordinates)]
plot(x, y, markershape = :o, aspect_ratio=:equal, legend=false)


nodal_displacements = [u[(1:3).+3*(i-1)] for i in eachindex(node_coordinates)]

ΔX = [nodal_displacements[i][1] for i in eachindex(nodal_displacements)]
ΔY = [nodal_displacements[i][2] for i in eachindex(nodal_displacements)]

scale=100.0
plot!(x.+ΔX*scale, y.+ΔY*scale, markershape = :o, markercolor=:red, legend=false)


axial_force = [-P[1][1]; [P[i][4] for i in eachindex(P)]]

moment = [0.0; [P[i][6] for i in eachindex(P)]]

shear = [-P[1][2]; [P[i][5] for i in eachindex(P)]]


A_node = [0.3496, 0.6992, 0.874, 2.7968, 3.496, 1.748, 0.874, 0.5244]

I_node = [0.027, 0.189, 0.335, 5.955, 10.395, 1.871, 0.335, 0.087]

depth_node = [0.25, 0.5, 1.25/2, 2.0, 2.5, 1.25, 1.25/2, 0.75/2] * 2

axial_stress = [fill(axial_force[i] ./ A_node[i], 3) for i in eachindex(axial_force)]

flexural_stress = [[-moment[i]*depth_node[i]/2/I_node[i], 0.0, moment[i]*depth_node[i]/2/I_node[i]] for i in eachindex(moment)]

shear_stress = [[0.0, 3/2*shear[i]/A_node[i], 0.0] for i in eachindex(moment)]

#Find Stress Tensor
function create_stress_tensor(σ_axial, σ_flexural, τ)
    T = zeros(Float64, 2, 2)

    T[1,1] = σ_axial + σ_flexural
    T[1,2] = τ
    T[2,1] = τ

    return T

end


T = Vector{Any}(undef, size(node_coordinates, 1))
principal_stresses = Vector{Any}(undef, size(node_coordinates, 1))
principal_stress_directions = Vector{Any}(undef, size(node_coordinates, 1))

for i in eachindex(T)

    T[i] = [create_stress_tensor(axial_stress[i][j], flexural_stress[i][j], shear_stress[i][j]) for j in eachindex(axial_stress[i])]

    principal_stresses[i] = [eigvals(T[i][j]) for j in eachindex(axial_stress[i])]

    principal_stress_directions[i] = [eigvecs(T[i][j]) for j in eachindex(axial_stress[i])]


end


kN/m^2

kPa

MPa   -4.3 #MPa

0.85 * 27