#EN.560.301 Structural Systems I 
#Lecture 6
#February 2023

using LinearAlgebra

#Calculate the deformation of a 2D truss element.

#Define section properties.
A = [1.0]

#Define material properties.
E = [29000.0]

#Define nodal coordinates.
node_coordinates = [[0.0, 0.0], [48.0, 0.0]]

#Define element connectivity.
element_connectivity = [(1, 2), (2,3), (3,4)]

#Define support degrees of freedom.
supports = [1, 2, 4]

#Define external loads.
num_dof = size(node_coordinates, 1) * 2
F = zeros(Float64, num_dof)
F[3] = 10.0 

####Calculations#####


#element length
L = norm(node_coordinates[2] - node_coordinates[1])

#element orientation
Δy = node_coordinates[2][2] - node_coordinates[1][2]
Δx = node_coordinates[2][1] - node_coordinates[1][1]
θ = atan(Δy, Δx)

#local element stiffness matrix
function define_k_local(E, A, L)
    
    k_local = zeros(Float64, 4, 4)
    k_local[1, 1] = E*A/L
    k_local[1, 3] = -E*A/L
    k_local[3, 1] = -E*A/L
    k_local[3, 3] = E*A/L

    return k_local

end

k_local = define_k_local(E[1], A[1], L)
                
#rotation matrix 
function define_rotation_matrix(θ)

    β = zeros(Float64, 4, 4)

    unit_rotation = [cos(θ) sin(θ)
                     -sin(θ) cos(θ)]

    β[1:2, 1:2] = unit_rotation
    β[3:4, 3:4] = unit_rotation

    return β

end

β = define_rotation_matrix(θ)

#define global element K matrix 
k_global = β' * k_local * β


#assemble 
K = k_global

#partition
p = [3]
s = [1, 2, 4]

Kpp = K[p, p]
Kss = K[s, s]
Kps = K[p, s]
Ksp = K[s, p]

Fp = F[p]

#define imposed displacements
u = zeros(Float64, num_dof)
us = u[s] 

#solve for displacements at free degrees of freedom
up = Kpp \ (Fp - Kps * us)