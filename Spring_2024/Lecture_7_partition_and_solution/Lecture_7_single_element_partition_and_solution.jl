#EN.560.301 Structural Systems I 
#Lecture 7
#February 2024

using LinearAlgebra

#Calculate the deformation of a 2D truss element.

#Define section properties.
A = [1.0, 2.3, 4.6]  #in^2

#Define material properties.
E = [29000.0, 15000.0, 1000.0] #ksi

#Define nodal coordinates.
node_coordinates = [[0.0, 0.0], [24.0, -48.0]]  #in

#Define element connectivity.  i to j
element_connectivity = [(1, 2)]

####Calculations#####


#element length
L = norm(node_coordinates[2] - node_coordinates[1])


node_coordinates = [[0.0, 0.0], [24.0, -48.0]]
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
k_e_global = β' * k_local * β

K = k_e_global

#partition
p = [4]  #free
s = [1, 2, 3]

Kpp = K[p, p]
Kss = K[s, s]
Kps = K[p, s]
Ksp = K[s, p]

num_dof = size(node_coordinates)[1] * 2
F = zeros(Float64, num_dof)
F[4] = -1000

Fp = F[p]

#define imposed displacements
u = zeros(Float64, num_dof)
us = u[s] 

#solve for displacements at free degrees of freedom
up = Kpp \ (Fp - Kps * us)