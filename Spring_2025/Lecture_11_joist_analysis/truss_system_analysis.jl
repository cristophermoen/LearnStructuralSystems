
using LinearAlgebra, CSV, DataFrames

include("./truss_helpers.jl")

#Analyze a truss, calculate displacement field when some loads are applied 

node_data = CSV.read("Spring_2025/Lecture_11_joist_analysis/nodes.csv", DataFrame)

element_data = CSV.read("Spring_2025/Lecture_11_joist_analysis/elements.csv", DataFrame)


#[X, Y]   #inches
nodes = Vector{Vector{Float64}}(undef, 0)
for i in eachindex(node_data.X)
    push!(nodes, [node_data.X[i], node_data.Y[i]])
end


#material properties    #inches lbf 
# library 
#                    E  steel     
material_properties = [29500000.0]

#define section properties
#library 

         #cross-sectional area A   inches^2    
         #chord, diagonal        
section_properties = [0.877, 0.30]

#[node i  to node j  material property, section properties] 
elements = Vector{Vector{Int}}(undef, 0)
for i in eachindex(element_data.node_i)
    push!(elements, [element_data.node_i[i], element_data.node_j[i], element_data.mat_prop[i], element_data.sect_prop[i]])
end


#Define fixed degrees of freedom
s = [1, 2, 24]

#Define any imposed displacements 
u = zeros(Float64, size(nodes)[1]*2)

#Define uplift uniform load on top chord 
L = [calculate_element_length(nodes[elements[i][1]], nodes[elements[i][2]]) for i in eachindex(elements)]
w = -300.0 / 12 #lbs/ft to lbs/in

F = zeros(Float64, size(nodes)[1]*2)

top_chord_elements = 1:11

for i in eachindex(top_chord_elements)

    F[2* i] += w * L[i]/2
    F[2* (i + 1)] += w * L[i]/2

end

###############

#get vector of cross-sectional areas 
A = [section_properties[elements[i][4]] for i in eachindex(elements)]
E = [material_properties[elements[i][3]] for i in eachindex(elements)]

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

#maximum uplift deflection
minimum(u[2:2:end])
argmin(u[2:2:end])

#SOLVE FOR REACTIONS
Fs = Ksp * up + Kss * us

####SOLVE FOR ELEMENT INTERNAL FORCES

#find global displacements at node i and node j of each element
u_element_global = [u[element_global_dof[i]] for i in eachindex(element_global_dof)]

#find element displacements in element local coordinate system 
δ = [β[i] * u_element_global[i] for i in eachindex(u_element_global)]

#solve for internal forces in each element 
P = [ke_local[i] * δ[i] for i in eachindex(ke_local)] 




