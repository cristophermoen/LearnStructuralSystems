
#containers for holding information about our analysis

struct member

    node_i::Int64
    node_j::Int64
    cross_section::String
    material::String

end

struct node

    x::Float64
    y::Float64

end


struct material

    name::String
    E::Float64

end

struct cross_section

    name::String
    A::Float64

end

#define those containers
#user defined...

materials = [material("steel", 29000.0), material("aluminium", 15000.0), material("concrete", 4000)]

nodes = [node(0.0, 0.0), node(9.0*12, 11.0*12), node(18.0*12, 0.0)]

cross_sections = [cross_section("I-beam", 20.0), cross_section("angle", 10.0)]

members = [member(1, 2, "I-beam", "steel"), member(2, 3, "angle", "steel")]

#define the local element stiffness matrix

function construct_local_element_k_matrix(E, A, L)

    k = E*A/L [1.0 0.0 -1.0 0.0
               0.0 0.0  0.0 0.0
               -1.0 0.0 1.0 0.0
               0.0 0.0 0.0 0.0]

    return k

end

#calculate member length

#define k element stiffness matrix in local coordinates

#loop for all elements

num_members = size(members)[1]

for i = 1:num_members
    
    k_element_local[i] = construct_local_element_k_matrix(materials[1].E, cross_sections[1].A, L)

end

#define the transformation matrix

k_element_global


#transform k element matrix from local to global (real world)

#assemble global stiffness matrix

#partition global stiffness matrix

#solve for displacements at free degrees of freedom




