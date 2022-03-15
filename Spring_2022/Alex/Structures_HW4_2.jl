
## Structs

struct node

    x::Float64
    y::Float64

end

mutable struct member

    nodei::node
    nodej::node
    E::Float64
    Area::Float64

end

global E_steel = 200 #Gpa
global downward_force = 0.0 #giganewtons

## Functions
function compile_points_to_nodes(text)
    i = 1
    Doc = open(text)
    points = readlines(Doc)
    close(Doc)
    nodes = Array{node, 1}(undef, parse(Int64, points[1]))
    while i <= parse(Int64, points[1])
        current_point = points[i+1]
        coords = split(current_point; limit=2)
        x_coord = parse(Float64, coords[1])
        y_coord = parse(Float64, coords[2])
        nodes[i] = node(x_coord, y_coord)
        i += 1
    end
    return nodes
end

function construct_elements_hw4(node_set)
    i = 1
    elements = Array{member, 1}(undef, size(node_set)[1] - 1)

    while i <= size(node_set)[1] - 1
        elements[i] = member(node_set[i], node_set[i+1], E_steel, (0.1*0.1))
        i += 1
    end

    lin_elements = Array{member, 1}(undef, size(node_set)[1] - 2)

    for j=1:(size(node_set)[1]-2)
        if iseven(j)
            lin_elements[j] = member(node_set[j], node_set[j+2], E_steel, (0.47 * 0.3))
        else
            lin_elements[j] = member(node_set[j], node_set[j+2], E_steel, (0.35 * 0.35))
        end
    end
    all_elements = vcat(elements, lin_elements)
    return all_elements
end

function find_member_length(member)

        x1 = member.nodei.x
        x2 = member.nodej.x
        y1 = member.nodei.y
        y2 = member.nodej.y

        member_length = sqrt((x1 - x2)^2 + (y1 - y2)^2)

    return member_length
end

function define_truss_element_orientations(members)
    theta = zeros(Float64, length(members))

    for i=1:length(members)

        x1 = members[i].nodei.x
        y1 = members[i].nodei.y
        x2 = members[i].nodej.x
        y2 = members[i].nodej.y

        element_vector = [(x2 - x1)
                          (y2 - y1)]

        theta[i] = atan(element_vector[2], element_vector[1])
    end


    return theta
end

function construct_member_matrix(member, member_set, theta)

    theta_index = findall(x -> x == member, member_set)[1]
    angle = theta[theta_index]
    L = find_member_length(member)
    E = member.E
    A = member.Area

    k = (E*A/L)* [1.0 0.0 -1.0 0.0
                 0.0 0.0 0.0 0.0
               -1.0 0.0 1.0 0.0
               0.0 0.0 0.0 0.0]

    T_matrix =           [cos(angle) -sin(angle) 0 0
                          sin(angle) cos(angle) 0 0
                          0 0 cos(angle) -sin(angle)
                          0 0 sin(angle) cos(angle)]



        full_matrix =  (T_matrix) * k * inv(T_matrix)

        return full_matrix
end

function assemble_global_stiffness_matrix(node_geometry, members, k_element_global, num_dof_per_node)
    num_nodes = size(node_geometry)[1]

    k_system_global = zeros(Float64, num_nodes * num_dof_per_node, num_nodes * num_dof_per_node)

    for i=1:length(members)
        node_i = members[i].nodei
        node_j = members[i].nodej

        node_i_index = findall(x -> x == node_i, node_geometry)[1]
        node_j_index = findall(x -> x == node_j, node_geometry)[1]

        node_i_dof = Array{Int64, 1}(undef, num_dof_per_node)
        node_j_dof = Array{Int64, 1}(undef, num_dof_per_node)

        ##will need to alter this dof assignment part for 3dof per node, not too hard

        node_i_dof[2] = node_i_index * 2
        node_i_dof[1] = (node_i_index * 2) - 1

        node_j_dof[2] = node_j_index * 2
        node_j_dof[1] = (node_j_index * 2) - 1

        ##

        global_dof = [node_i_dof; node_j_dof]

        k_system_global[global_dof, global_dof] += k_element_global[i]
    end
    return k_system_global

end
## Commands

# cd("C:\\Users\\Alex Echtler\\.atom\\Julia Projects") #= <- This is just to move the directory to
# the place where the text file is at. I didnt feel like putting all of the nodes in by hand so a
# text file of coordinates and a node compiler seemed easier and more easy to use in the future.
# Youre going to have to either move the text file into the directory already and delete this line,
# or alter the path to reach the file, or just change this to where you downloaded the file.
# =#

text = "/Users/crismoen/Documents/dev/LearnStructuralSystems/Spring_2022/Alex/hw4 points.txt"

node_set = compile_points_to_nodes(text)

x = [node_set[i].x for i in eachindex(node_set)]
y = [node_set[i].y for i in eachindex(node_set)]


using Plots
plot(x, y, markershape = :o, seriestype = :scatter)


element_set = construct_elements_hw4(node_set)

theta = define_truss_element_orientations(element_set)

k_global_matrix_array = Array{Array{Float64, 2}, 1}(undef, size(element_set)[1])

for i=1:size(element_set)[1]
    k_global_matrix_array[i] = construct_member_matrix(element_set[i], element_set, theta)
end

K_global_matrix = assemble_global_stiffness_matrix(node_set, element_set, k_global_matrix_array, 2)

length_set = Array{Float64, 1}(zeros(size(element_set)[1]))

for i=1:size(element_set)[1]
    length_set[i] = find_member_length(element_set[i])
end

volume_set = Array{Float64, 1}(zeros(size(element_set)[1]))

for i=1:size(element_set)[1]
volume_set[i] = length_set[i] * element_set[i].Area
end

for i=1:size(element_set)[1]
    global downward_force += - 9.8 * 7850/1000000000 * volume_set[i] #giganewtons
end

forces = Array{Float64, 1}(zeros((2*size(node_set)[1])))

for i=1:2*size(node_set)[1]
    if iseven(i)
        forces[i] = downward_force/(size(node_set)[1])
    end
end

displacement = inv(K_global_matrix) * forces #m^2 ?

#= not sure if this unit conversion is supposed to happen, but either way the numbers are way to big

for i=1:size(displacement)[1] #m^2 -> m
    if displacement[i] < 0
        displacement[i] = -1 * displacement[i]
        displacement[i] = sqrt(displacement[i])
        displacement[i] = -1 * displacement[i]
    else
        displacement[i] = sqrt(displacement[i])
    end
end

=#
