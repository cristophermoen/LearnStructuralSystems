
# #x y
# node_geometry = [0.0 0.0
#                  1.0 3.0
                 
#                  ....
#                  ]

#              #element# #node_start  #node_end #cross_section_type #material_type
# element_info = [1 1 2 1 1  #concrete chord 
#                 2 3 4 2 2
                
#                 ..........
                
#                 ] #steel diagonal

num_nodes = 51
node_geometry = zeros(Float64, (num_nodes, 2))

node_geometry[:, 1] .= 0.0:30/(num_nodes-1):30.0  #x
node_geometry[:, 2] .= node_geometry[:, 1] .* 0.5 #y


num_elements = 50
element_info=zeros(Int, (num_elements, 5))

#deal with the top chord
num_top_chord_elements = 15

#fill element number
element_info[1:15, 1] .= 1:num_top_chord_elements

#fill node_start number
element_info[1:15, 2] .= 1:num_top_chord_elements

#fill node_end number
element_info[1:15, 3] .= 2:num_top_chord_elements+1

#fill cross-section type=1 for concrete top chord 
element_info[1:15, 4] .= 1

#fill material type = 1 for concrete
element_info[1:15, 5] .= 1


#deal with the bottom chord



#deal with the steel diagonals 


#calculating member length


using LinearAlgebra

element_length = zeros(Float64, num_elements)
element_angle = zeros(Float64, num_elements)

for i=1:num_top_chord_elements

    node_start_number = element_info[i, 2]
    node_end_number = element_info[i, 3]

    node_start_coords = node_geometry[node_start_number, :] 
    node_end_coords = node_geometry[node_end_number, :]

    element_vector = node_end_coords - node_start_coords

    element_length[i] = norm(element_vector)

    element_angle[i] = atan(element_vector[1], element_vector[2]) #radians

end