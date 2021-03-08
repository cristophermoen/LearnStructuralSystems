#EN.560.301 Structural Systems I
#Spring 2021


#Perform a structural analysis of Eduardo Torroja's CASA aircraft manufacturing building.  Focus on a typical parabolic arch that support the doubly-curved thin concrete shell roof.   The arch has concrete top and bottom chords and a steel diagonal truss network.

#HELPER FUNCTIONS

function parabola(h, L, x)

    b = 4 .* h / L
    a = -b ./ L
    c = 0.0

    y = a .* x.^2 + b .* x .+ c

    dy_dx = 2 .* a .* x .+ b

    return y, dy_dx

end

#DEFINE ARCH GEOMETRY

#Define the arch bottom chord centerline span.
L_bottom_chord = 73000 #mm

#Define the arch bottom chord centerline height.
h_bottom_chord = 21000 + 350/2 - 4000 #mm

#Define the number of equal bottom chord arch segments.
num_bottom_chord_arch_segments = 36

#Define the centerline distance between top chord and bottom chord of arch.
chord_spacing = (22800 - 470/2) - (21000 - 350/2)

#Define the arch top chord centerline span.

#We need the arch slope at springing point of the arch.
x = [0.0]
y, dy_dx = parabola(h_bottom_chord, L_bottom_chord, x)
θ_springing_point = atan(dy_dx[1])
delta_x_bottom_to_top_chord = chord_spacing * cos(θ_springing_point)

L_top_chord = L_bottom_chord + 2 * delta_x_bottom_to_top_chord 

#Define the arch top chord centerline height.
h_top_chord = 22800 - 470/2 - 4000 #mm

#Define the number of equal top chord arch segments.
num_top_chord_arch_segments = 36



# DISCRETIZE BOTTOM ARCH CHORD

#Calculate the bottom arch centerline arc length.  


x = 0:L_bottom_chord/200:L_bottom_chord
y, dy_dx = parabola(h_bottom_chord, L_bottom_chord, x)
dx = diff(x)
dy = diff(y)
ds = sqrt.(dx.^2 + dy.^2)

s_bottom_chord = sum(ds)  #This is the length of the bottom chord.  

#Calculate the typical length of a bottom chord segment between truss diagonals.
length_bottom_chord_segment = s_bottom_chord / num_bottom_chord_arch_segments

#Define all the bottom chord segments along the bottom chord.
bottom_chord_arch_layout = [length_bottom_chord_segment/2 ; ones(Float64, num_bottom_chord_arch_segments - 1) .* length_bottom_chord_segment; length_bottom_chord_segment/2]

#Calculate the x-y locations of the nodes along the bottom chord that define the equal arch segments.

function calculate_arch_xy(arch_layout, h, L)
	
    x_index = 0.0

    num_segments = length(arch_layout)

    dx_segment = zeros(Float64, num_segments)
    
    num_nodes = num_segments + 1
    x_arch = zeros(Float64, num_nodes)
    y_arch = zeros(Float64, num_nodes)

	for i = 1:num_segments
	
		arch_segment = arch_layout[i]
		
        y, dy_dx = parabola(h, L, x_arch[i])
        
        θ = atan(dy_dx)
	
		dx_segment[i] = arch_segment * cos(θ)
		
        x_index = x_index + dx_segment[i]
        
        x_arch[i+1] = x_index
		
	end

    y_arch, dy_dx = parabola(h, L, x_arch)

	return x_arch, y_arch
	
end

x_bottom_chord, y_bottom_chord = calculate_arch_xy(bottom_chord_arch_layout, h_bottom_chord, L_bottom_chord)




# DISCRETIZE TOP ARCH CHORD

x = 0:L_top_chord/200:L_top_chord
y, dy_dx = parabola(h_top_chord, L_top_chord, x)
dx = diff(x)
dy = diff(y)
ds = sqrt.(dx.^2 + dy.^2)

s_top_chord = sum(ds)  #This is the length of the bottom chord.  

#Calculate the typical length of a bottom chord segment between truss diagonals.
length_top_chord_segment = s_top_chord / num_top_chord_arch_segments

#Define all the bottom chord segments along the bottom chord.
top_chord_arch_layout = ones(Float64, num_top_chord_arch_segments) * length_top_chord_segment

x_top_chord, y_top_chord = calculate_arch_xy(top_chord_arch_layout, h_top_chord, L_top_chord)

#Shift the top chord x-coordinates.
x_top_chord = x_top_chord .- delta_x_bottom_to_top_chord 

#Shift the top chord y-coordinates.
y_top_chord = y_top_chord .+ chord_spacing

# using Plots
# plot(x_bottom_chord, y_bottom_chord, markershape = :o, size=(600,600), legend = false, xlims = (-3000, L_bottom_chord), ylims = (0,L_bottom_chord))
# plot!(x_top_chord, y_top_chord, markershape = :o, size=(600,600), legend = false, xlims = (-3000, L_top_chord), ylims = (0,L_top_chord))


#CALCULATE SECTION PROPERTIES 

A_top_chord = 470.0 * 300.0   
A_bottom_chord = 350.0 * 350.0
A_diagonal = 100.0 * 100 - (100.0 - 5.0) * (100.0 - 5.0)

#DEFINE MATERIAL PROPERTIES

E_steel = 199947 #N/mm^2
E_concrete = 27579 #N/mm^2   #4000 ksi

#DEFINE DISCRETIZATION INFO
num_top_chord_nodes = length(x_top_chord)
num_bottom_chord_nodes = length(x_bottom_chord)
num_nodes = num_top_chord_nodes + num_bottom_chord_nodes

num_top_chord_elements = length(top_chord_arch_layout)
num_bottom_chord_elements = length(bottom_chord_arch_layout)
num_diagonal_elements = (num_bottom_chord_nodes - 2) * 2
num_elem = num_top_chord_elements + num_bottom_chord_elements + num_diagonal_elements


#DEFINE STRUCTURAL SYSTEM MODEL

#                 x  y
#Combine top chord and bottom chord nodes here into one matrix.
node_geometry = [[x_bottom_chord; x_top_chord] [y_bottom_chord; y_top_chord]]
                 
section_properties = [A_top_chord, A_bottom_chord, A_diagonal]

material_properties = [E_steel, E_concrete] 


#Initialize members array.
members = Array{NTuple{4,Int64},1}(undef, num_elem)

#Define bottom chord nodal connectivity, section properties, and material properties for bottom chord elements.
for i = 1:num_bottom_chord_elements
   
    members[i] = (i, i+1, 2, 2)

end

#Define top chord nodal connectivity, section properties, and material properties.
for i = 1:num_top_chord_elements
   
    index = num_bottom_chord_elements + i

    members[index] = (index + 1, index + 2, 1, 2)

end

#Define top diagonals nodal connectivity, section properties, and material properties for diagonal element pairs.

index = num_bottom_chord_elements + num_top_chord_elements + 1

for i = 1:(num_bottom_chord_nodes - 2)
   
    #Define the first diagonal element coming off of the bottom chord.
    node_i = i + 1
    node_j = num_bottom_chord_nodes + i
    
    members[index] = (node_i, node_j, 3, 1)

    index = index + 1

    #Define the second diagonal element coming off of the bottom chord.
    node_i = i + 1
    node_j = num_bottom_chord_nodes + i + 1 
    
    members[index] = (node_i, node_j, 3, 1)

    index = index + 1

end

#Define the degrees of freedom that are fixed.
supports = zeros(Int64, num_nodes * 2)

#Left arch springing point, bottom chord
supports[1] = 1   #x direction
supports[2] = 1   #y direction

#Left arch springing point, top chord
supports[num_bottom_chord_nodes * 2 + 1] = 1 #x direction
supports[num_bottom_chord_nodes * 2 + 2] = 1 #y direction

#Right arch springing point, bottom chord
supports[num_bottom_chord_nodes * 2 - 1] = 1   #x direction
supports[num_bottom_chord_nodes * 2] = 1   #y direction

#Right arch springing point, top chord
supports[num_nodes * 2 - 1] = 1 #x direction
supports[num_nodes * 2] = 1 #y direction

external_forces = zeros(Float64, num_nodes * 2)
#Apply a downward force at the middle of the arch.
external_forces[(num_bottom_chord_nodes + 19) * 2] = -10000.0 #N, y-direction

using StructuresKit

truss = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces)

truss = Truss.solve(truss)