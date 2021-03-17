### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 4968508c-859a-11eb-24d9-ade066fdd136
using StructuresKit

# ╔═╡ 5c43cba4-85a5-11eb-0d4b-e3306ed3cbb0
using Plots

# ╔═╡ b1cd585c-8103-11eb-24c8-bfa388e7c514
function parabola(h, L, x)

    b = 4 .* h / L
    a = -b ./ L
    c = 0.0

    y = a .* x.^2 + b .* x .+ c

    dy_dx = 2 .* a .* x .+ b

    return y, dy_dx

end

# ╔═╡ c72f0fcc-8103-11eb-38ae-5d1a5ca121a7
#Define the arch bottom chord centerline span.
L_bottom_chord = 73000 #mm

# ╔═╡ c6bb296c-8103-11eb-1bf4-e7fcb8e684d0
#Define the arch bottom chord centerline height.
h_bottom_chord = 21000 + 350/2 - 4000 #mm

# ╔═╡ d791d496-8103-11eb-0d20-5ba0e7ca81a8
#Define the number of equal bottom chord arch segments.
num_bottom_chord_arch_segments = 36

# ╔═╡ dfa5032e-8103-11eb-153d-85e4ace5f809
#Define the centerline distance between top chord and bottom chord of arch.
chord_spacing = (22800 - 470/2) - (21000 - 350/2)

# ╔═╡ f5475772-8103-11eb-0b92-b96dfac08b30
begin
#We need the arch slope at springing point of the arch.
	y_arch, dy_dx_arch = parabola(h_bottom_chord, L_bottom_chord, 0.0)
	θ_springing_point = atan(dy_dx_arch[1])
end


# ╔═╡ 0245c26c-8104-11eb-260e-01b030d6bbf7
delta_x_bottom_to_top_chord = chord_spacing * cos(θ_springing_point)

# ╔═╡ 06c40dc4-8104-11eb-1649-930bfe1a8663
delta_y_bottom_to_top_chord = chord_spacing * sin(θ_springing_point)

# ╔═╡ 1c2de6da-8104-11eb-2f9e-914fd3b047dd
L_top_chord = L_bottom_chord + 2 * delta_x_bottom_to_top_chord 

# ╔═╡ 209f686c-8104-11eb-392f-255f2f0c1413
#Define the arch top chord centerline height.
h_top_chord = 22800 - 470/2 - 4000 #mm

# ╔═╡ 25c4ac88-8104-11eb-02ff-adb7f1fdf70e
#Define the number of equal top chord arch segments.
num_top_chord_arch_segments = 36

# ╔═╡ 38f3cd20-8104-11eb-2247-69b7f1163402
begin
# DISCRETIZE BOTTOM ARCH CHORD

#Calculate the bottom arch centerline arc length.  


	x = 0:L_bottom_chord/200:L_bottom_chord
	y, dy_dx = parabola(h_bottom_chord, L_bottom_chord, x)
	dx = diff(x)
	dy = diff(y)
	ds = sqrt.(dx.^2 + dy.^2)

	s_bottom_chord = sum(ds)  #This is the length of the bottom chord. 
end

# ╔═╡ 7e2e109e-8104-11eb-245f-8d15ee9b20a0
#Calculate the typical length of a bottom chord segment between truss diagonals.
length_bottom_chord_segment = s_bottom_chord / num_bottom_chord_arch_segments

# ╔═╡ 87e8b846-8104-11eb-1038-894f04a0727e
#Define all the bottom chord segments along the bottom chord.
bottom_chord_arch_layout = [length_bottom_chord_segment/2 ; ones(Float64, num_bottom_chord_arch_segments - 1) .* length_bottom_chord_segment; length_bottom_chord_segment/2]

# ╔═╡ a4f4fc56-8104-11eb-17c1-4756a521364b

function calculate_arch_xy(arch_layout, h, L)
	
    x_index = 0.0

    num_segments = length(arch_layout)

    dx_segment = zeros(Float64, num_segments)
    
    num_nodes = num_segments + 1
    x_arch = zeros(Float64, num_nodes)
    y_arch = zeros(Float64, num_nodes)

	for i = 1:num_segments
	
		arch_segment_length = arch_layout[i]
		
        y, dy_dx = parabola(h, L, x_arch[i])
        
        θ = atan(dy_dx)
	
		dx_segment[i] = arch_segment_length * cos(θ)
		
        x_index = x_index + dx_segment[i]
        
        x_arch[i+1] = x_index
		
	end

    y_arch, dy_dx = parabola(h, L, x_arch)

	return x_arch, y_arch
	
end

# ╔═╡ af780a64-8105-11eb-1ec4-3118052de520
x_bottom_chord, y_bottom_chord = calculate_arch_xy(bottom_chord_arch_layout, h_bottom_chord, L_bottom_chord)

# ╔═╡ 6e8ab964-8107-11eb-1c57-01a7a1fd710a
# DISCRETIZE TOP ARCH CHORD
begin 
x_top_chord_disc = 0:L_top_chord/200:L_top_chord
y_top_chord_disc, dy_dx_top_chord = parabola(h_top_chord, L_top_chord, x_top_chord_disc)
dx_top_chord = diff(x_top_chord_disc)
dy_top_chord = diff(y_top_chord_disc)
ds_top_chord = sqrt.(dx_top_chord.^2 + dy_top_chord.^2)

s_top_chord = sum(ds_top_chord)  #This is the length of the top chord.
	
end

# ╔═╡ b0b7559c-8107-11eb-1669-3946daddfb7d
#Calculate the typical length of a top chord segment between truss diagonals.
length_top_chord_segment = s_top_chord / num_top_chord_arch_segments

# ╔═╡ b75e4dd6-8107-11eb-229a-d721ae3f3674
#Define all the top chord segments along the top chord.
top_chord_arch_layout = ones(Float64, num_top_chord_arch_segments) * length_top_chord_segment

# ╔═╡ be32003a-8107-11eb-1e0d-c9c96c1b75e1
begin
x_top_chord, y_top_chord = calculate_arch_xy(top_chord_arch_layout, h_top_chord, L_top_chord)

#Shift the top chord x-coordinates.
x_top_chord = x_top_chord .- delta_x_bottom_to_top_chord 

#Shift the top chord y-coordinates.
y_top_chord = y_top_chord .+ delta_y_bottom_to_top_chord 
end

# ╔═╡ 4fe85858-8108-11eb-2927-734257192557
begin
#CALCULATE SECTION PROPERTIES 

A_top_chord = 470.0 * 300.0   
A_bottom_chord = 350.0 * 350.0
A_diagonal = 100.0 * 100 - (100.0 - 5.0) * (100.0 - 5.0)
end

# ╔═╡ 94915568-8108-11eb-0ba2-8b316afb9495
begin
#DEFINE MATERIAL PROPERTIES

E_steel = 199947 #N/mm^2
E_concrete = 27579 #N/mm^2   #4000 ksi
	
end

# ╔═╡ b0645880-8108-11eb-30a2-8997109ff992
begin
	
#DEFINE DISCRETIZATION INFO
num_top_chord_nodes = length(x_top_chord)
num_bottom_chord_nodes = length(x_bottom_chord)
num_nodes = num_top_chord_nodes + num_bottom_chord_nodes

num_top_chord_elements = length(top_chord_arch_layout)
num_bottom_chord_elements = length(bottom_chord_arch_layout)
num_diagonal_elements = (num_bottom_chord_nodes - 2) * 2
num_elem = num_top_chord_elements + num_bottom_chord_elements + num_diagonal_elements
end

# ╔═╡ be7530c8-8108-11eb-3532-ff6e2e7614ab
node_geometry = [[x_bottom_chord; x_top_chord] [y_bottom_chord; y_top_chord]]

# ╔═╡ dd79cbb6-8108-11eb-3da5-d9e6fcd8fc68
section_properties = [A_top_chord, A_bottom_chord, A_diagonal]

# ╔═╡ e83eac5e-8108-11eb-19ca-0961c8cfd596
material_properties = [E_steel, E_concrete] 

# ╔═╡ 0166b12c-8109-11eb-005b-35f6ca83637b

#Initialize members array.
members = Array{NTuple{4,Int64},1}(undef, num_elem)

# ╔═╡ 0b46f78a-8109-11eb-139f-fb048c62638a
#Define bottom chord nodal connectivity, section properties, and material properties for bottom chord elements.
for i = 1:num_bottom_chord_elements
   
    members[i] = (i, i+1, 2, 2)

end

# ╔═╡ 662ee98c-8109-11eb-202a-0f72148a69aa
#Define top chord nodal connectivity, section properties, and material properties.
for i = 1:num_top_chord_elements
   
    index = num_bottom_chord_elements + i

    members[index] = (index + 1, index + 2, 1, 2)

end

# ╔═╡ 956803dc-8109-11eb-3dcf-839736da0635
begin

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
	
end

# ╔═╡ ed0636d8-8109-11eb-0258-adc809c96b6b

#Define the degrees of freedom that are fixed.
supports = zeros(Int64, num_nodes * 2)

# ╔═╡ 96e6fa82-810a-11eb-07bf-d3a2424cb786
begin
#Left arch springing point, bottom chord
supports[1] = 1   #x direction
supports[2] = 1   #y direction
end

# ╔═╡ b693a4fc-810a-11eb-0da0-7157b16b24c9
begin
#Left arch springing point, top chord
supports[num_bottom_chord_nodes * 2 + 1] = 1 #x direction
supports[num_bottom_chord_nodes * 2 + 2] = 1 #y direction
end

# ╔═╡ c84290e6-810a-11eb-235b-b9c6cbaac6d2
begin
#Right arch springing point, bottom chord
supports[num_bottom_chord_nodes * 2 - 1] = 1   #x direction
supports[num_bottom_chord_nodes * 2] = 1   #y direction
end

# ╔═╡ f820eaa6-810a-11eb-033e-83006c662ee5
begin
#Right arch springing point, top chord
supports[num_nodes * 2 - 1] = 1 #x direction
supports[num_nodes * 2] = 1 #y direction
end

# ╔═╡ 3904c5a6-810b-11eb-257d-410d94d98734
p_windward = 663.97/1000^2 #N/mm^2

# ╔═╡ 6aa671e2-810b-11eb-3ac6-855bf0b3b9f0
p_center = -1537.62/1000^2 #N/mm^2

# ╔═╡ 7f4259d4-810b-11eb-1672-771092efea9e
p_leeward = 759.69/1000^2  #N/mm^2

# ╔═╡ ab76bd56-810b-11eb-2992-5d5fba4177cf
trib_width = 8260. #mm

# ╔═╡ 91d839b8-810b-11eb-1ac1-a335bb8e7e7e
w_windward = p_windward * trib_width #N/mm

# ╔═╡ 352e4c58-81b1-11eb-2449-cd5c483cdfc2
w_center = p_center * trib_width

# ╔═╡ 3cee1400-81b1-11eb-1606-1df75151ba01
w_leeward = p_leeward * trib_width

# ╔═╡ 4a36bee6-81b1-11eb-02fd-fbbc7fe65b74
windward_node_index = findall(x->x<=L_top_chord/4, x_top_chord) 

# ╔═╡ f220e9a0-81b2-11eb-2f07-61b8f791d854
windward_node_range = windward_node_index .+ num_bottom_chord_nodes

# ╔═╡ 9b010cfa-81b1-11eb-24cf-85f003a6049e
center_node_index = findall(x->((x>L_top_chord/4) & (x<3 * L_top_chord/4)), x_top_chord) 

# ╔═╡ 07f6a288-81b3-11eb-114b-49cb5830f762
center_node_range = center_node_index .+ num_bottom_chord_nodes

# ╔═╡ f4804bce-81b1-11eb-2632-5328c3f417aa
leeward_node_index = findall(x->x> 3*L_top_chord/4, x_top_chord)

# ╔═╡ 1946fc72-81b3-11eb-24ca-f766b8e80c1f
leeward_node_range = leeward_node_index .+ num_bottom_chord_nodes

# ╔═╡ 894f99d6-810c-11eb-0b0e-1bdc56a2de41
P_windward = w_windward * length_top_chord_segment   #N

# ╔═╡ 7dfa2ab4-81b2-11eb-1f6a-29e4eeb4e60f
P_center = w_center * length_top_chord_segment  #N

# ╔═╡ 8d4120ae-81b2-11eb-0a23-45b509d1a40b
P_leeward = w_leeward * length_top_chord_segment #N

# ╔═╡ 6b63f3e8-81b3-11eb-21e6-fde00c5dcc29
P_wind = zeros(Float64, num_top_chord_nodes)

# ╔═╡ 9f0d935c-81b3-11eb-02da-ff1cd1751b6b
P_wind[windward_node_index] .= P_windward

# ╔═╡ fdebaabc-81b3-11eb-16c2-37a14a0f89c4
P_wind[center_node_index] .= P_center

# ╔═╡ 08e34eca-81b4-11eb-2de5-df2c13ab1bfb
P_wind[leeward_node_index] .= P_leeward

# ╔═╡ 1a6cec78-81b4-11eb-2434-f109162c84c3
not_used, dy_dx_top_chord_nodes = parabola(h_top_chord, L_top_chord, x_top_chord)

# ╔═╡ 7c5b9c86-81b4-11eb-0f18-652107f741b9
top_chord_angle = atan.(dy_dx_top_chord_nodes)

# ╔═╡ 6c9273ba-81b4-11eb-12a1-53d45f118b65
P_wind_x = P_wind .* cos.(π/2 .- top_chord_angle)

# ╔═╡ a5ea598e-81b4-11eb-1da8-bbba07b26ded
P_wind_y = -P_wind .* sin.(π/2 .- top_chord_angle)  #note negative sign here

# ╔═╡ 34dbd85e-810e-11eb-02b3-e928e9bcbc7d
external_forces = zeros(Float64, num_nodes * 2)

# ╔═╡ aaeb3b64-81b5-11eb-3bd7-739827b62939
external_forces[(windward_node_range .- 1) .* 2 .+ 1] .= P_wind_x[windward_node_index]

# ╔═╡ f79e4adc-8732-11eb-3dc1-735634b7879e
external_forces

# ╔═╡ 5528de4c-81b6-11eb-1ba4-5d5b816f3804
external_forces[(windward_node_range .- 1) .* 2 .+ 2] .= P_wind_y[windward_node_index]

# ╔═╡ 65cd8c2a-81b6-11eb-0a1d-cf523d856b8a
external_forces[(center_node_range .- 1) .* 2 .+ 1] .= P_wind_x[center_node_index]

# ╔═╡ 8388b212-81b6-11eb-3cc2-116509bde9b9
external_forces[(center_node_range .- 1) .* 2 .+ 2] .= P_wind_y[center_node_index]

# ╔═╡ 8c877e02-81b6-11eb-1a62-c97efbadf078
external_forces[(leeward_node_range .- 1) .* 2 .+ 1] .= P_wind_x[leeward_node_index]

# ╔═╡ 9b69ad5a-81b6-11eb-373c-192e0a6cf9d2
external_forces[(leeward_node_range .- 1) .* 2 .+ 2] .= P_wind_y[leeward_node_index]

# ╔═╡ cf2d7d30-8731-11eb-1498-4967ea7caa67
external_forces

# ╔═╡ 028f0e26-81b7-11eb-1bdb-1f93f33adbe3
begin
CASA = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces)
	
CASA = Truss.solve(CASA)
	
end

# ╔═╡ 50ac3eb4-872f-11eb-35eb-571cf2a42cc5
function plot_truss(truss_object, node_geometry, members, scale_x, scale_y)
	
	
	δx = truss_object.u[1:2:end]
	δy = truss_object.u[2:2:end]
	
	num_nodes = floor(Int, length(truss_object.u)/2)
	
	g = zeros(Int64, (num_nodes, num_nodes))
	
	start_nodes = [x[1] for x in members]
	end_nodes = [x[2] for x in members]

	for i = 1:num_nodes
	
		index = findall(x->x==i, start_nodes)
	
		g[i, end_nodes[index]] .= 1
	
	end
	
	graphplot(g,
          x=node_geometry[:,1], y=node_geometry[:,2],
          nodeshape=:circle, nodesize=2,
          axis_buffer=0.01,
          curves=false,
          color=:red,
          linewidth=1)
	
	graphplot!(g,
          x=node_geometry[:,1] .+ δx .* scale_x, y=node_geometry[:,2] .+ δy .* 	scale_y,
          nodeshape=:circle, nodesize=2,
          axis_buffer=0.01,
          curves=false,
          color=:black,
          linewidth=1)
	
	quiver!(node_geometry[:,1],node_geometry[:,2],quiver=(truss_object.external_forces[1:2:end]./10, truss_object.external_forces[2:2:end]./10))
	
	return current()
	
end
	

# ╔═╡ ae3b28c6-81b8-11eb-1c90-1fdccbc827b3
δx = CASA.u[1:2:end]

# ╔═╡ f46166ee-81b8-11eb-05be-a5d55e9dac2a
δy = CASA.u[2:2:end]

# ╔═╡ 0173266a-81b9-11eb-1343-e5fa5cd85298
scale_x = 200.0

# ╔═╡ 0c493566-81b9-11eb-00d1-afd586cb640f
scale_y = 200.0

# ╔═╡ 0a5888f2-81bd-11eb-2f3b-a504d4ae781f
#Define graph of arch elements and nodes.

# ╔═╡ 18f8a626-81bd-11eb-371d-8b09229fdeae
g = zeros(Int64, (num_nodes, num_nodes))

# ╔═╡ c9d71a14-81bc-11eb-09d6-976ba1d71ee0
begin
	using GraphRecipes

graphplot(g,
          x=node_geometry[:,1], y=node_geometry[:,2],
          nodeshape=:circle, nodesize=2,
          axis_buffer=0.01,
          curves=false,
          color=:red,
          linewidth=1)
	
graphplot!(g,
          x=node_geometry[:,1] .+ δx .* scale_x, y=node_geometry[:,2] .+ δy .* scale_y,
          nodeshape=:circle, nodesize=2,
          axis_buffer=0.01,
          curves=false,
          color=:black,
          linewidth=1)
	
	quiver!(node_geometry[:,1],node_geometry[:,2],quiver=(CASA.external_forces[1:2:end]./10, CASA.external_forces[2:2:end]./10))
	
end

# ╔═╡ 3d35e698-81bd-11eb-3a99-0f4d644cc8ee
begin
start_nodes = [x[1] for x in members]
end_nodes = [x[2] for x in members]

for i = 1:num_nodes
	
	index = findall(x->x==i, start_nodes)
	
	g[i, end_nodes[index]] .= 1
	
end
end

# ╔═╡ 46d71620-81bf-11eb-3a55-af708b8bf3d0
g

# ╔═╡ 80c15a46-85a0-11eb-21b1-279843cc2978
line_weight = ones(num_nodes)*3

# ╔═╡ fc89b1e6-85a0-11eb-256c-6fc163965faf
axial_force = [ x[3] for x in CASA.f_element]

# ╔═╡ 4e1dee06-85a8-11eb-3799-194980ccce1b
axial_force

# ╔═╡ 6f3abc6e-85a2-11eb-043b-cf3df249a404
normalized_axial_force = axial_force ./ maximum(abs.(axial_force))

# ╔═╡ a2da93be-85a2-11eb-2631-b3dee12a26e4
function plot_axial_force(truss_object, node_geometry, members, scale)
	
	axial_force = [ x[3] for x in truss_object.f_element]
	
	normalized_axial_force = axial_force ./ maximum(abs.(axial_force)) .* scale

	p = plot(node_geometry[:,1], node_geometry[:,2], seriestype = scatter, legend=false, size = (600,200))
	
	for i = 1:size(members)[1]
	
		node_i = members[i][1]
		node_j = members[i][2]
		
		if normalized_axial_force[i] <= 0.0
		
			plot!(p, [node_geometry[node_i, 1], node_geometry[node_j, 1]], [node_geometry[node_i, 2], node_geometry[node_j, 2]], color = :blue, linewidth = abs(normalized_axial_force[i]), legend=false)
			
		else
			
			plot!(p, [node_geometry[node_i, 1], node_geometry[node_j, 1]], [node_geometry[node_i, 2], node_geometry[node_j, 2]], color = :red, linewidth = abs(normalized_axial_force[i]), legend=false)

		end
		
	end
	
	return current()

end
	


# ╔═╡ 9895be68-85a6-11eb-216e-73f313b179bb
plot_axial_force(CASA, node_geometry, members, 5)

# ╔═╡ 783f3652-859d-11eb-0a9b-6f49dc81a0ec
md""" Calculate self weight of arch truss."""

# ╔═╡ 890395d4-859d-11eb-14a9-ebdf7374fe97
γ_concrete = 23563/1000^3 #N/mm^3

# ╔═╡ 6f427bda-859e-11eb-1a85-f1a095e5fae3
γ_steel = 77287/1000^3 #N/mm^3

# ╔═╡ 86504064-859e-11eb-198a-6dff40c3a926
W_top_chord = s_top_chord * A_top_chord * γ_concrete

# ╔═╡ a12efb1e-859e-11eb-3238-fdbe298d58da
W_bottom_chord = s_bottom_chord * A_bottom_chord * γ_concrete

# ╔═╡ b1379282-859e-11eb-0e7c-85f9c946e610
W_roof_shell = s_top_chord * 50 * trib_width/2 * γ_concrete + s_bottom_chord * 50 * trib_width/2 * γ_concrete

# ╔═╡ 20cbffc0-859f-11eb-21aa-c1f9e40a5bc4
s_diagonals = sum(CASA.L[(num_bottom_chord_elements+num_top_chord_elements+1):end])

# ╔═╡ 7a0151c6-859f-11eb-0b36-73cc5b5f8962
W_steel_diagonals = s_diagonals * A_diagonal * γ_steel

# ╔═╡ 902cbe0e-859f-11eb-1421-df0d95e2e39e
W_total = W_top_chord + W_bottom_chord + W_roof_shell + W_steel_diagonals

# ╔═╡ cbd910ce-85b3-11eb-1439-9d05fae9cc4a
w_self_weight = W_total/L_top_chord

# ╔═╡ 448727ec-872e-11eb-1caa-b9798941be52
P_self_weight = W_total/(num_bottom_chord_nodes + num_top_chord_nodes) #N

# ╔═╡ 7d4579e4-872e-11eb-3aac-8986c0e056e9
external_forces_self_weight = zeros(Float64, num_nodes * 2)

# ╔═╡ ac1fa276-872e-11eb-29f2-33da90c45532
external_forces_self_weight[2:2:end] .= -P_self_weight   #note negative sign, for gravity

# ╔═╡ 3f421abe-8736-11eb-27e8-a9cced210440
external_forces_self_weight

# ╔═╡ e8657170-872e-11eb-11cb-dd39b183a73f
begin
CASA_self_weight = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces_self_weight)
	
CASA_self_weight = Truss.solve(CASA_self_weight)
	
end

# ╔═╡ 568570e0-8736-11eb-3d1b-f309605350f1
CASA_self_weight.external_forces

# ╔═╡ 299a42fa-8730-11eb-3fb9-d3300ea99fc1
plot_truss(CASA_self_weight, node_geometry, members, 1000.0, 1000.0)

# ╔═╡ 9f6e3540-8730-11eb-3420-57d6b0fc2d87
maximum(abs.(CASA_self_weight.u[2:2:end]))

# ╔═╡ 7208e75c-8731-11eb-3f5a-53081c021f1d
plot_axial_force(CASA_self_weight, node_geometry, members, 5)

# ╔═╡ 98005832-8731-11eb-1fe0-f9cbcc04a0dd
md"""DL+W"""

# ╔═╡ a574199a-8731-11eb-28db-6d29a9d03c11
begin
	external_forces_DL_W = zeros(Float64, num_nodes * 2)
	external_forces_DL_W = CASA.external_forces .+ CASA_self_weight.external_forces
end

# ╔═╡ f4a38132-8735-11eb-036c-29088205f4ff


# ╔═╡ d17694ec-8735-11eb-06d4-179dc79b6ff6
CASA.external_forces

# ╔═╡ d7e504bc-8735-11eb-1ba6-d5ae58dbeb2d
CASA_self_weight.external_forces

# ╔═╡ 652653e8-8735-11eb-0fa9-09f380360819
begin
CASA_DL_W = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces_DL_W)
	
CASA_DL_W = Truss.solve(CASA_DL_W)
	
end

# ╔═╡ 75da8206-8735-11eb-142b-653d35df9e2d
plot_truss(CASA_DL_W, node_geometry, members, 200.0, 200.0)

# ╔═╡ 71bbca34-8737-11eb-0a01-c3f03351bb58
plot_axial_force(CASA_DL_W, node_geometry, members, 5)

# ╔═╡ 93bd5ae8-8738-11eb-159f-bd9cd82796a3
axial_forces_DL_W = [x[3] for x in CASA_DL_W.f_element]

# ╔═╡ ed445f58-8738-11eb-1f69-3390a0af3cc4
axial_forces_DL_W[1:37]

# ╔═╡ c83ce05e-8738-11eb-3cf7-5f07f7c7c242
max_bottom_chord_tension = maximum(axial_forces_DL_W[1:37])

# ╔═╡ 29b8a5fc-8739-11eb-1169-b5e8ce061630
max_bottom_chord_tensile_stress = max_bottom_chord_tension / A_bottom_chord

# ╔═╡ 6139b6f6-8739-11eb-14f3-bb79f9d4f149
A_bottom_chord

# ╔═╡ Cell order:
# ╠═b1cd585c-8103-11eb-24c8-bfa388e7c514
# ╠═c72f0fcc-8103-11eb-38ae-5d1a5ca121a7
# ╠═c6bb296c-8103-11eb-1bf4-e7fcb8e684d0
# ╠═d791d496-8103-11eb-0d20-5ba0e7ca81a8
# ╠═dfa5032e-8103-11eb-153d-85e4ace5f809
# ╠═f5475772-8103-11eb-0b92-b96dfac08b30
# ╠═0245c26c-8104-11eb-260e-01b030d6bbf7
# ╠═06c40dc4-8104-11eb-1649-930bfe1a8663
# ╠═1c2de6da-8104-11eb-2f9e-914fd3b047dd
# ╠═209f686c-8104-11eb-392f-255f2f0c1413
# ╠═25c4ac88-8104-11eb-02ff-adb7f1fdf70e
# ╠═38f3cd20-8104-11eb-2247-69b7f1163402
# ╠═7e2e109e-8104-11eb-245f-8d15ee9b20a0
# ╠═87e8b846-8104-11eb-1038-894f04a0727e
# ╠═a4f4fc56-8104-11eb-17c1-4756a521364b
# ╠═af780a64-8105-11eb-1ec4-3118052de520
# ╠═6e8ab964-8107-11eb-1c57-01a7a1fd710a
# ╠═b0b7559c-8107-11eb-1669-3946daddfb7d
# ╠═b75e4dd6-8107-11eb-229a-d721ae3f3674
# ╠═be32003a-8107-11eb-1e0d-c9c96c1b75e1
# ╠═4fe85858-8108-11eb-2927-734257192557
# ╠═94915568-8108-11eb-0ba2-8b316afb9495
# ╠═b0645880-8108-11eb-30a2-8997109ff992
# ╠═be7530c8-8108-11eb-3532-ff6e2e7614ab
# ╠═dd79cbb6-8108-11eb-3da5-d9e6fcd8fc68
# ╠═e83eac5e-8108-11eb-19ca-0961c8cfd596
# ╠═0166b12c-8109-11eb-005b-35f6ca83637b
# ╠═0b46f78a-8109-11eb-139f-fb048c62638a
# ╠═662ee98c-8109-11eb-202a-0f72148a69aa
# ╠═956803dc-8109-11eb-3dcf-839736da0635
# ╠═ed0636d8-8109-11eb-0258-adc809c96b6b
# ╠═96e6fa82-810a-11eb-07bf-d3a2424cb786
# ╠═b693a4fc-810a-11eb-0da0-7157b16b24c9
# ╠═c84290e6-810a-11eb-235b-b9c6cbaac6d2
# ╠═f820eaa6-810a-11eb-033e-83006c662ee5
# ╠═3904c5a6-810b-11eb-257d-410d94d98734
# ╠═6aa671e2-810b-11eb-3ac6-855bf0b3b9f0
# ╠═7f4259d4-810b-11eb-1672-771092efea9e
# ╠═ab76bd56-810b-11eb-2992-5d5fba4177cf
# ╠═91d839b8-810b-11eb-1ac1-a335bb8e7e7e
# ╠═352e4c58-81b1-11eb-2449-cd5c483cdfc2
# ╠═3cee1400-81b1-11eb-1606-1df75151ba01
# ╠═4a36bee6-81b1-11eb-02fd-fbbc7fe65b74
# ╠═f220e9a0-81b2-11eb-2f07-61b8f791d854
# ╠═9b010cfa-81b1-11eb-24cf-85f003a6049e
# ╠═07f6a288-81b3-11eb-114b-49cb5830f762
# ╠═f4804bce-81b1-11eb-2632-5328c3f417aa
# ╠═1946fc72-81b3-11eb-24ca-f766b8e80c1f
# ╠═894f99d6-810c-11eb-0b0e-1bdc56a2de41
# ╠═7dfa2ab4-81b2-11eb-1f6a-29e4eeb4e60f
# ╠═8d4120ae-81b2-11eb-0a23-45b509d1a40b
# ╠═6b63f3e8-81b3-11eb-21e6-fde00c5dcc29
# ╠═9f0d935c-81b3-11eb-02da-ff1cd1751b6b
# ╠═fdebaabc-81b3-11eb-16c2-37a14a0f89c4
# ╠═08e34eca-81b4-11eb-2de5-df2c13ab1bfb
# ╠═1a6cec78-81b4-11eb-2434-f109162c84c3
# ╠═7c5b9c86-81b4-11eb-0f18-652107f741b9
# ╠═6c9273ba-81b4-11eb-12a1-53d45f118b65
# ╠═a5ea598e-81b4-11eb-1da8-bbba07b26ded
# ╠═34dbd85e-810e-11eb-02b3-e928e9bcbc7d
# ╠═aaeb3b64-81b5-11eb-3bd7-739827b62939
# ╠═f79e4adc-8732-11eb-3dc1-735634b7879e
# ╠═5528de4c-81b6-11eb-1ba4-5d5b816f3804
# ╠═65cd8c2a-81b6-11eb-0a1d-cf523d856b8a
# ╠═8388b212-81b6-11eb-3cc2-116509bde9b9
# ╠═8c877e02-81b6-11eb-1a62-c97efbadf078
# ╠═9b69ad5a-81b6-11eb-373c-192e0a6cf9d2
# ╠═cf2d7d30-8731-11eb-1498-4967ea7caa67
# ╠═4968508c-859a-11eb-24d9-ade066fdd136
# ╠═028f0e26-81b7-11eb-1bdb-1f93f33adbe3
# ╠═50ac3eb4-872f-11eb-35eb-571cf2a42cc5
# ╠═ae3b28c6-81b8-11eb-1c90-1fdccbc827b3
# ╠═f46166ee-81b8-11eb-05be-a5d55e9dac2a
# ╠═0173266a-81b9-11eb-1343-e5fa5cd85298
# ╠═0c493566-81b9-11eb-00d1-afd586cb640f
# ╠═0a5888f2-81bd-11eb-2f3b-a504d4ae781f
# ╠═18f8a626-81bd-11eb-371d-8b09229fdeae
# ╠═3d35e698-81bd-11eb-3a99-0f4d644cc8ee
# ╠═46d71620-81bf-11eb-3a55-af708b8bf3d0
# ╠═c9d71a14-81bc-11eb-09d6-976ba1d71ee0
# ╠═80c15a46-85a0-11eb-21b1-279843cc2978
# ╠═fc89b1e6-85a0-11eb-256c-6fc163965faf
# ╠═4e1dee06-85a8-11eb-3799-194980ccce1b
# ╠═6f3abc6e-85a2-11eb-043b-cf3df249a404
# ╠═5c43cba4-85a5-11eb-0d4b-e3306ed3cbb0
# ╠═a2da93be-85a2-11eb-2631-b3dee12a26e4
# ╠═9895be68-85a6-11eb-216e-73f313b179bb
# ╠═783f3652-859d-11eb-0a9b-6f49dc81a0ec
# ╠═890395d4-859d-11eb-14a9-ebdf7374fe97
# ╠═6f427bda-859e-11eb-1a85-f1a095e5fae3
# ╠═86504064-859e-11eb-198a-6dff40c3a926
# ╠═a12efb1e-859e-11eb-3238-fdbe298d58da
# ╠═b1379282-859e-11eb-0e7c-85f9c946e610
# ╠═20cbffc0-859f-11eb-21aa-c1f9e40a5bc4
# ╠═7a0151c6-859f-11eb-0b36-73cc5b5f8962
# ╠═902cbe0e-859f-11eb-1421-df0d95e2e39e
# ╠═cbd910ce-85b3-11eb-1439-9d05fae9cc4a
# ╠═448727ec-872e-11eb-1caa-b9798941be52
# ╠═7d4579e4-872e-11eb-3aac-8986c0e056e9
# ╠═ac1fa276-872e-11eb-29f2-33da90c45532
# ╠═3f421abe-8736-11eb-27e8-a9cced210440
# ╠═e8657170-872e-11eb-11cb-dd39b183a73f
# ╠═568570e0-8736-11eb-3d1b-f309605350f1
# ╠═299a42fa-8730-11eb-3fb9-d3300ea99fc1
# ╠═9f6e3540-8730-11eb-3420-57d6b0fc2d87
# ╠═7208e75c-8731-11eb-3f5a-53081c021f1d
# ╠═98005832-8731-11eb-1fe0-f9cbcc04a0dd
# ╠═a574199a-8731-11eb-28db-6d29a9d03c11
# ╠═f4a38132-8735-11eb-036c-29088205f4ff
# ╠═d17694ec-8735-11eb-06d4-179dc79b6ff6
# ╠═d7e504bc-8735-11eb-1ba6-d5ae58dbeb2d
# ╠═652653e8-8735-11eb-0fa9-09f380360819
# ╠═75da8206-8735-11eb-142b-653d35df9e2d
# ╠═71bbca34-8737-11eb-0a01-c3f03351bb58
# ╠═93bd5ae8-8738-11eb-159f-bd9cd82796a3
# ╠═ed445f58-8738-11eb-1f69-3390a0af3cc4
# ╠═c83ce05e-8738-11eb-3cf7-5f07f7c7c242
# ╠═29b8a5fc-8739-11eb-1169-b5e8ce061630
# ╠═6139b6f6-8739-11eb-14f3-bb79f9d4f149
