### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

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

# ╔═╡ 18dc1602-810c-11eb-300f-d7dd63be8b35


# ╔═╡ 657229dc-810c-11eb-2457-47382dac4cc2
center_node_range = 49:49+18

# ╔═╡ 7636074a-810c-11eb-1bab-ab37276a377e
leeward_node_range = 68:75

# ╔═╡ 894f99d6-810c-11eb-0b0e-1bdc56a2de41
P_windward = w_windward * length_top_chord_segment

# ╔═╡ de9168de-810c-11eb-0064-cb0d9e164dbe
#for left springing point, first node
not_used, dy_dx_1 = parabola(h_top_chord, L_top_chord, 0.0)

# ╔═╡ 8093ec88-810d-11eb-24ab-bb36186ae7ac
face_angle = atan(dy_dx_1)

# ╔═╡ 34dbd85e-810e-11eb-02b3-e928e9bcbc7d
external_forces = zeros(Float64, num_nodes * 2)

# ╔═╡ b6bdbe06-810d-11eb-0b03-f55f968e6602
external_forces[38*2 + 1] = P_windward/2 * cos(face_angle)

# ╔═╡ ef4a7d68-810d-11eb-3f56-f93dd251f8c6
external_forces[38*2 + 2] = -P_windward/2 * sin(face_angle)

# ╔═╡ f2a279a4-810b-11eb-19c3-79e560d45b43
top_chord_arch_layout

# ╔═╡ f23a685a-810b-11eb-186b-259e129ca440


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
# ╠═18dc1602-810c-11eb-300f-d7dd63be8b35
# ╠═657229dc-810c-11eb-2457-47382dac4cc2
# ╠═7636074a-810c-11eb-1bab-ab37276a377e
# ╠═894f99d6-810c-11eb-0b0e-1bdc56a2de41
# ╠═de9168de-810c-11eb-0064-cb0d9e164dbe
# ╠═8093ec88-810d-11eb-24ab-bb36186ae7ac
# ╠═34dbd85e-810e-11eb-02b3-e928e9bcbc7d
# ╠═b6bdbe06-810d-11eb-0b03-f55f968e6602
# ╠═ef4a7d68-810d-11eb-3f56-f93dd251f8c6
# ╠═f2a279a4-810b-11eb-19c3-79e560d45b43
# ╠═f23a685a-810b-11eb-186b-259e129ca440
