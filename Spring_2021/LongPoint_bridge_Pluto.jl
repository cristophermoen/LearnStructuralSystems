### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 0a029b9c-708f-44b6-bcc5-e70d01b9f723
using Plots

# ╔═╡ 22efeb6e-ac97-417d-aaac-debfd7d2369c
using StructuresKit

# ╔═╡ 08b08d8d-b0d7-4309-a190-69a9c9ce12da
using GraphRecipes

# ╔═╡ bb2c2ab4-923c-11eb-084a-01e34131f9eb
md""" EN.560.301 Structural Systems I

Spring 2021

"""


# ╔═╡ 5ed85b6a-924e-11eb-02ae-af70f433eef8
md"""Perform a structural analysis of the LongPoint Bridge, a timber pedestrian bridge concept connecting Brooklyn and Queens."""

# ╔═╡ 9a838252-603e-4d4a-9116-d3adf6c4ac79
L = 125.0   #half the span length

# ╔═╡ 46a18b3b-3f44-4a60-904d-f72c14c5e02d
h = 40.0   #tower height

# ╔═╡ 583fcc57-a283-4f58-9adf-8cef61a09080
m = h/(L - L/3) #top chord slope, center span

# ╔═╡ db031a8d-8f43-4f9d-ad59-08c2e16f3f98
left_half_nodes = [0.0 0.0;
	             L/3 0.0;
	             2L/3 0.0;
	             2L/3+L/9 0.0;
				 2L/3+2L/9 0.0;
	             2L/3+3L/9 0.0;
	             L/6 h/2;
	             L/3 h;
	             L/3+L/6 h-m*(L/6);
				 2L/3+L/18 h-m*(L/3+L/18);
				 2L/3+2L/9-L/18 h-m*(L/3+2L/9-L/18);
	             2L/3+2L/9+L/18 h-m*(L/3+2L/9+L/18)]
	             

# ╔═╡ 9b67f085-0fb5-4108-aeeb-df99db1e461e
right_half_nodes = [2*L .- left_half_nodes[:,1] left_half_nodes[:, 2]]

# ╔═╡ 730e96c9-0cc2-4282-9c98-2c4ea0429a86
right_half_nodes[1:6, 1] = reverse(right_half_nodes[1:6, 1])

# ╔═╡ 84d0e238-7b86-44be-9199-c921c742b2b4
right_half_nodes[7:12, 1] = reverse(right_half_nodes[7:12, 1])

# ╔═╡ 8dab2709-a8ac-471c-b3e6-d5d9c9156f6c
right_half_nodes[7:12, 2] = reverse(right_half_nodes[7:12, 2])

# ╔═╡ fd070741-4735-4a0a-befc-3620715504e7
right_half_nodes

# ╔═╡ 4df16b6d-9cb5-4d90-92fb-b699502c144a
left_half_nodes


# ╔═╡ ae56f8a1-741a-4dca-b17c-790398e00d56
bridge_node_geometry = [left_half_nodes[1:6, :]; right_half_nodes[2:6, :]; left_half_nodes[7:end,:]; right_half_nodes[7:end,:]] .* 12  #convert ft to in.

# ╔═╡ d621a754-2c72-469f-9383-3267e71ca912
pontoon_nodes = [[bridge_node_geometry[3, 1] -36]; [bridge_node_geometry[4, 1] -36]; [bridge_node_geometry[5, 1] -36]; [bridge_node_geometry[7, 1] -36]; [bridge_node_geometry[8, 1] -36]; [bridge_node_geometry[9, 1] -36]]

# ╔═╡ c1e552ab-072a-43e3-8c46-ad39722c26b7
node_geometry = [bridge_node_geometry; pontoon_nodes]

# ╔═╡ 96db487a-421f-454a-9104-0f70e7ca46a6
plot(node_geometry[:, 1], node_geometry[:, 2], seriestype = scatter)

# ╔═╡ dbf428c5-f2e6-4801-a86d-1396f3ff258a
members = [(1, 2, 1, 1),
			(2, 3, 1, 1),
	       (3, 4, 1, 1),
			(4, 5, 1, 1),
			(5, 6, 1, 1),
			(6, 7, 1, 1),
	        (7, 8, 1, 1),
			(8, 9, 1, 1),
	        (9, 10, 1, 1),
	        (10, 11, 1, 1),
	        (1, 12, 2, 1),
	        (12, 13, 2, 1),
			(13, 14, 2, 1),
			(14, 15, 2, 1),
			(15, 16, 2, 1),
	        (16, 17, 2, 1),
	        (17, 6, 2, 1),
			(6, 18, 2, 1),
			(18, 19, 2, 1),
			(19, 20, 2, 1),
			(20, 21, 2, 1),
			(21, 22, 2, 1),
			(22, 23, 2, 1),
			(23, 11, 2, 1),
			(2, 12, 3, 1),
	        (2, 13, 3, 1),
	        (2, 14, 3, 1),
	        (3, 14, 3, 1),
	        (3, 15, 3, 1),
	        (4, 15, 3, 1),
	        (4, 16, 3, 1),
	        (5, 16, 3, 1),
	        (5, 17, 3, 1),
			(7, 18, 3, 1),
	        (7, 19, 3, 1),
	        (8, 19, 3, 1),
	        (8, 20, 3, 1),
	        (9, 20, 3, 1),
	        (9, 21, 3, 1),
	        (10, 21, 3, 1),
	        (10, 22, 3, 1),
	        (10, 23, 3, 1),
            (24, 3, 6, 3),
            (25, 4, 6, 3),
            (26, 5, 6, 3),
            (27, 7, 6, 3),
            (28, 8, 6, 3),
            (29, 9, 6, 3) ]

	
	        
		     

# ╔═╡ c79371ea-94f0-4004-9fde-b79df2f270f1
length(members)

# ╔═╡ ab66669b-c7e9-4037-a2a9-a3f77188441e
md"""Consider the steel mast."""

# ╔═╡ f2834f11-84a5-41c5-bca2-0fd013c88c96
mast_outer_diameter = 2 *12 # in.

# ╔═╡ b915c99c-e653-462f-80a4-fd64a1b0f0e1
mast_wall_thickness = 0.5 #in.

# ╔═╡ de88029a-1e6a-46b9-8d9c-fd28aee964c5
A_mast = pi * (mast_outer_diameter/2)^2 - pi * (mast_outer_diameter/2 - mast_wall_thickness)^2   #in^2

# ╔═╡ 953b9025-0874-409d-b6c8-51c1f7107208
md""" Add floating pontoon supports to model."""

# ╔═╡ c690dc94-8580-4745-87f0-a6abd4b4087a
L_pontoon = 15 * 12 #in.

# ╔═╡ 03382882-c427-4222-a5ef-2050d2c1f10b
w_pontoon = 3 * 12 #in.

# ╔═╡ da590853-bd78-4158-be21-4a55c385575e
h_pontoon = 2 * 12 #in.

# ╔═╡ 7fbb35ac-2b65-49fb-9931-4701a8d4d3a5
γ_water = 62.4 / 12^3 / 1000 #kips/in^3

# ╔═╡ ded542aa-303e-42bf-8944-24b1aa5cc856
k_pontoon = L_pontoon * w_pontoon * γ_water   #kips/in.

# ╔═╡ cd33a935-b766-45d0-93d8-7f51adf8ca71
max_pontoon_force = k_pontoon * h_pontoon  #kips

# ╔═╡ 8bf9ea21-8fdf-4892-928a-1f44e2e47892
water_surface_offset = 36.0  #in.

# ╔═╡ 93165c27-f0fc-4611-9035-087c2682e801
A_pontoon_spring = k_pontoon * water_surface_offset   #in^2, 

# ╔═╡ 85538480-d528-4824-b64e-c7d4a0d10630
section_properties = [10.0, 20.0, 30.0, A_mast, A_pontoon_spring/2, A_pontoon_spring]   #in^2

# ╔═╡ 067142ec-308f-4490-aa4e-7e76a90f680b
material_properties = [1600.0, 29000.0, 1.0]  #ksi

# ╔═╡ 1fa002a0-b43b-4e8b-b9fc-098c912503f4
num_nodes = size(node_geometry)[1]

# ╔═╡ 3a05eab0-da7e-4453-832f-6fe0a3e24dbc
external_forces = zeros(num_nodes*2)

# ╔═╡ 5fbffe0b-f19b-4801-ba50-e91646b245fa
bridge_width = 15*12  #in.

# ╔═╡ 39ee3a94-ccdb-45bc-a79a-63c80a405dbc
p_AASHTO_ped = 90/1000/144 #ksi

# ╔═╡ 33d985b4-d74c-434f-b227-fde4cb48d93c
w_AASHTO_ped = p_AASHTO_ped * bridge_width  #kips/in.

# ╔═╡ 07a43f08-c41d-4395-96f2-2858b7ff70f9
deck_element_lengths = diff(node_geometry[1:11,1])

# ╔═╡ e1535f10-7edb-4bd1-8a0b-071d933cb51a
trib_width = zeros(Float64, 11)

# ╔═╡ 2855e2e1-f92d-4b81-9e6a-5aaf086ad7b6
begin
	for i = 1:11
		if i==1
			trib_width[i] = deck_element_lengths[1]/2
		elseif i==11
			trib_width[i] = deck_element_lengths[end]/2
		else
			trib_width[i] = deck_element_lengths[i-1]/2 + deck_element_lengths[i]/2
		end
	end
end

# ╔═╡ 3595a2c5-e7c7-4cec-853c-c4969c2ad977
trib_width

# ╔═╡ 0d9887a9-f6b3-438a-a2c7-2d4ccd9832cb
ped_node_forces = w_AASHTO_ped .* trib_width #kips

# ╔═╡ 4ac690c6-771d-44e4-a2b1-f143afdd2b50
external_forces[2:2:22] .= -ped_node_forces

# ╔═╡ 925e3176-4f4f-49e5-80c7-00a4ad7a5554
supports = zeros(Int, num_nodes*2)

# ╔═╡ ea18a406-d9b5-436a-a4d2-e60fd9cd1fe2
supports[[1, 2, 22, 26, 44, 24*2, 25*2, 26*2, 27*2, 28*2, 29*2, (24*2-1), (25*2-1), (26*2-1), (27*2-1), (28*2-1), (29*2-1)]] .= 1

# ╔═╡ e646246b-9ccc-4c52-98c5-f2c30711d82a
begin

	bridge = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces)
	
	bridge = Truss.solve(bridge)
	
end

# ╔═╡ 6f628faa-f4bb-4427-9e60-44e531e5e825
bridge.u[12]

# ╔═╡ 01143829-d3ab-4978-811e-8bac8dd3f5b5
bridge.f_element[end-5]

# ╔═╡ af383036-798f-4468-a189-bc7369d7ed17
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
          linewidth=0.5,
	      arrowlengthfrac=0)
	
	graphplot!(g,
          x=node_geometry[:,1] .+ δx .* scale_x, y=node_geometry[:,2] .+ δy .* 	scale_y,
          nodeshape=:circle, nodesize=2,
          axis_buffer=0.01,
          curves=false,
          color=:black,
          linewidth=2,
	      arrowlengthfrac=0)
	
	quiver!(node_geometry[:,1],node_geometry[:,2],quiver=(truss_object.external_forces[1:2:end]./10, truss_object.external_forces[2:2:end]./10))
	
	return current()
	
end

# ╔═╡ cd7e2974-4520-4a9e-9d42-e3b3e2d1f224
δx = bridge.u[1:2:end]

# ╔═╡ 8d3250ca-7b5b-4924-9634-51b1cee236fc
δy = bridge.u[2:2:end]

# ╔═╡ 408fbbf1-8758-4f03-97e4-a295f1f99634
scale_x = 1.0

# ╔═╡ bcdafcf1-42a8-4cf7-b1b6-33dfad4bd862
scale_y = 1.0

# ╔═╡ 0d8d8f6a-5af8-467a-abfd-8bb6cc07621a
plot_truss(bridge, node_geometry, members, scale_x, scale_y)

# ╔═╡ 78549747-e431-4840-86a4-5a63fac2e2ff
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

# ╔═╡ 36aa3a8f-549c-4d79-b3ce-b0ed11d75034
plot_axial_force(bridge, node_geometry, members, 5)

# ╔═╡ Cell order:
# ╠═bb2c2ab4-923c-11eb-084a-01e34131f9eb
# ╠═5ed85b6a-924e-11eb-02ae-af70f433eef8
# ╠═9a838252-603e-4d4a-9116-d3adf6c4ac79
# ╠═46a18b3b-3f44-4a60-904d-f72c14c5e02d
# ╠═583fcc57-a283-4f58-9adf-8cef61a09080
# ╠═db031a8d-8f43-4f9d-ad59-08c2e16f3f98
# ╠═0a029b9c-708f-44b6-bcc5-e70d01b9f723
# ╠═9b67f085-0fb5-4108-aeeb-df99db1e461e
# ╠═730e96c9-0cc2-4282-9c98-2c4ea0429a86
# ╠═84d0e238-7b86-44be-9199-c921c742b2b4
# ╠═8dab2709-a8ac-471c-b3e6-d5d9c9156f6c
# ╠═fd070741-4735-4a0a-befc-3620715504e7
# ╠═4df16b6d-9cb5-4d90-92fb-b699502c144a
# ╠═ae56f8a1-741a-4dca-b17c-790398e00d56
# ╠═d621a754-2c72-469f-9383-3267e71ca912
# ╠═c1e552ab-072a-43e3-8c46-ad39722c26b7
# ╠═96db487a-421f-454a-9104-0f70e7ca46a6
# ╠═dbf428c5-f2e6-4801-a86d-1396f3ff258a
# ╠═c79371ea-94f0-4004-9fde-b79df2f270f1
# ╠═ab66669b-c7e9-4037-a2a9-a3f77188441e
# ╠═f2834f11-84a5-41c5-bca2-0fd013c88c96
# ╠═b915c99c-e653-462f-80a4-fd64a1b0f0e1
# ╠═de88029a-1e6a-46b9-8d9c-fd28aee964c5
# ╠═953b9025-0874-409d-b6c8-51c1f7107208
# ╠═c690dc94-8580-4745-87f0-a6abd4b4087a
# ╠═03382882-c427-4222-a5ef-2050d2c1f10b
# ╠═da590853-bd78-4158-be21-4a55c385575e
# ╠═7fbb35ac-2b65-49fb-9931-4701a8d4d3a5
# ╠═ded542aa-303e-42bf-8944-24b1aa5cc856
# ╠═cd33a935-b766-45d0-93d8-7f51adf8ca71
# ╠═8bf9ea21-8fdf-4892-928a-1f44e2e47892
# ╠═93165c27-f0fc-4611-9035-087c2682e801
# ╠═85538480-d528-4824-b64e-c7d4a0d10630
# ╠═067142ec-308f-4490-aa4e-7e76a90f680b
# ╠═1fa002a0-b43b-4e8b-b9fc-098c912503f4
# ╠═3a05eab0-da7e-4453-832f-6fe0a3e24dbc
# ╠═5fbffe0b-f19b-4801-ba50-e91646b245fa
# ╠═39ee3a94-ccdb-45bc-a79a-63c80a405dbc
# ╠═33d985b4-d74c-434f-b227-fde4cb48d93c
# ╠═07a43f08-c41d-4395-96f2-2858b7ff70f9
# ╠═e1535f10-7edb-4bd1-8a0b-071d933cb51a
# ╠═2855e2e1-f92d-4b81-9e6a-5aaf086ad7b6
# ╠═3595a2c5-e7c7-4cec-853c-c4969c2ad977
# ╠═0d9887a9-f6b3-438a-a2c7-2d4ccd9832cb
# ╠═4ac690c6-771d-44e4-a2b1-f143afdd2b50
# ╠═925e3176-4f4f-49e5-80c7-00a4ad7a5554
# ╠═ea18a406-d9b5-436a-a4d2-e60fd9cd1fe2
# ╠═22efeb6e-ac97-417d-aaac-debfd7d2369c
# ╠═e646246b-9ccc-4c52-98c5-f2c30711d82a
# ╠═6f628faa-f4bb-4427-9e60-44e531e5e825
# ╠═01143829-d3ab-4978-811e-8bac8dd3f5b5
# ╠═08b08d8d-b0d7-4309-a190-69a9c9ce12da
# ╠═af383036-798f-4468-a189-bc7369d7ed17
# ╠═cd7e2974-4520-4a9e-9d42-e3b3e2d1f224
# ╠═8d3250ca-7b5b-4924-9634-51b1cee236fc
# ╠═408fbbf1-8758-4f03-97e4-a295f1f99634
# ╠═bcdafcf1-42a8-4cf7-b1b6-33dfad4bd862
# ╠═0d8d8f6a-5af8-467a-abfd-8bb6cc07621a
# ╠═78549747-e431-4840-86a4-5a63fac2e2ff
# ╠═36aa3a8f-549c-4d79-b3ce-b0ed11d75034
