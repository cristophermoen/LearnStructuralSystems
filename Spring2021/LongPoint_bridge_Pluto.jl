### A Pluto.jl notebook ###
# v0.14.0

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
node_geometry = [left_half_nodes[1:6, :]; right_half_nodes[2:6, :]; left_half_nodes[7:end,:]; right_half_nodes[7:end,:]] .* 12  #convert ft to in.

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
	        (10, 23, 3, 1)]

	
	        
		     

# ╔═╡ 85538480-d528-4824-b64e-c7d4a0d10630
section_properties = [10.0, 20.0, 30.0]   #in^2

# ╔═╡ 067142ec-308f-4490-aa4e-7e76a90f680b
material_properties = [1600]  #ksi

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
supports[[1, 2, 22, 26, 44]] .= 1

# ╔═╡ e646246b-9ccc-4c52-98c5-f2c30711d82a
begin

	bridge = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces)
	
	bridge = Truss.solve(bridge)
	
end

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

# ╔═╡ 15fd59f7-2866-44f4-a5ea-3472b4a5fa05
g = zeros(Int64, (num_nodes, num_nodes))

# ╔═╡ be49b7d3-4c70-4970-86e4-a4a6609e34b8
begin
start_nodes = [x[1] for x in members]
end_nodes = [x[2] for x in members]

for i = 1:num_nodes
	
	index = findall(x->x==i, start_nodes)
	
	g[i, end_nodes[index]] .= 1
	
end
end


# ╔═╡ 0d8d8f6a-5af8-467a-abfd-8bb6cc07621a
plot_truss(bridge, node_geometry, members, scale_x, scale_y)

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
# ╠═96db487a-421f-454a-9104-0f70e7ca46a6
# ╠═dbf428c5-f2e6-4801-a86d-1396f3ff258a
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
# ╠═08b08d8d-b0d7-4309-a190-69a9c9ce12da
# ╠═af383036-798f-4468-a189-bc7369d7ed17
# ╠═cd7e2974-4520-4a9e-9d42-e3b3e2d1f224
# ╠═8d3250ca-7b5b-4924-9634-51b1cee236fc
# ╠═408fbbf1-8758-4f03-97e4-a295f1f99634
# ╠═bcdafcf1-42a8-4cf7-b1b6-33dfad4bd862
# ╠═15fd59f7-2866-44f4-a5ea-3472b4a5fa05
# ╠═be49b7d3-4c70-4970-86e4-a4a6609e34b8
# ╠═0d8d8f6a-5af8-467a-abfd-8bb6cc07621a
