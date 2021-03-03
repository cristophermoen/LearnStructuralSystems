### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 47306786-6656-11eb-072f-ebc214b9b8e9
using LinearAlgebra

# ╔═╡ 8904b762-64ac-11eb-30d0-a7cc5cf93812
md""" EN.560.301 Structural Systems I

Spring 2021

"""

# ╔═╡ ebc991e0-6661-11eb-00c8-1747883147cb
md""" In this notebook we start playing with the stiffness method on the braced frame for Caribbean Town. Our goal is to calculate deflections and internal forces. We are idealizing the frame as a two bar truss system."""

# ╔═╡ 35818d78-6a41-11eb-23a3-71ebb88b7af5
md""" Units are in mm and N."""

# ╔═╡ 84a1e98c-6653-11eb-06b8-3daf7a963bba
md""" Define nodal geometry of our structural system.


"""

# ╔═╡ ea4e82a4-64ad-11eb-05f9-0dcb561d1268
#                 x  y
node_geometry = [0.0 0.0;
	             3100/tan(deg2rad(48)) 3100.0;
	             3100/tan(deg2rad(48))*2 0.0]  #mm
	             

# ╔═╡ 4d2e00b4-6654-11eb-1885-b17c28c298e9
md""" Define section properties."""

# ╔═╡ c7b26d94-6653-11eb-3bc4-0710205a79f0
section_properties = [6693.0]  #area A, mm^2

# ╔═╡ 2f049adc-6655-11eb-2f70-0fdab1a9b4ec
md""" Define material properties."""

# ╔═╡ a2bd25dc-6654-11eb-2628-b112c1bfa78d
material_properties = [199947.96] # elastic modulus E, N/mm^2

# ╔═╡ 33d2dcae-6655-11eb-28a3-b133c2b40e83
md""" Declare member connectivity, section properties, material properties."""

# ╔═╡ da5af3c0-6654-11eb-18f9-394a7785cf65
members = [(1, 2, 1, 1),
	        (2, 3, 1, 1)]   #node i, node j, A, E

# ╔═╡ 3f4b0a72-6655-11eb-2462-f3f2c16810c6
md""" Write a function that defines a local element stiffness matrix."""

# ╔═╡ 2cb42f42-6655-11eb-0925-455be7020f84
function k_local(E, A, L)
	
	k_local = E * A/ L * [1  0 -1 0
		                  0  0 0 0
		                  -1 0  1 0
	                      0  0  0 0]
	
	return k_local
	
end
	
	

# ╔═╡ daa7c6f6-6658-11eb-1ce2-093213e16a3e
md""" Calculate member lengths."""

# ╔═╡ 451deec8-6656-11eb-0f46-97946b2c930e
function get_member_lengths(members, node_geometry)

    L=zeros(length(members))

    for i=1:length(members)

        iNode=members[i][1]
        jNode=members[i][2]

        iNodePosition=node_geometry[iNode,:]
        jNodePosition=node_geometry[jNode,:]

        L[i] =norm(iNodePosition-jNodePosition)

    end

    return L

end

# ╔═╡ 46db0c78-6656-11eb-13dc-df7acdb32838
L = get_member_lengths(members, node_geometry)

# ╔═╡ 1634e6d6-6659-11eb-260f-5f0b748917e6
md""" Calculate element local stiffness matrices."""

# ╔═╡ 4664e8a4-6656-11eb-1e3b-6172ffe9a164
E = material_properties[1][1]

# ╔═╡ 3be15d94-6657-11eb-293c-556ecdb63a68
A = section_properties[1][1]

# ╔═╡ 4aa08e72-6657-11eb-3955-2fdd4feeb7f5

k_local_1 = k_local(E, A, L[1])

# ╔═╡ e43f31e0-6658-11eb-0fbe-6703d1f52c57
k_local_2 = k_local(E, A, L[2])

# ╔═╡ 5bc6e882-6659-11eb-2e4a-9b900693ebd1
md"""Calculate element stiffness matrices in global coordinates."""

# ╔═╡ 840a269e-6659-11eb-302a-330cf25cc147
md"""First find the member orientation angles relative to the horizon."""

# ╔═╡ c635393c-6659-11eb-2342-39c12aec0516
md"""Write a function to get this for us."""

# ╔═╡ 9357b8aa-6659-11eb-2825-cf34caa5a6f2
function get_member_orientations(members, node_geometry)

    θ=zeros(length(members))

    for i=1:length(members)

       θ[i] = atan((node_geometry[members[i][2],2]-node_geometry[members[i][1],2])/(node_geometry[members[i][2],1]-node_geometry[members[i][1],1]))
    end

    return θ

end

# ╔═╡ db9afec4-6659-11eb-1cb5-4b86b1601c97
θ = get_member_orientations(members, node_geometry)   #angle from the horizon in radians, positive angle is counterclockwise from the horizon

# ╔═╡ 5991f9b2-6660-11eb-02e2-a15c29c228bb
md""" Write a function that generates a vector rotation matrix that we will use to rotate the element displacement vector from local to global coordinates."""

# ╔═╡ 76026b0e-6660-11eb-2952-3f479ea21b4e
function rotation(θ)
	
	
	T = [cos(θ) sin(θ) 0 0
		 -sin(θ) cos(θ) 0 0
		 0 0 cos(θ) sin(θ)
		 0 0  -sin(θ) cos(θ)]
	
	return T
	
end
	
	

# ╔═╡ 4ce954b6-6661-11eb-05ef-87940936150d
md""" Get the rotation matrix for element 1.  See Sennett Section 2.2 for more details."""

# ╔═╡ 2f195012-6661-11eb-2573-adb484cce039
T_1 = rotation(θ[1])

# ╔═╡ 55ace586-6661-11eb-25f1-db3d5586af62
md"""Get rotation matrix for element 2."""

# ╔═╡ 5e316a56-6661-11eb-08c8-41df87131cad
T_2 = rotation(θ[2])

# ╔═╡ 6cd33c74-6661-11eb-3e0b-81843314eb9d
md"""Calculate the element 1 local stiffness matrix in global coordinates."""

# ╔═╡ 809d4e54-6661-11eb-2e59-65f130953ddc
k_global_1 = T_1' * k_local_1 * T_1

# ╔═╡ 95260cb0-6661-11eb-3568-d3b4b8513d5f
k_global_2 = T_2' * k_local_2 * T_2

# ╔═╡ a52e5996-6661-11eb-2ee5-8dd542b1a27c
md"""Now these element stiffness matrices are ready for assembly into the global stiffness matrix.  See Sennett Section 2.3 for some examples."""

# ╔═╡ 52f154ae-6a43-11eb-1917-51b32a568eee
md"""Initialize the global stiffness matrix."""

# ╔═╡ 459b5aac-6a43-11eb-06a8-eda347c520aa
K = zeros(Float64, (6, 6))

# ╔═╡ 8c772a1e-6a43-11eb-2f89-05c638646ea0
md"""Add element 1 stiffness matrix to the global stiffness matrix."""

# ╔═╡ 8173f692-6a43-11eb-0816-5d301ae62150
K[1:4, 1:4] .= k_global_1

# ╔═╡ b2a13b4e-6a43-11eb-3439-6118f910dd3b
K

# ╔═╡ da256a84-6a46-11eb-0d75-eda5c49c59d9
md"""Add element 2 stiffness matrix to the global stiffness matrix."""

# ╔═╡ b949b692-6a43-11eb-03b7-6f310efe1ed3
K[3:6, 3:6] .= K[3:6, 3:6] .+ k_global_2

# ╔═╡ f1346f34-6a43-11eb-3ddc-2fc85bfa0a1f
K

# ╔═╡ e56927be-6a46-11eb-28e4-474faf96958a
md"""Define the support conditions."""

# ╔═╡ f4c179be-6a43-11eb-30fe-a77ccd4ceef8
supports = [1, 1, 0, 0, 1, 1]

# ╔═╡ faded9b8-6a46-11eb-259b-c3eefa2b3511
md"""Define the free degrees of freedom."""

# ╔═╡ c6cb7f1a-6a45-11eb-1007-1143ee4e8819
free_dof = findall(x->x==0, supports)

# ╔═╡ 0a80b9ae-6a47-11eb-1d8d-53d3fe747312
md"""Initialize the part of the stiffness matrix that pertains just to the free degrees of freedom."""

# ╔═╡ ffe512f4-6a45-11eb-1a4f-25e2b178cdb4
Kff = zeros(Float64, (2,2))

# ╔═╡ 230f5186-6a47-11eb-081b-55a21b8ecdec
md"""Define this part of the stiffness matrix."""

# ╔═╡ eb3357ec-6a45-11eb-397a-2d91eda3ddef
Kff .= K[free_dof, free_dof]

# ╔═╡ 36bcb798-6a47-11eb-0c23-0d40b2e97bdc
md"""Define the external forces on the structural system."""

# ╔═╡ fb53ee98-6a45-11eb-1433-654b02ed169c
F = [0.0, 0.0, 500000.0, 0.0, 0.0, 0.0]

# ╔═╡ 42233ea4-6a47-11eb-1997-d3df5a0c7e50
md"""Partition the external force vector to only find the free degrees of freedom."""

# ╔═╡ 3f17da74-6a46-11eb-07bd-872fdea65998
Ff = zeros(Float64, 2)

# ╔═╡ 33675d6a-6a46-11eb-1828-27e3fae21556
Ff .= F[free_dof]

# ╔═╡ 5715746c-6a47-11eb-0edc-affdece3d65d
md"""Calculate the displacements at the free degrees of freedom."""

# ╔═╡ 2eee3ccc-6a46-11eb-3943-4bf3e186c5bf
uf = Kff^-1*Ff

# ╔═╡ e939c8b0-6bb5-11eb-1546-cbc8dc218f37
num_nodes = size(node_geometry)[1]

# ╔═╡ 223423ae-6bb6-11eb-2c9e-c5b76c736771
num_dof_per_node = 2  #for a two force element

# ╔═╡ 2c7dc964-6bb6-11eb-248c-db03b71db892
num_dof = num_nodes * num_dof_per_node

# ╔═╡ 4d02893e-6bb6-11eb-1b5a-6fa0f4f89e26
u = zeros(Float64, num_dof)

# ╔═╡ 56a675ee-6bb6-11eb-2f77-f17cce08e28d
u[free_dof] = uf

# ╔═╡ 627c84a6-6bb6-11eb-1616-07a1b31fae3a
u

# ╔═╡ 64733c1c-6a46-11eb-31d6-3bc96f3235a5
fixed_dof = findall(x->x==1, supports)

# ╔═╡ d1fecb82-6bb5-11eb-021f-c1af7bb6fa78
s = fixed_dof

# ╔═╡ d824da10-6bb5-11eb-1a41-1ff2741aec0e
f = free_dof

# ╔═╡ 891e52a6-6bb6-11eb-25b8-8bb84fb6563e
Ksf = zeros(Float64, (length(s), length(f)))

# ╔═╡ de1998ac-6bb5-11eb-207d-9d9d31e7107b
Ksf .= K[s, f]

# ╔═╡ bdeca2ec-6fc3-11eb-2eda-ada0ae2bd56c
K

# ╔═╡ ac4608c8-6bb6-11eb-1c46-fba5d87d7974
Fs = Ksf * uf

# ╔═╡ 17ec726c-6bd3-11eb-0454-f9f79fcafd13
md"""Define global displacement vector for element 1."""

# ╔═╡ 831d4adc-6bb7-11eb-069a-2f51cbf327cf
u_global_e1 = u[1:4]

# ╔═╡ 9eb8df04-6bb7-11eb-2a33-6f9324c55428
u_local_e1 = T_1 * u_global_e1

# ╔═╡ 06c94c9a-6bb9-11eb-2d49-010baf4e0622
f_local_e1 = k_local_1 * u_local_e1

# ╔═╡ Cell order:
# ╟─8904b762-64ac-11eb-30d0-a7cc5cf93812
# ╠═ebc991e0-6661-11eb-00c8-1747883147cb
# ╟─35818d78-6a41-11eb-23a3-71ebb88b7af5
# ╟─84a1e98c-6653-11eb-06b8-3daf7a963bba
# ╠═ea4e82a4-64ad-11eb-05f9-0dcb561d1268
# ╟─4d2e00b4-6654-11eb-1885-b17c28c298e9
# ╠═c7b26d94-6653-11eb-3bc4-0710205a79f0
# ╟─2f049adc-6655-11eb-2f70-0fdab1a9b4ec
# ╠═a2bd25dc-6654-11eb-2628-b112c1bfa78d
# ╟─33d2dcae-6655-11eb-28a3-b133c2b40e83
# ╠═da5af3c0-6654-11eb-18f9-394a7785cf65
# ╟─3f4b0a72-6655-11eb-2462-f3f2c16810c6
# ╠═2cb42f42-6655-11eb-0925-455be7020f84
# ╟─daa7c6f6-6658-11eb-1ce2-093213e16a3e
# ╠═451deec8-6656-11eb-0f46-97946b2c930e
# ╠═47306786-6656-11eb-072f-ebc214b9b8e9
# ╠═46db0c78-6656-11eb-13dc-df7acdb32838
# ╟─1634e6d6-6659-11eb-260f-5f0b748917e6
# ╠═4664e8a4-6656-11eb-1e3b-6172ffe9a164
# ╠═3be15d94-6657-11eb-293c-556ecdb63a68
# ╠═4aa08e72-6657-11eb-3955-2fdd4feeb7f5
# ╠═e43f31e0-6658-11eb-0fbe-6703d1f52c57
# ╟─5bc6e882-6659-11eb-2e4a-9b900693ebd1
# ╟─840a269e-6659-11eb-302a-330cf25cc147
# ╟─c635393c-6659-11eb-2342-39c12aec0516
# ╠═9357b8aa-6659-11eb-2825-cf34caa5a6f2
# ╠═db9afec4-6659-11eb-1cb5-4b86b1601c97
# ╟─5991f9b2-6660-11eb-02e2-a15c29c228bb
# ╠═76026b0e-6660-11eb-2952-3f479ea21b4e
# ╟─4ce954b6-6661-11eb-05ef-87940936150d
# ╠═2f195012-6661-11eb-2573-adb484cce039
# ╟─55ace586-6661-11eb-25f1-db3d5586af62
# ╠═5e316a56-6661-11eb-08c8-41df87131cad
# ╟─6cd33c74-6661-11eb-3e0b-81843314eb9d
# ╠═809d4e54-6661-11eb-2e59-65f130953ddc
# ╠═95260cb0-6661-11eb-3568-d3b4b8513d5f
# ╠═a52e5996-6661-11eb-2ee5-8dd542b1a27c
# ╠═52f154ae-6a43-11eb-1917-51b32a568eee
# ╠═459b5aac-6a43-11eb-06a8-eda347c520aa
# ╟─8c772a1e-6a43-11eb-2f89-05c638646ea0
# ╠═8173f692-6a43-11eb-0816-5d301ae62150
# ╠═b2a13b4e-6a43-11eb-3439-6118f910dd3b
# ╠═da256a84-6a46-11eb-0d75-eda5c49c59d9
# ╠═b949b692-6a43-11eb-03b7-6f310efe1ed3
# ╠═f1346f34-6a43-11eb-3ddc-2fc85bfa0a1f
# ╠═e56927be-6a46-11eb-28e4-474faf96958a
# ╠═f4c179be-6a43-11eb-30fe-a77ccd4ceef8
# ╠═faded9b8-6a46-11eb-259b-c3eefa2b3511
# ╠═c6cb7f1a-6a45-11eb-1007-1143ee4e8819
# ╠═0a80b9ae-6a47-11eb-1d8d-53d3fe747312
# ╠═ffe512f4-6a45-11eb-1a4f-25e2b178cdb4
# ╠═230f5186-6a47-11eb-081b-55a21b8ecdec
# ╠═eb3357ec-6a45-11eb-397a-2d91eda3ddef
# ╠═36bcb798-6a47-11eb-0c23-0d40b2e97bdc
# ╠═fb53ee98-6a45-11eb-1433-654b02ed169c
# ╠═42233ea4-6a47-11eb-1997-d3df5a0c7e50
# ╠═3f17da74-6a46-11eb-07bd-872fdea65998
# ╠═33675d6a-6a46-11eb-1828-27e3fae21556
# ╠═5715746c-6a47-11eb-0edc-affdece3d65d
# ╠═2eee3ccc-6a46-11eb-3943-4bf3e186c5bf
# ╠═e939c8b0-6bb5-11eb-1546-cbc8dc218f37
# ╠═223423ae-6bb6-11eb-2c9e-c5b76c736771
# ╠═2c7dc964-6bb6-11eb-248c-db03b71db892
# ╠═4d02893e-6bb6-11eb-1b5a-6fa0f4f89e26
# ╠═56a675ee-6bb6-11eb-2f77-f17cce08e28d
# ╠═627c84a6-6bb6-11eb-1616-07a1b31fae3a
# ╠═64733c1c-6a46-11eb-31d6-3bc96f3235a5
# ╠═d1fecb82-6bb5-11eb-021f-c1af7bb6fa78
# ╠═d824da10-6bb5-11eb-1a41-1ff2741aec0e
# ╠═891e52a6-6bb6-11eb-25b8-8bb84fb6563e
# ╠═de1998ac-6bb5-11eb-207d-9d9d31e7107b
# ╠═bdeca2ec-6fc3-11eb-2eda-ada0ae2bd56c
# ╠═ac4608c8-6bb6-11eb-1c46-fba5d87d7974
# ╟─17ec726c-6bd3-11eb-0454-f9f79fcafd13
# ╠═831d4adc-6bb7-11eb-069a-2f51cbf327cf
# ╠═9eb8df04-6bb7-11eb-2a33-6f9324c55428
# ╠═06c94c9a-6bb9-11eb-2d49-010baf4e0622
