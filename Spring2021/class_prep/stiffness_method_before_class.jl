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
md""" In this lecture we start playing with the stiffness method on the braced frame for Caribbean Town.  We are idealizing the frame as a two bar truss."""

# ╔═╡ fcd98e56-6a22-11eb-0534-df45ca65520b
md"""Units are in N and mm."""

# ╔═╡ 84a1e98c-6653-11eb-06b8-3daf7a963bba
md""" Define nodal geometry of our structural system.


"""

# ╔═╡ ea4e82a4-64ad-11eb-05f9-0dcb561d1268
node_geometry = [0.0 0.0;
	             3100/tan(deg2rad(48)) 3100.0;
	             3100/tan(deg2rad(48))*2 0.0]  #mm
	             

# ╔═╡ 4d2e00b4-6654-11eb-1885-b17c28c298e9
md""" Define section properties."""

# ╔═╡ c7b26d94-6653-11eb-3bc4-0710205a79f0
section_properties = [6693.0]  #mm^2

# ╔═╡ 2f049adc-6655-11eb-2f70-0fdab1a9b4ec
md""" Define material properties."""

# ╔═╡ a2bd25dc-6654-11eb-2628-b112c1bfa78d
material_properties = [199947.96] # elastic modulus E, N/mm^2

# ╔═╡ 33d2dcae-6655-11eb-28a3-b133c2b40e83
md""" Declare member connectivity, section properties, material properties."""

# ╔═╡ da5af3c0-6654-11eb-18f9-394a7785cf65
members = [(1, 2, 1, 1),
	        (2, 3, 1, 1)]   #node i, node j, A, E

# ╔═╡ 8d8db00a-6a23-11eb-3985-2372a841c5a5
md""" Define the support conditions for this structural system."""

# ╔═╡ 9f8ce41a-6a23-11eb-020f-9d4158aa7522
supports = [1 1 0 0 1 1]   #0 is free, 1 is fixed for global dof

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

# ╔═╡ 3d0c896c-6a23-11eb-0d70-f9a69a1d4654
md""" Now get the element 2 local stiffness matrix in global coordinates."""

# ╔═╡ 95260cb0-6661-11eb-3568-d3b4b8513d5f
k_global_2 = T_2' * k_local_2 * T_2

# ╔═╡ a52e5996-6661-11eb-2ee5-8dd542b1a27c
md"""Now these element stiffness matrices are ready for assembly into the global stiffness matrix.  See Sennett Section 2.3 for some examples."""

# ╔═╡ 0cca81e2-6a27-11eb-2700-b933a59d9d8f
md"""First let's initialize the global stiffness matrix. """

# ╔═╡ 39eae748-6a27-11eb-1de7-a93fe77e0c5b
md"""We need to figure out how big the stiffness matrix should be.

# ╔═╡ 09856d58-6a27-11eb-3c8a-93199ab8eb56
K = zeros(Float64, (

# ╔═╡ da39d0e6-6a23-11eb-2038-23c5eeda62dc
md"""Let's get organized by defining the mapping from the local dof to the global dof for Element 1. """

# ╔═╡ 515a08ca-6a25-11eb-209a-e3a3283bc286
element_index = 1

# ╔═╡ e9a0d6f0-6a24-11eb-29be-e151e2fc0799
node_i = members[1][element_index]

# ╔═╡ 2273e530-6a25-11eb-0bcf-6167915abdc6
node_j = members[1][element_index]

# ╔═╡ b377a35a-6a25-11eb-217e-b1369f9bbcac
num_dof_per_node = 2

# ╔═╡ cc03b782-6a26-11eb-0547-7105f6e053c8
num_nodes_per_element = 2

# ╔═╡ 346a1a98-6a25-11eb-371a-6f4d1ff07117
global_dof_node_i = [1 2] .+ (element_index - 1) * num_dof_per_node * num_nodes_per_element

# ╔═╡ 792f9478-6a25-11eb-0910-651cae3ead6c
global_dof_node_j = [3 4] .+ (element_index - 1) * num_dof_per_node * num_nodes_per_element

# ╔═╡ Cell order:
# ╟─8904b762-64ac-11eb-30d0-a7cc5cf93812
# ╟─ebc991e0-6661-11eb-00c8-1747883147cb
# ╟─fcd98e56-6a22-11eb-0534-df45ca65520b
# ╟─84a1e98c-6653-11eb-06b8-3daf7a963bba
# ╠═ea4e82a4-64ad-11eb-05f9-0dcb561d1268
# ╟─4d2e00b4-6654-11eb-1885-b17c28c298e9
# ╠═c7b26d94-6653-11eb-3bc4-0710205a79f0
# ╟─2f049adc-6655-11eb-2f70-0fdab1a9b4ec
# ╠═a2bd25dc-6654-11eb-2628-b112c1bfa78d
# ╟─33d2dcae-6655-11eb-28a3-b133c2b40e83
# ╠═da5af3c0-6654-11eb-18f9-394a7785cf65
# ╠═8d8db00a-6a23-11eb-3985-2372a841c5a5
# ╠═9f8ce41a-6a23-11eb-020f-9d4158aa7522
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
# ╟─9357b8aa-6659-11eb-2825-cf34caa5a6f2
# ╠═db9afec4-6659-11eb-1cb5-4b86b1601c97
# ╟─5991f9b2-6660-11eb-02e2-a15c29c228bb
# ╠═76026b0e-6660-11eb-2952-3f479ea21b4e
# ╟─4ce954b6-6661-11eb-05ef-87940936150d
# ╠═2f195012-6661-11eb-2573-adb484cce039
# ╟─55ace586-6661-11eb-25f1-db3d5586af62
# ╠═5e316a56-6661-11eb-08c8-41df87131cad
# ╟─6cd33c74-6661-11eb-3e0b-81843314eb9d
# ╠═809d4e54-6661-11eb-2e59-65f130953ddc
# ╟─3d0c896c-6a23-11eb-0d70-f9a69a1d4654
# ╠═95260cb0-6661-11eb-3568-d3b4b8513d5f
# ╠═a52e5996-6661-11eb-2ee5-8dd542b1a27c
# ╠═0cca81e2-6a27-11eb-2700-b933a59d9d8f
# ╠═39eae748-6a27-11eb-1de7-a93fe77e0c5b
# ╠═09856d58-6a27-11eb-3c8a-93199ab8eb56
# ╠═da39d0e6-6a23-11eb-2038-23c5eeda62dc
# ╠═515a08ca-6a25-11eb-209a-e3a3283bc286
# ╠═e9a0d6f0-6a24-11eb-29be-e151e2fc0799
# ╠═2273e530-6a25-11eb-0bcf-6167915abdc6
# ╠═b377a35a-6a25-11eb-217e-b1369f9bbcac
# ╠═cc03b782-6a26-11eb-0547-7105f6e053c8
# ╠═346a1a98-6a25-11eb-371a-6f4d1ff07117
# ╠═792f9478-6a25-11eb-0910-651cae3ead6c
