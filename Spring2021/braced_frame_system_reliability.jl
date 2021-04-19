### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ c650bc45-c97a-49e4-b300-77b664e36632
using StructuresKit

# ╔═╡ 98dcd238-8040-43bf-83fc-2304479a2129
using Distributions

# ╔═╡ c6036863-9cf9-43cf-8f06-84750aa02e3c
using Plots

# ╔═╡ 0a9cc3e9-7334-4dfa-8ec9-e7d16cb9bca9
using Interpolations

# ╔═╡ 8904b762-64ac-11eb-30d0-a7cc5cf93812
md""" EN.560.301 Structural Systems I

Spring 2021

"""

# ╔═╡ 408727a9-55d3-48fb-86fc-5ed2392087eb
md""" # Calculate the probability of failure of a Caribbean Town 1st floor ordinary concentrically brace frame."""

# ╔═╡ 35818d78-6a41-11eb-23a3-71ebb88b7af5
md""" Units are in mm and N."""

# ╔═╡ 84a1e98c-6653-11eb-06b8-3daf7a963bba
md""" ## Define the structural model.


"""

# ╔═╡ ea4e82a4-64ad-11eb-05f9-0dcb561d1268
#                 x  y
node_geometry = [0.0 0.0;
	             3100/tan(deg2rad(48)) 3100.0;
	             3100/tan(deg2rad(48))*2 0.0]  #mm
	             

# ╔═╡ c7b26d94-6653-11eb-3bc4-0710205a79f0
section_properties = [6693.0]  #area A, mm^2

# ╔═╡ a2bd25dc-6654-11eb-2628-b112c1bfa78d
material_properties = [199947.96] # elastic modulus E, N/mm^2

# ╔═╡ 01b500e5-fe25-465c-a160-5c6b7c033aa4
members = [(1, 2, 1, 1),
	        (2, 3, 1, 1)]   #node i, node j, A, E

# ╔═╡ f4c179be-6a43-11eb-30fe-a77ccd4ceef8
supports = [1, 1, 0, 0, 1, 1]

# ╔═╡ 38d7f2d4-1b32-4f38-8b8b-512bba9a5cad
md"""Define a unit lateral load on the frame."""

# ╔═╡ 532dca2f-0e71-4ac3-821b-8e498ad75ef0
external_forces = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

# ╔═╡ fc0b735c-6685-48a5-ae78-5629524b1520
begin

	frame = Truss.define(members, section_properties, material_properties, node_geometry, supports, external_forces)
	
	frame = Truss.solve(frame)
	
end

# ╔═╡ 8b507824-19b2-49b2-86e0-adbabfc555f0
Truss.show(frame, 1000000000.0, 10000000000.0)

# ╔═╡ a3ef3d72-a600-48c8-b461-295c4d4c9776
md""" ## Define nominal capacities for each member strength limit state."""

# ╔═╡ b3df1c6f-bf1f-457a-9bb0-ade55a0c419b
Pn_tensile_yielding = 2307 #kN

# ╔═╡ afeb659b-8010-4ad9-8880-e753eea024bc
Pn_tensile_rupture = 1879 #kN

# ╔═╡ 7070216e-8aa8-49f3-8b73-ce508f70de38
Pn_compression = 1465 #kN

# ╔═╡ 496c70a6-2cbb-4a82-86e7-5961900d17cd
md""" ## Define strength limit state statistics and distributions using NBS 577 Table 3.3."""

# ╔═╡ c0ef507e-bd1b-4dd8-bc9b-d2221f836b91
μ_R_tensile_yielding = Pn_tensile_yielding * 1.05

# ╔═╡ 3b94f125-a9b5-40b1-b2f4-0e45d9ade87e
σ_R_tensile_yielding = 0.11 * μ_R_tensile_yielding

# ╔═╡ 10ab60d8-b7a5-4ee8-9456-7546a077556b
prob_dist_tensile_yielding = Normal(μ_R_tensile_yielding, σ_R_tensile_yielding)

# ╔═╡ 8386c690-82b5-4de5-aead-7b8c46046426
μ_R_tensile_rupture = Pn_tensile_rupture * 1.10

# ╔═╡ 32f6379f-4f37-4376-ab16-c693d04384a2
σ_R_tensile_rupture = 0.11 * μ_R_tensile_rupture

# ╔═╡ 3ff216d1-08c1-42fe-94ca-d49a519389de
prob_dist_tensile_rupture = Normal(μ_R_tensile_rupture, σ_R_tensile_rupture)

# ╔═╡ 58177b62-a677-420f-9055-78f2f0be58b9
μ_R_compression = Pn_compression * 1.08

# ╔═╡ f5b86970-b43d-4bfe-8f92-85eec8ea0642
σ_R_compression = 0.14 * μ_R_compression

# ╔═╡ 94bb064b-f4b4-4651-a488-5d039db07702
prob_dist_compression = Normal(μ_R_compression, σ_R_compression)

# ╔═╡ 359b7deb-4329-4aae-bd63-ae1b5c38a4a8
md""" ## Define cumulative distribution function of the earthquake loading."""

# ╔═╡ dbe43db0-1280-4b3c-b3ef-b16d0b16132b
function TypeII(a, a10)
	
	F = exp(-(a/(0.38 * a10))^(-2.3))
	
end

# ╔═╡ 235ca216-5e63-4fe1-8ccb-e3b1fed25402
a10 = 0.4 #g

# ╔═╡ 3e60400b-86c8-4bdd-b399-cdfdb8f367f4
a=0.001:0.001:2.0

# ╔═╡ ed73aa58-dfc7-4201-9049-94ac2d3ab17d
begin
	num_a = length(a)
	F_eq = zeros(Float64, num_a)
	for i=1:num_a
		F_eq[i] = TypeII(a[i], a10)
	end
end

# ╔═╡ 10c92939-51a2-4efb-abb5-9da87a901ea0
plot(a, F_eq, legend = false)

# ╔═╡ 22ebbe15-50b6-4c61-a73d-098c566f2ff3
a_interp = LinearInterpolation(F_eq, a, extrapolation_bc=Line())

# ╔═╡ f217dac9-ee45-4441-a203-37582253d437
md""" ## Calculate random instances."""

# ╔═╡ 7c6ba7d8-3c0a-490a-a9b0-f69937a4e2ab
instances = 10000000

# ╔═╡ eaac73b7-dc06-4c2d-81cf-03f69ed8f2cc
R_tensile_yielding = rand(prob_dist_tensile_yielding, instances)

# ╔═╡ 993fbf82-c690-400f-9b9a-55231357a683
histogram(R_tensile_yielding)

# ╔═╡ f671c7de-521c-44e2-9a4b-4a723fa291f2
R_tensile_rupture = rand(prob_dist_tensile_rupture, instances)

# ╔═╡ bf669f56-6a1c-4b47-9425-245ef8f32abe
R_compression = rand(prob_dist_compression, instances)

# ╔═╡ f5add394-2437-4f8d-b546-fcf2a8214756
md"""Take earthquake loading random samples."""

# ╔═╡ 98fa213d-628e-4137-a420-c7329fa34cd6
F_eq_samples = rand(instances)

# ╔═╡ 150b7d96-bcaa-42f1-aa4b-86a09629ee2a
a_eq_samples = a_interp(F_eq_samples)

# ╔═╡ 19b3d84d-ebf8-4c22-bdfd-ad3b3a049440
histogram(a_eq_samples)

# ╔═╡ 3d32199b-8aa3-4a04-939b-e000b75428f3
seismic_mass = 3063000/4 #kg, per brace spine

# ╔═╡ 10fdd21a-ceee-4c30-84ba-6e80d5c64af2
R = 3.25

# ╔═╡ 54702fd2-1653-4f77-bf54-21979fed5c15
V = (2/3 .* a_eq_samples .* seismic_mass .* 9.81) ./ R #N

# ╔═╡ b1d46572-2872-483c-85c3-324be790092c
Q = frame.f_element[1][3] .* V /1000 #kN

# ╔═╡ 61400532-3f62-42cc-acfa-a9f6a1a100ff
histogram(Q)

# ╔═╡ 7d45de4c-fb40-42c3-bb0f-1338f63eb671
md""" ## Check demand-to-capacity ratios for each limit state."""

# ╔═╡ 8c1ddca8-6f55-463d-b752-53dc1cacec70
DC_tensile_yielding = Q ./ R_tensile_yielding

# ╔═╡ 3a01038d-493b-4c59-ab21-b067de426958
DC_tensile_rupture = Q ./ R_tensile_rupture

# ╔═╡ 70428671-029e-4899-bbc4-d544eacaecec
DC_compression = Q ./ R_compression

# ╔═╡ 20a8d148-d433-46d4-94c7-1368dfa40a50
DC_min = maximum([DC_tensile_yielding DC_tensile_rupture DC_compression], dims=2)

# ╔═╡ 42f8e7c6-0ee5-4536-8deb-ae82e3210681
md""" ## Calculate system probability of failure."""

# ╔═╡ 41cd0099-2390-48a1-9876-36c78936c199
failures = findall(x->x>=1.0, DC_min)

# ╔═╡ b600eb9c-10fb-4548-9604-0ebaf735ddd5
num_failures = length(failures)

# ╔═╡ 7328403a-dbb7-4924-8f78-740c4c0ba8a4
system_probability_of_failure = num_failures/instances

# ╔═╡ 4395361d-37dc-4a77-8009-c0d1af0ef4c1
num_tension_yielding_failures = length(findall(x->x>=1.0, DC_tensile_yielding))

# ╔═╡ e0324e49-af30-4045-8ee7-fd4eea9d2f09
num_tension_rupture_failures = length(findall(x->x>=1.0, DC_tensile_rupture))

# ╔═╡ a726dbeb-62f4-461e-9894-bdf2f2fa092e
num_compression_failures = length(findall(x->x>=1.0, DC_compression))

# ╔═╡ Cell order:
# ╟─8904b762-64ac-11eb-30d0-a7cc5cf93812
# ╠═408727a9-55d3-48fb-86fc-5ed2392087eb
# ╟─35818d78-6a41-11eb-23a3-71ebb88b7af5
# ╠═84a1e98c-6653-11eb-06b8-3daf7a963bba
# ╠═ea4e82a4-64ad-11eb-05f9-0dcb561d1268
# ╠═c7b26d94-6653-11eb-3bc4-0710205a79f0
# ╠═a2bd25dc-6654-11eb-2628-b112c1bfa78d
# ╠═01b500e5-fe25-465c-a160-5c6b7c033aa4
# ╠═f4c179be-6a43-11eb-30fe-a77ccd4ceef8
# ╠═38d7f2d4-1b32-4f38-8b8b-512bba9a5cad
# ╠═532dca2f-0e71-4ac3-821b-8e498ad75ef0
# ╠═c650bc45-c97a-49e4-b300-77b664e36632
# ╠═fc0b735c-6685-48a5-ae78-5629524b1520
# ╠═8b507824-19b2-49b2-86e0-adbabfc555f0
# ╠═a3ef3d72-a600-48c8-b461-295c4d4c9776
# ╠═b3df1c6f-bf1f-457a-9bb0-ade55a0c419b
# ╠═afeb659b-8010-4ad9-8880-e753eea024bc
# ╠═7070216e-8aa8-49f3-8b73-ce508f70de38
# ╠═496c70a6-2cbb-4a82-86e7-5961900d17cd
# ╠═c0ef507e-bd1b-4dd8-bc9b-d2221f836b91
# ╠═3b94f125-a9b5-40b1-b2f4-0e45d9ade87e
# ╠═98dcd238-8040-43bf-83fc-2304479a2129
# ╠═10ab60d8-b7a5-4ee8-9456-7546a077556b
# ╠═8386c690-82b5-4de5-aead-7b8c46046426
# ╠═32f6379f-4f37-4376-ab16-c693d04384a2
# ╠═3ff216d1-08c1-42fe-94ca-d49a519389de
# ╠═58177b62-a677-420f-9055-78f2f0be58b9
# ╠═f5b86970-b43d-4bfe-8f92-85eec8ea0642
# ╠═94bb064b-f4b4-4651-a488-5d039db07702
# ╠═359b7deb-4329-4aae-bd63-ae1b5c38a4a8
# ╠═dbe43db0-1280-4b3c-b3ef-b16d0b16132b
# ╠═235ca216-5e63-4fe1-8ccb-e3b1fed25402
# ╠═3e60400b-86c8-4bdd-b399-cdfdb8f367f4
# ╠═ed73aa58-dfc7-4201-9049-94ac2d3ab17d
# ╠═c6036863-9cf9-43cf-8f06-84750aa02e3c
# ╠═10c92939-51a2-4efb-abb5-9da87a901ea0
# ╠═0a9cc3e9-7334-4dfa-8ec9-e7d16cb9bca9
# ╠═22ebbe15-50b6-4c61-a73d-098c566f2ff3
# ╠═f217dac9-ee45-4441-a203-37582253d437
# ╠═7c6ba7d8-3c0a-490a-a9b0-f69937a4e2ab
# ╠═eaac73b7-dc06-4c2d-81cf-03f69ed8f2cc
# ╠═993fbf82-c690-400f-9b9a-55231357a683
# ╠═f671c7de-521c-44e2-9a4b-4a723fa291f2
# ╠═bf669f56-6a1c-4b47-9425-245ef8f32abe
# ╠═f5add394-2437-4f8d-b546-fcf2a8214756
# ╠═98fa213d-628e-4137-a420-c7329fa34cd6
# ╠═150b7d96-bcaa-42f1-aa4b-86a09629ee2a
# ╠═19b3d84d-ebf8-4c22-bdfd-ad3b3a049440
# ╠═3d32199b-8aa3-4a04-939b-e000b75428f3
# ╠═10fdd21a-ceee-4c30-84ba-6e80d5c64af2
# ╠═54702fd2-1653-4f77-bf54-21979fed5c15
# ╠═b1d46572-2872-483c-85c3-324be790092c
# ╠═61400532-3f62-42cc-acfa-a9f6a1a100ff
# ╠═7d45de4c-fb40-42c3-bb0f-1338f63eb671
# ╠═8c1ddca8-6f55-463d-b752-53dc1cacec70
# ╠═3a01038d-493b-4c59-ab21-b067de426958
# ╠═70428671-029e-4899-bbc4-d544eacaecec
# ╠═20a8d148-d433-46d4-94c7-1368dfa40a50
# ╠═42f8e7c6-0ee5-4536-8deb-ae82e3210681
# ╠═41cd0099-2390-48a1-9876-36c78936c199
# ╠═b600eb9c-10fb-4548-9604-0ebaf735ddd5
# ╠═7328403a-dbb7-4924-8f78-740c4c0ba8a4
# ╠═4395361d-37dc-4a77-8009-c0d1af0ef4c1
# ╠═e0324e49-af30-4045-8ee7-fd4eea9d2f09
# ╠═a726dbeb-62f4-461e-9894-bdf2f2fa092e
