### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 5ea7d0a8-7c3f-11eb-0d4a-35ba92e372ad
using Plots

# ╔═╡ 8904b762-64ac-11eb-30d0-a7cc5cf93812
md""" EN.560.301 Structural Systems I

Spring 2021

"""

# ╔═╡ 35818d78-6a41-11eb-23a3-71ebb88b7af5
md""" Units are in mm and N."""

# ╔═╡ b35e496e-7c3b-11eb-0351-1b0aacd9223d
md""" Perform a structural analysis of the Eduardo Torroja CASA aircraft manufacturing building.  Focus on a typical parabolic arch that support the doubly-curved thin concrete shell roof.   The arch has concrete top and bottom chords and a steel diagonal truss network."""

# ╔═╡ 20c6a79e-801f-11eb-173a-01464071310d
md""" ## Define the bottom chord geometry."""

# ╔═╡ 432b5b4a-801f-11eb-30b3-e79187e70ed4
md""" We need the length of the arch along the parabola."""

# ╔═╡ 5ef33884-801f-11eb-0e53-5910f6b83bbd
md""" Define some points along the x-direction (horizontal)."""

# ╔═╡ efebebb4-7c3b-11eb-1250-e362a93d4f7a
x_bottom_chord = 0:200:73000 #mm

# ╔═╡ 89b4975c-801f-11eb-028b-99897b274f19
md"""Define the centerline height of the bottom chord."""

# ╔═╡ 8c645c66-7c3d-11eb-1955-f7e59714787b
h_bottom_chord = 21000 + 350/2 #mm

# ╔═╡ d95d4c22-801f-11eb-14a2-833e11a65f3d
md"""Define the centerline bottom chord arch span."""

# ╔═╡ a1f706d2-7c3d-11eb-3fbe-e3872d192c05
L_bottom_chord = 73000

# ╔═╡ 0f1f308c-8020-11eb-3168-01946ee4b4d1
md"""Define the coefficients of a parabola for the bottom arch, $y=ax^2+bx+c$."""

# ╔═╡ 6bcda4c6-7c3d-11eb-29b3-bbdb4044d641
b_bottom_chord = 4 * h_bottom_chord/ L_bottom_chord

# ╔═╡ a4d50c28-7c3d-11eb-0cfb-27d78f900844
a_bottom_chord = -b_bottom_chord/L_bottom_chord

# ╔═╡ 0545ce32-7c46-11eb-3ff4-bf024fb73811
c_bottom_chord= 0.0

# ╔═╡ 7b570fb2-7c3e-11eb-09f8-5df49e7f5a9b
y_bottom_chord = a_bottom_chord .* x_bottom_chord.^2 .+ b_bottom_chord .* x_bottom_chord .+ c_bottom_chord

# ╔═╡ 5a2eeb1a-7c3f-11eb-1e9e-8db2708de320
plot(x_bottom_chord, y_bottom_chord, markershape = :o, size=(600,600), legend = false, xlims = (0, L_bottom_chord), ylims = (0,L_bottom_chord))

# ╔═╡ ac2c05c0-8021-11eb-054a-d5a01e53b98e


# ╔═╡ 7d8b9fe6-7c40-11eb-1603-d7f60457bd55
# dx = diff(x)

# ╔═╡ 81664378-8021-11eb-3f3b-195bd4e83371


# ╔═╡ 8343b30c-7c42-11eb-3ae3-1ffe6a12f9d9
# dy = diff(y)

# ╔═╡ 8ce16c56-7c42-11eb-253d-9d1d6809379d
# ds = sqrt.(dx.^2 + dy.^2)

# ╔═╡ adb3ebfc-7c42-11eb-10b6-898daac90b45
# s = sum(ds)

# ╔═╡ b896fb5e-7c42-11eb-214b-9d348d1c2dba
# num_arch_segments = 36

# ╔═╡ 730b3aa2-7c40-11eb-249b-655ca7caa1b9
# length_arch_segment = s / num_arch_segments

# ╔═╡ 3d9d470e-7c43-11eb-39e8-dd88a86846fd
# arch_layout = [length_arch_segment/2 ; ones(Float64, num_arch_segments - 1) .* length_arch_segment; length_arch_segment/2]

# ╔═╡ 315336f4-7c5a-11eb-0a00-83e875d3c31c


# ╔═╡ 345ea3d0-7c44-11eb-0068-f7f694943b54
function arch_angle(a, b, x)
	
	θ = atan(2 * a * x + b)
	
	return θ
	
end

# ╔═╡ c220d69e-7c43-11eb-05dd-21b150b8d436
# num_bottom_chord_nodes = length(arch_layout) + 1

# ╔═╡ 230d0878-7c45-11eb-0f2c-976edf65ab04
# num_bottom_chord_segments = length(arch_layout)

# ╔═╡ 2cb6793e-7c45-11eb-1b36-1b179a6b1168
# dx_segment = zeros(Float64, num_bottom_chord_segments)

# ╔═╡ 5d355172-7c40-11eb-1334-9fb53014f9a7
function calculate_arch_dx_segment(num_bottom_chord_segments, arch_layout)
	
	x_index = 0.0

	for i = 1:num_bottom_chord_segments
	
		arch_segment = arch_layout[i]
		
		θ = arch_angle(a, b, x_index)
	
		dx_segment[i] = arch_segment * cos(θ)
		
		x_index = x_index + dx_segment[i]
		
	end
	
	return dx_segment
	
end
	

# ╔═╡ ba6f73fe-7c4d-11eb-3a8a-912e71e751ce
# dx_arch = calculate_arch_dx_segment(num_bottom_chord_segments, arch_layout)

# ╔═╡ f8eea1ce-7c4d-11eb-2a18-21ae5d63457e
# x_bottom_chord = [0.0; cumsum(dx_arch)]

# ╔═╡ 83aaf5a6-7c5a-11eb-24ed-379e47866ce7
# y_bottom_chord = a .* x_bottom_chord.^2 .+ b .* x_bottom_chord .+ c

# ╔═╡ 9b21b178-7c5a-11eb-1ac7-ffd7ea79a5a1
# plot(x_bottom_chord, y_bottom_chord, markershape = :o, size=(600,600), legend = false, xlims = (0, L), ylims = (0,L))

# ╔═╡ b4cd33e0-7c5a-11eb-1e70-ef71af2e3e2f


# ╔═╡ f087d26e-7c4d-11eb-20e0-ff49a8cb1514


# ╔═╡ df8801be-7c4d-11eb-3c88-4705e4dff1c2


# ╔═╡ d4fd6f8e-7c4d-11eb-0512-49222239c4ce


# ╔═╡ c28accfc-7c4d-11eb-0eb1-bffa0f2b0887


# ╔═╡ 81599eca-7c4d-11eb-0af0-b3b02e43f112


# ╔═╡ 76d5d6bc-7c4d-11eb-16ea-cf98397e0f90


# ╔═╡ 5aae7ad4-7c4d-11eb-1d94-a5a3ada7ecd6


# ╔═╡ 38700098-7c46-11eb-034c-d787179ef393


# ╔═╡ fc8fff88-7c45-11eb-11c4-c95c7d35d388


# ╔═╡ 57101fac-7c40-11eb-3af8-557ad12ca677


# ╔═╡ 51d9b4bc-7c40-11eb-1412-5fae0ff8268d


# ╔═╡ e156ecdc-7c3f-11eb-00c6-c35ad2ca2468


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

# ╔═╡ Cell order:
# ╠═8904b762-64ac-11eb-30d0-a7cc5cf93812
# ╟─35818d78-6a41-11eb-23a3-71ebb88b7af5
# ╠═b35e496e-7c3b-11eb-0351-1b0aacd9223d
# ╠═20c6a79e-801f-11eb-173a-01464071310d
# ╠═432b5b4a-801f-11eb-30b3-e79187e70ed4
# ╠═5ef33884-801f-11eb-0e53-5910f6b83bbd
# ╠═efebebb4-7c3b-11eb-1250-e362a93d4f7a
# ╠═89b4975c-801f-11eb-028b-99897b274f19
# ╠═8c645c66-7c3d-11eb-1955-f7e59714787b
# ╠═d95d4c22-801f-11eb-14a2-833e11a65f3d
# ╠═a1f706d2-7c3d-11eb-3fbe-e3872d192c05
# ╠═0f1f308c-8020-11eb-3168-01946ee4b4d1
# ╠═6bcda4c6-7c3d-11eb-29b3-bbdb4044d641
# ╠═a4d50c28-7c3d-11eb-0cfb-27d78f900844
# ╠═0545ce32-7c46-11eb-3ff4-bf024fb73811
# ╠═7b570fb2-7c3e-11eb-09f8-5df49e7f5a9b
# ╠═5ea7d0a8-7c3f-11eb-0d4a-35ba92e372ad
# ╠═5a2eeb1a-7c3f-11eb-1e9e-8db2708de320
# ╠═ac2c05c0-8021-11eb-054a-d5a01e53b98e
# ╠═7d8b9fe6-7c40-11eb-1603-d7f60457bd55
# ╠═81664378-8021-11eb-3f3b-195bd4e83371
# ╠═8343b30c-7c42-11eb-3ae3-1ffe6a12f9d9
# ╠═8ce16c56-7c42-11eb-253d-9d1d6809379d
# ╠═adb3ebfc-7c42-11eb-10b6-898daac90b45
# ╠═b896fb5e-7c42-11eb-214b-9d348d1c2dba
# ╠═730b3aa2-7c40-11eb-249b-655ca7caa1b9
# ╠═3d9d470e-7c43-11eb-39e8-dd88a86846fd
# ╠═315336f4-7c5a-11eb-0a00-83e875d3c31c
# ╠═345ea3d0-7c44-11eb-0068-f7f694943b54
# ╠═c220d69e-7c43-11eb-05dd-21b150b8d436
# ╠═230d0878-7c45-11eb-0f2c-976edf65ab04
# ╠═2cb6793e-7c45-11eb-1b36-1b179a6b1168
# ╠═5d355172-7c40-11eb-1334-9fb53014f9a7
# ╠═ba6f73fe-7c4d-11eb-3a8a-912e71e751ce
# ╠═f8eea1ce-7c4d-11eb-2a18-21ae5d63457e
# ╠═83aaf5a6-7c5a-11eb-24ed-379e47866ce7
# ╠═9b21b178-7c5a-11eb-1ac7-ffd7ea79a5a1
# ╠═b4cd33e0-7c5a-11eb-1e70-ef71af2e3e2f
# ╠═f087d26e-7c4d-11eb-20e0-ff49a8cb1514
# ╠═df8801be-7c4d-11eb-3c88-4705e4dff1c2
# ╠═d4fd6f8e-7c4d-11eb-0512-49222239c4ce
# ╠═c28accfc-7c4d-11eb-0eb1-bffa0f2b0887
# ╠═81599eca-7c4d-11eb-0af0-b3b02e43f112
# ╠═76d5d6bc-7c4d-11eb-16ea-cf98397e0f90
# ╠═5aae7ad4-7c4d-11eb-1d94-a5a3ada7ecd6
# ╠═38700098-7c46-11eb-034c-d787179ef393
# ╠═fc8fff88-7c45-11eb-11c4-c95c7d35d388
# ╠═57101fac-7c40-11eb-3af8-557ad12ca677
# ╠═51d9b4bc-7c40-11eb-1412-5fae0ff8268d
# ╠═e156ecdc-7c3f-11eb-00c6-c35ad2ca2468
# ╠═ea4e82a4-64ad-11eb-05f9-0dcb561d1268
# ╟─4d2e00b4-6654-11eb-1885-b17c28c298e9
# ╠═c7b26d94-6653-11eb-3bc4-0710205a79f0
# ╟─2f049adc-6655-11eb-2f70-0fdab1a9b4ec
# ╠═a2bd25dc-6654-11eb-2628-b112c1bfa78d
