
using LinearAlgebra

σx = 0.0
σy = 0.0
τxy = 10.0

σ = [σx  τxy
     τxy σy]

principal_stresses = eigvals(σ)


principal_stress_directions = eigvecs(σ)

