
using LinearAlgebra

σx = 0.0
σy = 10.0
τxy = 0.0

σ = [σx  τxy
     τxy σy]

principal_stresses = eigvals(σ)


principal_stresses = eigvecs(σ)

