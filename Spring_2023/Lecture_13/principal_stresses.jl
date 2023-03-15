using LinearAlgebra

σ = [0.0 20.0
     20.0 0.0]

principal_stresses = eigvals(σ) 

principal_stress_directions = eigvecs(σ) 