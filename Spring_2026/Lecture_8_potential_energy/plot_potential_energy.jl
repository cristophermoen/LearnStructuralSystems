
using CairoMakie 

E = 29000000.0
A = 1.0 
L = 100.0 
P = 10000.0 

Δ_solution = P * L / (A * E)


Δ = range(-0.5, 0.5, 20)

k = E * A / L

U = 1/2 * k .* Δ.^2
V = Δ .* P 

Π = U .- V 

scatterlines(Δ, Π)


