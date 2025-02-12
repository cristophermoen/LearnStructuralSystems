using ForwardDiff

E = 29000000.0
I = 1.0
L = 60.0
a = L / 1000
P = 10000.0

# x_range = range(0.0, L, 11)

# yo = a * sin.(π .* x / L)

y_basis = sin.(π * x_range ./ 60.0)

f(x::Real) = sin(π * x / 60.0);  # returns a vector
ddy_basis = [ForwardDiff.derivative(x -> ForwardDiff.derivative(f, x), x_range[i]) for i in eachindex(x_range)]

Dxx = ddy_basis[2] / y_basis[2]

Δ = (-P * a) ./ (E * I * Dxx .+ P)


Pe = π^2 * E * I / L^2
Δ_theory = a - a / (1 - P/Pe) 




