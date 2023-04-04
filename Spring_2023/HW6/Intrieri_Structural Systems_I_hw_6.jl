#EN.560.301 Structural Systems I 
#Lecture 17
#April 2023

#Calculate the probability of failure of the Agora Institute.

using LinearAlgebra, Distributions


num_samples = 10000000

#Find the resistance of a column
A = 13.6 #in^2
E = 29000000 #psi
L = 12*12 #in
Iyy = 262 #in^4
r = sqrt(Iyy / A)
Fe = π^2 * E / (L/r)^2 #psi
Fy_nominal = 50000 #psi
fractional_power = Fy_nominal / Fe

Fn_nominal = (0.658 ^ fractional_power) * Fy_nominal

μ_f_c = Fn_nominal * A #psi
σ_f_c = 0.2 * μ_f_c
f_c_dist = Normal(μ_f_c, σ_f_c)
f_c = rand(f_c_dist, num_samples)

#Find the load placed on a column
μ_DL = 20 * 1000/12 #psf*sf = pounds
σ_DL = 1 * 1000/12

μ_LL = 50 * 1000/12
σ_LL = 2.5 * 1000/12


μ_DL_LL = μ_DL + μ_LL
σ_DL_LL = sqrt(σ_DL^2 + σ_LL^2)

DL_LL_dist = Normal(μ_DL_LL, σ_DL_LL)
load_scale = rand(DL_LL_dist, num_samples)

function calculate_probability_of_failure(f_c, load_scale)
    number_of_failures = 0
    for i in eachindex(f_c)
        if f_c[i] < load_scale[i]
            number_of_failures += 1
        end
    end
    Pf = number_of_failures / num_samples
    return Pf
end

results = calculate_probability_of_failure(f_c, load_scale)