#EN.560.301 Structural Systems I 
#Lecture 15
#March 2023

#Calculate the probability of failure of the Salginatoble bridge.

using LinearAlgebra, Distributions

include("helper_functions.jl")

num_samples = 100000

μ_f_c = 15.0 * 6.9 * 1000  #kPa 
σ_f_c = 0.2 * μ_f_c
f_c_dist = Normal(μ_f_c, σ_f_c)
f_c = rand(f_c_dist, num_samples)

μ_DL_LL = 1.0
σ_DL_LL = 0.10
DL_LL_dist = Normal(μ_DL_LL, σ_DL_LL)
load_scale = rand(DL_LL_dist, num_samples)

results = calculate_probability_of_failure(f_c, load_scale)



