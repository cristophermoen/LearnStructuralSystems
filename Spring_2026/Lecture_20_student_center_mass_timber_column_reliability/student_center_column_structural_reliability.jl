using Distributions, CairoMakie

#Bloomberg Student Center 
#Mass timber column probability of failure 

num_of_samples = 100000000

crack_factor = 0.50
A_column = 16 * 16 * crack_factor #in^2


#column demand

#live load 
μ = 40.0 #psf #1 floor 
σ = μ * 0.15
dist = Normal(μ, σ)
p_LL = rand(dist, num_of_samples)

#dead load 
μ = 20.0 #psf #1 floor  
σ = μ * 0.05
dist = Normal(μ, σ)
p_DL = rand(dist, num_of_samples)

#total load 
p = p_LL .+ p_DL

#column tributary area 
trib_area = 40.0 * 40.0 #square ft 

#column demand 
Pu = p .* trib_area 


#Pine mass timber strength parallel to grain 
μ = 5000.0 #psi 
σ = μ * 0.18 
dist = Normal(μ, σ)
Fc_nominal = rand(dist, num_of_samples) * 0.5

#column capacity 
Pn = A_column * Fc_nominal

#calculate column probability of failure 
pf = length(findall(DC -> DC > 1.0, Pu ./ Pn)) / num_of_samples


#show demand and capacity histograms 
f = Figure()
ax = Axis(f[1, 1])
hist!(ax, Pu, label="column demand Pu")
hist!(ax, Pn, label="column capacity Pn")
axislegend(ax)
f


#ASCE 7-22
pf_goal_per_year = 2.5E-7 
pf_goal_100_years = 1 - exp(-pf_goal_per_year * 100)


