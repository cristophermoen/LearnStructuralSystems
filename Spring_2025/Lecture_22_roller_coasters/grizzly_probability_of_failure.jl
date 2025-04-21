using Distributions, CairoMakie

num_of_samples = 100000000

#strut demand 
μ = 3400.0 #lbs 
σ= 200.0 #lbs 
dist = Normal(μ, σ)
strut_demand = rand(dist, num_of_samples)

#strut capacity 
strut_area = 1.5 * 7.5 * 6 #in^2 
Fc = 1000.0 #psi, parallel to grain, eastern pine No. 1 

#wet service factor Cm 
μ = 0.75  
σ= 0.1 
dist = Normal(μ, σ)
CM = rand(dist, num_of_samples)
strut_capacity = strut_area * Fc * CM 


#calculate strut probability of failure 
pf = length(findall(DC -> DC > 1.0, strut_demand ./ strut_capacity)) / num_of_samples


#show demand and capacity histograms 
f = Figure()
ax = Axis(f[1, 1])
hist!(ax, strut_capacity, color=:red)
hist!(ax, strut_demand, color=:blue)
f


#ASCE 7-22
pf_goal_per_year = 1E-7 
pf_goal_100_years = 1 - exp(-pf_goal_per_year * 100)


