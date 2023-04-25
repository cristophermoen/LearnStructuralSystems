#EN.560.301 Structural Systems I 
#World Trade Center study
#April 2023



using LinearAlgebra, Distributions, AISIS100



floors_above = 22
surface_area = 43680.0 - 11000.0 #sq. ft 
num_perimeter_columns = 236

A_column = 20.0*20.0 - 19*19  #in^2
I_column = 2473.0 * 2 #in^4
L_column = 14.0 *12 #in
E = 29000.0 #ksi 


#yield_stress 
μ_f_y = 36  #ksi
σ_f_y = 5 #ksi
f_y_dist = Normal(μ_f_y, σ_f_y)


function calculate_column_demand(DL, LL, surface_area, floors_above, num_perimeter_columns)
  
    Pu = (DL + LL) * surface_area * floors_above / 2 / num_perimeter_columns / 1000  #divide by 2, half goes to core, half goes to perimeter 

    return Pu

end


function calculate_column_strengths_at_floor(Fcre, f_y, A_column)

    Pn = zeros(Float64, length(f_y))

    for i in eachindex(f_y)
        Pn[i], not_used = AISIS100.v16.e2(Fcre=Fcre, Fy=f_y[i], Ag=A_column, design_code="nominal")
    end

    return Pn 

end


function calculate_gravity_system_failure_probability(f_DL, f_LL, f_y, num_samples, surface_area, floors_above, num_perimeter_columns, A_column, E, I_column, L_column, num_failed_columns_limit)

    count = 0

    for i=1:num_samples

        Pu = calculate_column_demand(f_DL[i], f_LL[i], surface_area, floors_above, num_perimeter_columns)

        Fcre = π^2 * E * I_column/(L_column^2 * A_column)
        Pn = calculate_column_strengths_at_floor(Fcre, f_y, A_column)

        Pn = A_column .* f_y

        column_failure_count = length(findall(R->R < Pu, Pn))

        if column_failure_count >= num_failed_columns_limit

            count += 1

        end

    end

    probability_of_failure = count / num_samples

    return probability_of_failure

end

num_samples = 100000

#dead load
μ_DL = 20  #psf
σ_DL = μ_DL * 0.05 #psf
f_DL_dist = Normal(μ_DL, σ_DL)
f_DL = rand(f_DL_dist, num_samples)

#live load
μ_LL = 40  #psf
σ_LL = 0.15 * μ_LL #psf
f_LL_dist = Normal(μ_LL, σ_LL)
f_LL = rand(f_LL_dist, num_samples)


f_y = rand(f_y_dist, num_perimeter_columns)

#probability of gravity system failure at 80 stories up 
#no damage, no fire 

num_failed_columns_limit = 1
pf_1 = calculate_gravity_system_failure_probability(f_DL, f_LL, f_y, num_samples, surface_area, floors_above, num_perimeter_columns, A_column, E, I_column, L_column, num_failed_columns_limit)


#probability of gravity system failure at 80 stories up 
#29 columns removed 
num_perimeter_columns = 236 - 29

num_failed_columns_limit = 1
pf_2 = calculate_gravity_system_failure_probability(f_DL, f_LL, f_y, num_samples, surface_area, floors_above, num_perimeter_columns, A_column, E, I_column, L_column, num_failed_columns_limit)


#probability of gravity system failure at 80 stories up 
#29 columns removed + fire  
num_perimeter_columns = 236 - 29
E = 29000.0 * 0.05
f_y = f_y * 0.05

num_failed_columns_limit = 1
pf_3 = calculate_gravity_system_failure_probability(f_DL, f_LL, f_y, num_samples, surface_area, floors_above, num_perimeter_columns, A_column, E, I_column, L_column, num_failed_columns_limit)


(pf_1, pf_2, pf_3)