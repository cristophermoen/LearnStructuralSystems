#EN.560.301 Structural Systems I 
#HW6 Probability of Failure of Agora
#April 2023

#Calculate the probability of failure of a column in the Agora.

using LinearAlgebra, Distributions

num_samples = 100000000


#define random variables

#yield stress
μ_f_y = 50  #ksi
σ_f_y = 10 #ksi
f_y_dist = Normal(μ_f_y, σ_f_y)
f_y = rand(f_y_dist, num_samples)

#dead load
μ_DL = 20  #psf
σ_DL = 2 #psf
f_DL_dist = Normal(μ_DL, σ_DL)
f_DL = rand(f_DL_dist, num_samples)

#live load
μ_LL = 50  #psf
σ_LL = 12.5 #psf
f_LL_dist = Normal(μ_LL, σ_LL)
f_LL = rand(f_LL_dist, num_samples)

function solve_Pu(f_DL, f_LL, num_samples)
    Pu = zeros(Float64, num_samples)

    for i in 1:num_samples
        Pu[i] = ( f_DL[i] + f_LL[i] ) * 154.167 * 2.5 #consider two levels for building + a mezzannine
    end

    return Pu

end

Pu = solve_Pu(f_DL, f_LL, num_samples)

function solve_Pn(f_y, num_samples)
    num = zeros(Float64, num_samples)
    Fn = zeros(Float64, num_samples)
    Pn = zeros(Float64, num_samples)

    for i in 1:num_samples
        num[i] = f_y[i] / 186.9254335

        if num[i] <= 2.25
            Fn[i] = (0.658^(num[i]))*f_y[i]
            Pn[i] = Fn[i] * 11.1 * 1000
        else
            Fn[i] = 0.877 * 186.9254335
            Pn[i] = Fn[i] * 11.1 * 1000
        end
    end
    return Pn
end

Pn = solve_Pn(f_y, num_samples)


# count = []

function get_failure_count(Pn, Pu, num_samples)

    count = 0

    for i=1:num_samples
        if Pn[i] < Pu[i]
            count +=1 
        end
    end

    return count

end

w = get_failure_count(Pn, Pu, num_samples)

probability_of_failure = (w / num_samples) * 100