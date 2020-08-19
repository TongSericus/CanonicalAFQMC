function pilot_sampling(num_PilotSamples)
    shift_array_SpinUp = Vector{Float64}()
    shift_array_SpinDn = Vector{Float64}()
    output_file = open("result","a+")
    for N = 1 : num_PilotSamples
        eff_spec_SpinUp, eff_spec_SpinDn = matrix_product()
        Z_up, Z_dn = recursion_partition_function(eff_spec_SpinUp, eff_spec_SpinDn)
        write(output_file, string(Z_up[num_SpinUp] * Z_dn[num_SpinDn]), "\n")
    end
    close(output_file)
end

"Functions for Importance Sampling"

"function roulette_wheel_generator(γ, num_sites)
    roulette_wheel = Vector{Float64}()
    cumulative_array = Vector{Float64}()
    W = 0 # stands for total weight
    for i = 0 : num_sites
        selection_length = exp( -2 * i * γ + 2 * (num_sites - i) * γ)
        W += selection_length
        push!(roulette_wheel, selection_length)
        push!(cumulative_array, W)
    end
    roulette_wheel /= W
    cumulative_array /= W
    return W, roulette_wheel, cumulative_array
end

function roulette_selection(r, W, cumulative_array) # Binary Search
    L = length(cumulative_array)
    i_index, f_index = 1, L
    p_index = floor(Int64, (i_index + f_index) / 2) # pointer index
    while true
        if cumulative_array[p_index] <= r <= cumulative_array[p_index + 1]
            break
        elseif p_index == 1 && r <= cumulative_array[p_index]
            return p_index
        end
        if r < cumulative_array[p_index]
            f_index = p_index
            p_index = floor(Int64, (i_index + f_index) / 2)
        else
            i_index = p_index
            p_index = floor(Int64, (i_index + f_index) / 2)
        end
    end
    p_index + 1
end"

function σlistGenerator(num_sites)
    σlist = Array{Int64,1}[]
    for i = 0 : num_sites
        σarray = 2 * vcat(zeros(Int16, i), ones(Int16, num_sites - i)) .- 1
        push!(σlist, σarray)
    end
    σlist
end

function σfieldGenerator(k, σlist, num_sites) # i takes value from 0 to num_sites
    σfield = Random.shuffle(σlist[k + 1])
    WeightPerStep = (num_sites + 1) * binomial(num_sites, k) / 2 ^ num_sites
    σfield, WeightPerStep
end
