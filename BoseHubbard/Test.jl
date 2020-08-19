function test_density_matrix()
    effective_spectrum, U_array = sub_partition_weight(L)
    Z = recursion_partition_function(effective_spectrum)
    D = complex(zeros(num_sites, num_sites))
    for i = 1 : num_sites
        for j = 1 : num_sites
            n_i = recursion_OccNum(Z, effective_spectrum)
            D[i, j] = sum(U_array[i, j, :] .* n_i[num_bosons, :])
        end
    end
    return Z[num_bosons], D
end

function test_density_matrix_sampling(num_iterations) # tests for one-body density matrices
    partition_func_array = Vector{Complex{Float64}}()
    density_mat_array = Matrix{Complex{Float64}}[]
    for N = 1 : num_iterations
        partition_func_sample, density_mat_sample = test_density_matrix()
        push!(partition_func_array, partition_func_sample)
        push!(density_mat_array, density_mat_sample)
    end
    Z_mean = mean(partition_func_array)
    Z_error_bar = stdm(partition_func_array, Z_mean) / sqrt(num_iterations)
    D_mean = mean(density_mat_array) / Z_mean
    return Z_mean, Z_error_bar, D_mean
end

function test_matrix_product(num_iterations)
    mat_product_array = Vector{Array{Complex{Float64},1}}()
    for N = 1 : num_iterations
        mat_product_sample = sub_partition_weight(L)
        push!(mat_product_array, mat_product_sample)
    end
    mat_product_mean = mean(mat_product_array)
end

function test_partition_function(num_iterations)
    partition_func_array = Vector{Complex{Float64}}()
    for i = 1 : num_iterations
        effective_spectrum, U_array = sub_partition_weight(L)
        Z = recursion_partition_function(effective_spectrum)
        push!(partition_func_array, Z[num_bosons])
    end
    Z_mean = mean(partition_func_array)
    Z_error_bar = stdm(partition_func_array, Z_mean) / num_iterations
    return Z_mean, Z_error_bar
end