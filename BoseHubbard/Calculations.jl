using LinearAlgebra

function time_slice_spectrum()
    H_matrix = kinetic_matrix_imgtime * auxiliary_field_matrix() * kinetic_matrix_imgtime
end

function sub_partition_weight(L)
    matrix_product = Matrix{Complex{Float64}}(I, num_sites, num_sites)
    for i = 1 : L
        matrix_product *= time_slice_spectrum()
    end
    #return matrix_product
    #if flag_background_subtraction == 1
        #matrix_product *= exp(sum(β * (U/2) * background_dist.^2))
        #return matrix_product
    #end
    effective_spectrum, u = eigen(matrix_product) # spectrum element: e^(-βϵᵢ)
    u_inv = inv(u)
    U_array = complex(zeros(num_sites, num_sites, num_sites))
    for i = 1 : num_sites
        for j = 1 : num_sites
            for k = 1 : num_sites
                U_array[i, j, k] = u_inv[k, i] * u[j, k]
            end
        end
    end
    return effective_spectrum, U_array
end

function recursion_partition_function(effective_spectrum)
    Z = complex(zeros(num_bosons))
    zₖ = complex(zeros(num_bosons))
    for k = 1 : num_bosons
        zₖ[k] = sum(effective_spectrum .^ k)
    end
    for N = 1 : num_bosons
        for M = 1 : N
            if M == N
                Z[N] += 1 * zₖ[M]
            else
                Z[N] += Z[N-M] * zₖ[M]
            end
        end
        Z[N] = Z[N] / N
    end
    return Z
end

function recursion_OccNum(Z, effective_spectrum) # Calculating <nᵢ>
    n_i = complex(zeros(num_bosons, num_sites))
    n_i[1, :] = effective_spectrum
    for N = 2 : num_bosons
        n_i[N, :] = effective_spectrum .* (Z[N-1] .+ n_i[N-1, :])
    end
    return n_i
end

function recursion_double_OccNum_eq(Z, effective_spectrum, n_i) # Calculating <nᵢ^2>
    n²_i = complex(zeros(num_bosons, num_sites)) # <nᵢ^2> array
    n²_i[1, :] = effective_spectrum
    for N = 2 : num_bosons
        n²_i[N, :] = effective_spectrum .* (Z[N-1] .+ 2 * n_i[N-1, :] .+ n²_i[N-1, :])
    end
    return n²_i[num_bosons, :]
end

function recursion_double_OccNum_neq(Z, effective_spectrum, n_i) # Calculating <nᵢnⱼ>
    n²_i = complex(zeros(num_bosons, num_sites, num_sites)) # <nᵢnⱼ> array
    for i = 1 : num_sites
        for j = i + 1 : num_sites
            for N = 2 : num_bosons
                n²_i[N, i, j] = effective_spectrum[j] / 2 * (n_i[N-1, i] + n²_i[N-1, i, j]) +
                            + effective_spectrum[i] / 2 * (n_i[N-1, j] + n²_i[N-1, i, j])
            end
            n²_i[num_bosons, j, i] = n²_i[num_bosons, i, j]
        end
    end
    return n²_i[num_bosons, :, :]
end

#function measure_condensate_fraciton()
#end

function measure_TotEnergy(num_iterations)
    partition_func_array = Vector{Complex{Float64}}()
    Ek_array = Vector{Complex{Float64}}()
    Ep_array = Vector{Complex{Float64}}()
    E_tot_array = Vector{Complex{Float64}}()
    for N = 1 : num_iterations
        effective_spectrum, U_array = sub_partition_weight(L)
        Z = recursion_partition_function(effective_spectrum)
        n_i = recursion_OccNum(Z, effective_spectrum)
        n²_i_diag = recursion_double_OccNum_eq(Z, effective_spectrum, n_i)
        n²_i_offdiag = recursion_double_OccNum_neq(Z, effective_spectrum, n_i)
        E_onebody, E_twobody = 0, 0
        # one-body energy calculation
        for i = 1 : num_sites
            for j = 1 : num_sites
                if abs(i - j) == 1 || abs(i - j) == num_sites - 1
                    D_ij = sum(U_array[i, j, :] .* n_i[num_bosons, :])
                    E_onebody += -t * D_ij
                end
            end
        end
        # two-body energy calculation
        for i = 1 : num_sites
            U_matrix = complex(zeros(num_sites, num_sites))
            D_iiii = 0
            for p = 1 :num_sites
                for r = 1 : num_sites
                    if p != r
                        U_matrix[p, r] = U_array[i, i, p] * U_array[i, i, r]
                        D_iiii += U_matrix[p, r] * n_i[num_bosons, p]
                    end
                end
                D_iiii += U_array[i, i, p]^2 * n²_i_diag[p]
            end
            D_iiii += sum( 2 * U_matrix .* n²_i_offdiag)
            E_twobody += U/2 * D_iiii
        end
        push!(partition_func_array, Z[num_bosons])
        push!(Ek_array, E_onebody)
        push!(Ep_array, E_twobody)
        push!(E_tot_array, E_onebody + E_twobody)
    end
    Z_mean = mean(partition_func_array)
    Ek_mean = mean(Ek_array) / Z_mean
    Ek_error = std(Ek_array) / (sqrt(num_iterations) * real(Z_mean))
    Ep_mean = mean(Ep_array) / Z_mean
    Ep_error = std(Ep_array) / (sqrt(num_iterations) * real(Z_mean))
    E_tot_mean = mean(E_tot_array) / Z_mean
    E_tot_error_bar = std(E_tot_array) / (sqrt(num_iterations) * real(Z_mean))
    return Ek_mean, Ek_error, Ep_mean, Ep_error, E_tot_mean, E_tot_error_bar
end
