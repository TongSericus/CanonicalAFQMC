function time_slice_spectrum()
    if flag_background_subtraction == 1
        residue_term = exp(-Δτ * sum(background_dist .^2) * U)
        H_matrix = kinetic_matrix_imgtime * auxiliary_field_matrix() * kinetic_matrix_imgtime * residue_term
    else
        H_matrix = kinetic_matrix_imgtime * auxiliary_field_matrix() * kinetic_matrix_imgtime
    end
    return eigvals(H_matrix)
end

function sub_partition_weight(L)
    zₖ = complex(zeros(num_bosons))
    element_product = complex(ones(num_sites))
    for i = 1 : L
        element_array = time_slice_spectrum()
        for j = 1 : num_sites
            element_product[j] *= element_array[j]
        end
    end
    for k = 1 : num_bosons
        zₖ[k] = sum(element_product .^ k)
    end
    return zₖ
end

function partition_recursive_calculation(L)
    zₖ = sub_partition_weight(L)
    Z = complex(zeros(num_bosons))
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
    return Z[num_bosons]
end
