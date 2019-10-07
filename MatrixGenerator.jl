function kinetic_matrix()   # one-body matrix, periodic boundary condition
    K_matrix = zeros(num_sites, num_sites)
    for i =  1 : num_sites
        for j = 1 : num_sites
            if abs(i - j) == 1 || abs(i - j) == num_sites - 1
                K_matrix[i, j] = -t
            end
        end
    end
    return K_matrix
end

function k_background_subtraction(kinetic_matrix)
    if flag_background_subtraction == 1
        kinetic_matrix += 2 * Diagonal(background_dist) * U
    else
    return kinetic_matrix
    end
end

function auxiliary_field_matrix()
    ϕ = rand(Normal(), num_sites)
    coefficent_vector = ϕ * AF_coefficient
    AF_matrix = complex(zeros(num_sites, num_sites))
    for i = 1 : num_sites
        AF_matrix[i, i] = coefficent_vector[i]
    end
    if flag_background_subtraction == 1
        background_exponent = exp.(-coefficent_vector .* background_dist)
        background_weight = prod(background_exponent)
        return background_weight * exp(AF_matrix)
    else
        return exp(AF_matrix)
    end
end
