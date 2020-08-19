using LinearAlgebra
using Distributions

function kinetic_matrix(num_sites, t)   # one-body matrix, periodic boundary condition
    K_matrix = zeros(num_sites, num_sites)
    for i =  1 : num_sites
        for j = 1 : num_sites
            if abs(i - j) == 1 || abs(i - j) == num_sites - 1
                K_matrix[i, j] = -t
            end
        end
    end
    if flag_background_subtraction == 1
        K_matrix += Diagonal(background_dist) * U
    end
    return K_matrix
end

function auxiliary_field_matrix()
    ϕ = rand(Normal(), num_sites)
    coefficent_vector = ϕ * AF_coefficient
    AF_matrix = complex(zeros(num_sites, num_sites))
    for i = 1 : num_sites
        AF_matrix[i, i] = coefficent_vector[i]
    end
    if flag_background_subtraction == 1
        background_exponent = sum(-coefficent_vector .* background_dist)
        background_weight = exp(background_exponent)
        return background_weight * exp(AF_matrix)
    elseif flag_background_subtraction == 0
        return exp(AF_matrix)
    end
end
