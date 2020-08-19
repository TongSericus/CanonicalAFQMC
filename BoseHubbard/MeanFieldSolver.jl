using LinearAlgebra

include("MatrixGenerator.jl")

function mean_field_solver(β, t, U, num_sites)
    convergence_tolerance = 1e-6
    kinetic_val = kinetic_matrix(num_sites, t)
    n_set = ones(num_sites) / num_sites
    # each element in the array represents the density of particles on each site
    trial_matrix = kinetic_val + Diagonal(n_set * U)
    n_set_old = zeros(num_sites)
    while sum(abs.(n_set_old - n_set)) > convergence_tolerance
        energy_spectrum, basis_vector = eigen(trial_matrix)
        density_matrix = zeros(num_sites, num_sites)
        n_set_old = n_set
        Z = sum(exp.(- β * energy_spectrum))
        for i = 1:num_sites
            density_matrix += (exp.(- β * energy_spectrum[i])  / Z) * (basis_vector[:, i] * basis_vector[:, i]')
        end
        n_set = diag(density_matrix)
        trial_matrix = kinetic_val + Diagonal(n_set * U)
    end
    return n_set
end
