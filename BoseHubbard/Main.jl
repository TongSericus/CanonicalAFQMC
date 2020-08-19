using Statistics
using DelimitedFiles

include("MatrixGenerator.jl")
include("Calculations.jl")
include("MeanFieldSolver.jl")
include("Test.jl")

# get inputs, initialization
lines = readlines("InputBoseHubbard")
# inverse temperature
temp_name, β_str = split(lines[1])
const β = parse(Float64, β_str)
# number of sites
temp_name, num_sites_str = split(lines[2])
const num_sites = parse(Int64, num_sites_str)
# number of bosons
temp_name, num_bosons_str = split(lines[3])
const num_bosons = parse(Int64, num_bosons_str)
# length of the imaginary time interval
temp_name, Δτ_str = split(lines[4])
const Δτ = parse(Float64, Δτ_str)
# hopping constant
temp_name, t_str = split(lines[5])
const t = parse(Float64, t_str)
# on-site repulsion constant
temp_name, U_str = split(lines[6])
const U = parse(Float64, U_str)
# number of iterations
temp_name, num_iterations_str = split(lines[7])
const num_iterations = convert(Int64, parse(Float64, num_iterations_str))
# flag for the background subtraction
temp_name, flag_background_subtraction_str = split(lines[8])
const flag_background_subtraction = parse(Int16, flag_background_subtraction_str)
# flag for calculating the density matrix
temp_name, flag_density_matrix_str = split(lines[9])
const flag_density_matrix = convert(Int64, parse(Int16, flag_density_matrix_str))

const background_dist = ones(num_sites) / num_sites #  Equally distributed among sites
const residue_term = sum(β * (U/2) * background_dist.^2)

# Constants for calculations
# Kinetic matrix in each imaginary time slice is fixed throughout the calculation.
const kinetic_matrix_imgtime = complex(exp(-kinetic_matrix(num_sites, t) * Δτ / 2))
# number of imaginary time steps
const L = convert(Int64, round(β / Δτ))
# coefficient of the H-S transformed potential term
const AF_coefficient = sqrt(complex(-Δτ * U))

#mean_value, error_bar = Z_iteration(num_iteration, L)

#write_to_file(mean_value, error_bar, Δτ)

#function write_to_file(mean_value, error_bar, Δτ)
    #output_file = open("result","a+")
    #write(output_file, string(β), " ", string(Δτ), " ", string(mean_value), " ", string(error_bar), "\n")
    #close(output_file)
#end

Ek_mean, Ek_error, Ep_mean, Ep_error, E_tot_mean, E_tot_error_bar = measure_TotEnergy(num_iterations)
println("U = ", U)
println("Beta = ", β)
println("Ns = ", num_sites)
println("NBosons = ", num_bosons)
println("Ek = ", Ek_mean, ", Error = ", Ek_error)
println("Ep = ", Ep_mean, ", Error = ", Ep_error)
println("Etot = ", E_tot_mean, ", Error = ", E_tot_error_bar)