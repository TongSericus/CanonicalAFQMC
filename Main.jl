using Distributions
using Statistics
using LinearAlgebra
using DelimitedFiles

include("MatrixGenerator.jl")
include("Calculation.jl")
include("SaplingAndOutput.jl")
include("MeanFieldSolver.jl")

# get inputs, initialization
lines = readlines("input_hubbard_boson")
temp_name, β_str = split(lines[1])
const β = parse(Float64, β_str)
temp_name, num_sites_str = split(lines[2])
const num_sites = parse(Int16, num_sites_str)
temp_name, num_bosons_str = split(lines[3])
const num_bosons = parse(Int16, num_bosons_str)
temp_name, Δτ_str = split(lines[4])
const Δτ = parse(Float64, Δτ_str)
temp_name, t_str = split(lines[5])
const t = parse(Float64, t_str)
temp_name, U_str = split(lines[6])
const U = parse(Float64, U_str)
temp_name, flag_background_subtraction_str = split(lines[7])
const flag_background_subtraction = parse(Int16, flag_background_subtraction_str)

const background_dist = mean_field_solver()

const kinetic_matrix_imgtime = complex(exp(k_background_subtraction(kinetic_matrix()) * Δτ / 2))
# Kinetic matrix in each imaginary time slice is fixed through the calculation.
const L = convert(Int64, round(β / Δτ))
# number of imaginary time steps
const AF_coefficient = sqrt(complex(-Δτ * U))
# coefficient of the mean-field-like term after H-S transformation

mean_value, error_bar = Z_iteration(L)

write_to_file(mean_value, error_bar, Δτ)
