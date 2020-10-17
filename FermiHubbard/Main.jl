using Statistics
using DelimitedFiles
using Random
using Printf

struct Lattice
    nsites::Int64
    nup::Int64
    ndn::Int64
    ntot::Int64
    t::Float64
    U::Float64
end

struct QMC
    nblocks::Int64
    samples::Int64
    Δτ::Float64
    L::Int64
    γ::Float64
    Kmatrix::Array{Float64,2}
end

struct Opt  # Optimization
    stab_interval::Int64
    nburn_in::Int64
end

const FloatType = Union{Float64, Complex{Float64}}

include("./MatrixGenerator.jl")
include("./Calculations.jl")

function getinputs()
    lines = readlines("./InputFermiHubbard")
    
    # basic parameters
    β = parse(Float64, split(lines[1])[2])
    num_sites = parse(Int64, split(lines[2])[2])
    num_spinup = parse(Int64, split(lines[3])[2])
    num_spindn = parse(Int64, split(lines[4])[2])
    num_elec = num_spinup + num_spinup

    Δτ = parse(Float64, split(lines[5])[2])
    t = parse(Float64, split(lines[6])[2])
    U = parse(Float64, split(lines[7])[2])
    num_blocks = convert(Int64, parse(Float64, split(lines[8])[2]))
    num_samples = convert(Int64, parse(Float64, split(lines[9])[2]))
    
    # optimization parameters
    stab_interval = parse(Int64, split(lines[12])[2])
    num_burnins = convert(Int64, parse(Float64, split(lines[13])[2]))

    # constants for calculations
    L = convert(Int64, round(β / Δτ))
    γ = atanh(sqrt(tanh(Δτ * U / 4)))

    lattice = Lattice(num_sites, num_spinup, num_spindn, num_elec, t, U)

    K_matrix = kinetic_matrix(lattice)
    qmc = QMC(num_blocks, num_samples, Δτ, L, γ, exp(-K_matrix * Δτ / 2))

    opt = Opt(stab_interval, num_burnins)
    return lattice, qmc, opt
end

function mainbody()
    lattice, qmc, opt = getinputs()
    SubEk_array = Vector{Float64}()
    SubEp_array = Vector{Float64}()
    SubE_array = Vector{Float64}()
    SubEk_error_array = Vector{Float64}()
    SubEp_error_array = Vector{Float64}()
    SubE_error_array = Vector{Float64}()
    filename = string("output_Ns", lattice.nsites, "_U", lattice.U, "_beta", qmc.Δτ * qmc.L)
    outputfile = open(filename, "a+")
    @printf(outputfile, "       Canonical Ensemble AFQMC for Hubbard Model     \n\n")
    @printf(outputfile, "Number of Electrons: %d\n", lattice.ntot)
    @printf(outputfile, "Number of Blocks: %d\n", qmc.nblocks)
    @printf(outputfile, "Length of Img Timestep: %.4f\n", qmc.Δτ)
    @printf(outputfile, "Stablization Interval: %d\n\n", opt.stab_interval)
    for i = 1 : qmc.nblocks
        SubE, SubE_error, SubEk, SubEk_error, SubEp, SubEp_error, sgn_mean = MCIntegration(lattice, qmc, opt)
        SubE = SubE / sgn_mean
        SubEk = SubEk / sgn_mean
        SubEp = SubEp / sgn_mean

        @printf(outputfile, "Block # %d out of %d\n", i, qmc.nblocks)
        @printf(outputfile, "--------------------\n")
        @printf(outputfile, "Average Sign: %.2f\n", sgn_mean)
        @printf(outputfile, "Ek            Ep            Etot\n")
        @printf(outputfile, "%.5f     %.5f       %.5f\n", SubEk, SubEp, SubE)
        @printf(outputfile, "Error:\n")
        @printf(outputfile, "%.5f     %.5f       %.5f\n\n", SubEk_error, SubEp_error, SubE_error)

        push!(SubEk_array, SubEk)
        push!(SubEk_error_array, SubEk_error)
        push!(SubEp_array, SubEp)
        push!(SubEp_error_array, SubEp_error)
        push!(SubE_array, SubE)
        push!(SubE_error_array, SubE_error)
    end
    Ek_mean, Ek_error = mean(SubEk_array), mean(SubEk_error_array) / sqrt(qmc.nblocks)
    Ep_mean, Ep_error = mean(SubEp_array), mean(SubEp_error_array) / sqrt(qmc.nblocks)
    E_mean, E_error = mean(SubE_array), mean(SubE_error_array) / sqrt(qmc.nblocks)

    @printf(outputfile, "--------------------\n")
    @printf(outputfile, "Results for this run:\n")
    @printf(outputfile, "Ek            Ep            Etot\n")
    @printf(outputfile, "%.5f     %.5f       %.5f\n", Ek_mean, Ep_mean, E_mean)
    @printf(outputfile, "Error:\n")
    @printf(outputfile, "%.5f     %.5f       %.5f\n\n", Ek_error, Ep_error, E_error)
    close(outputfile)
end
