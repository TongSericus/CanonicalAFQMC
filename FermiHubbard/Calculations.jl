using LinearAlgebra

function recursion_partition_function(num_elec, eff_spec_SpinUp, eff_spec_SpinDn)
    Z_up, Z_dn = zeros(num_elec), zeros(num_elec)
    zₖ_up, zₖ_dn = zeros(num_elec), zeros(num_elec)
    for k = 1 : num_elec
        zₖ_up[k] = sum(eff_spec_SpinUp .^ k)
        zₖ_dn[k] = sum(eff_spec_SpinDn .^ k)
    end
    for N = 1 : num_elec
        for k = 1 : N
            sgn_k = 2 * isodd(k) - 1
            if k == N
                Z_up[N] += sgn_k * zₖ_up[k]
                Z_dn[N] += sgn_k * zₖ_dn[k]
            else
                Z_up[N] += sgn_k * Z_up[N-k] * zₖ_up[k]
                Z_dn[N] += sgn_k * Z_dn[N-k] * zₖ_dn[k]
            end
        end
        Z_up[N] = Z_up[N] / N
        Z_dn[N] = Z_dn[N] / N
    end
    return Z_up, Z_dn
end

function recursion_OccNum(Z_up, eff_spec_SpinUp, Z_dn, eff_spec_SpinDn, num_elec, num_sites) # Calculating <nᵢ> * Z_N
    n_up_mean = zeros(num_elec, num_sites)
    n_dn_mean = zeros(num_elec, num_sites)
    n_up_mean[1, :], n_dn_mean[1, :] = eff_spec_SpinUp, eff_spec_SpinDn
    for N = 2 : num_elec
        n_up_mean[N, :] = eff_spec_SpinUp .* (Z_up[N-1] .- n_up_mean[N-1, :])
        n_dn_mean[N, :] = eff_spec_SpinDn .* (Z_dn[N-1] .- n_dn_mean[N-1, :])
    end
    return n_up_mean ./ Z_up, n_dn_mean ./ Z_dn
end

function measure_tot_energy(num_sites, t, U, U_array_SpinUp, U_array_SpinDn, n_up_mean, n_dn_mean, num_SpinUp, num_SpinDn, Z_up, Z_dn)
    E_onebody, E_twobody = 0, 0
    for i = 1 : num_sites
        for j = 1 : num_sites
            if abs(i - j) == 1 || abs(i - j) == num_sites - 1
                D_ij_up = sum(U_array_SpinUp[i, j, :] .* n_up_mean[num_SpinUp, :])
                D_ij_dn = sum(U_array_SpinDn[i, j, :] .* n_dn_mean[num_SpinDn, :])
                E_onebody += -t * ( D_ij_up + D_ij_dn )
            end
        end
    end
    for i = 1 : num_sites
        U_matrix = U_array_SpinUp[i, i, :] * U_array_SpinDn[i, i, :]'
        n_matrix = n_up_mean[num_SpinUp, :] * n_dn_mean[num_SpinDn, :]'
        D_iiii = sum(U_matrix .* n_matrix)
        E_twobody += U * D_iiii
    end
    return E_onebody, E_twobody
end

# Perfrom a truncated version of recursion; Added 10/06/2020
function SpectrumTruncation(nsites::Int64, nelec::Int64, expϵ::Array{Float64,1})
    expϵ_ratio = expϵ[2 : nsites] ./ expϵ[1 : nsites - 1]
    γᵤ, γₗ = nsites - nelec + 1, nsites - nelec + 1
    expgapᵤ, expgapₗ = 1, 1

    while true
        expgapᵤ *= expϵ_ratio[γᵤ - 1]
        if expgapᵤ > 1e4
            break
        elseif γᵤ == 2
            γᵤ -= 1
            break
        end
        γᵤ -= 1
    end
    
    while true
        expgapₗ *= expϵ_ratio[γₗ]
        if expgapₗ > 1e4
            break
        elseif γₗ + 1 == nsites
            γₗ += 1
            break
        end
        γₗ += 1
    end

    return γᵤ, γₗ
end

function RecursionPFSgn(nelecs::Int64, expϵ::Array{Float64,1})
    Z, zₖ = zeros(nelecs), zeros(nelecs)
    for k = 1 : nelecs
        zₖ[k] = sum(expϵ .^ k)
    end
    for N = 1 : nelecs
        for k = 1 : N
            sgnk = 2 * isodd(k) - 1
            if k == N
                Z[N] += sgnk * zₖ[k]
            else
                Z[N] += sgnk * Z[N-k] * zₖ[k]
            end
        end
        Z[N] = Z[N] / N
    end
    return Rsign(Z[nelecs])
end

# A new recursion function that deals with the truncated spectrum;
function TruncatedRecursion(nelecs::Int64, expϵ::Array{Float64,1})
    NSites = length(expϵ)
    n_mean = zeros(Float64, NSites, NSites)
    Zₙratio = zeros(Float64, nelecs)

    Zₙratio[1] = sum(expϵ)
    n_mean[:, 1] = expϵ / Zₙratio[1]

    if nelecs == 1
        return n_mean, Zₙratio
    else
        for N = 2 : nelecs
            Zₙratio[N] = sum(expϵ .* (1 .- n_mean[:, N-1])) / N
            n_mean[:, N] = expϵ .* (1 .- n_mean[:, N - 1]) / Zₙratio[N]
        end
    end
    return n_mean, Zₙratio
end

function SpectrumRemodification(nsites::Int64, nspins::Int64, expϵ::Array{Float64,1})
    expϵperm = sortperm(abs.(expϵ))
    absexpϵ = abs.(expϵ)[expϵperm]

    γᵤ, γₗ = SpectrumTruncation(nsites, nspins, absexpϵ)
    n_mean, Zₙratio = TruncatedRecursion(nspins - (nsites - γₗ), absexpϵ[γᵤ : γₗ])
    
    if nsites > γₗ
        trun_n_mean, trun_Zₙratio = TruncatedRecursion(nsites - γₗ, absexpϵ)
        Zₙratio = vcat(trun_Zₙratio, Zₙratio)
    end

    n_mean = vcat(zeros(γᵤ - 1, γₗ - γᵤ + 1), n_mean)
    n_mean = hcat(zeros(γₗ, nsites - γₗ), n_mean)
    n_mean = vcat(n_mean, ones(nsites - γₗ, nsites - γᵤ + 1))

    if nsites - γₗ == 0
        SgnZₙ = 1
    else
        SgnZₙ = RecursionPFSgn(nsites - γₗ, expϵ ./ maximum(expϵ))
    end
    SgnZₙRest = RecursionPFSgn(nspins - (nsites - γₗ), expϵ[expϵperm][γᵤ : γₗ] ./ maximum(expϵ[expϵperm][γᵤ : γₗ]))

    #return n_mean[:, nspins][sortperm(expϵperm)], Zₙratio, SgnZₙ * SgnZₙRest
    return n_mean
end

function measure_energy(
    self::Lattice, Uup_array::Array{Float64,3}, Udn_array::Array{Float64,3},
    nup_mean::Array{Float64,1}, ndn_mean::Array{Float64,1}
)
    EOneBody, ETwoBody = 0, 0

    for i = 1 : self.nsites
        for j = 1 : self.nsites
            if abs(i - j) == 1 || abs(i - j) == self.nsites - 1
                Dij_up = sum(Uup_array[:, i, j] .* nup_mean)
                Dij_dn = sum(Udn_array[:, i, j] .* ndn_mean)
                EOneBody += -self.t * ( Dij_up + Dij_dn )
            end
        end
    end

    for i = 1 : self.nsites
        Umatrix = Uup_array[:, i, i] * Udn_array[:, i, i]'
        nmatrix = nup_mean * ndn_mean'
        Diiii = sum(Umatrix .* nmatrix)
        ETwoBody += self.U * Diiii
    end

    return EOneBody, ETwoBody
end

# determine the sign of the accept_ratio
@inline function Rsign(R::Float64)
    Rsgn = 2 * (R > 0) - 1
end

function MCIntegration(lattice::Lattice, qmc::QMC, opt::Opt)

    # initialize vectors that store MC data
    Ek_array = Vector{Float64}() # kinetic energy
    Ep_array = Vector{Float64}() # potential energy
    E_array = Vector{Float64}() # total energy
    sgn_array = Vector{Int64}() # sign

    # initialize a random configuration
    σlist = 2 * (rand(Float64, (lattice.nsites, qmc.L)) .< 0.5) .- 1
    EigensUp, EigensDn = matrix_product_QR(σlist, lattice, qmc, opt.stab_interval)
    # Spin Section (with truncation)
    nup_mean, Zₙupratio, Zₙupsgn = SpectrumRemodification(lattice.nsites, lattice.nup, real(EigensUp.values))
    ndn_mean, Zₙdnratio, Zₙdnsgn = SpectrumRemodification(lattice.nsites, lattice.ndn, real(EigensDn.values))
    sgn = Zₙupsgn * Zₙdnsgn

    for n = 1 : qmc.samples
        σlist_new = copy(σlist)

        # random flip a spin
        flip = rand(1 : qmc.L * lattice.nsites)
        σlist_new[flip] *= -1
        EigensUp_new, EigensDn_new = matrix_product_QR(σlist_new, lattice, qmc, opt.stab_interval)
        # New Spin Section (with truncation)
        nup_mean_new, Zₙupratio_new, Zₙupsgn_new = SpectrumRemodification(lattice.nsites, lattice.nup, real(EigensUp_new.values))
        ndn_mean_new, Zₙdnratio_new, Zₙdnsgn_new = SpectrumRemodification(lattice.nsites, lattice.ndn, real(EigensDn_new.values))
        # test
        if nup_mean_new[lattice.nsites] > 1.1
            println(EigensUp.values)
        elseif ndn_mean_new[lattice.nsites] > 1.1
            println(EigensDn.values)
        end

        Rup = prod(Zₙupratio_new ./ Zₙupratio)
        Rdn = prod(Zₙdnratio_new ./ Zₙdnratio)
        sgn_new = Zₙupsgn_new * Zₙdnsgn_new
        R = abs(Rup * Rdn)
        accept_ratio = R / (1 + R)

        if rand() <= accept_ratio   # accept

            σlist = σlist_new
            Zₙupratio, Zₙdnratio = Zₙupratio_new, Zₙdnratio_new
            EigensUp, EigensDn = EigensUp_new, EigensDn_new
            nup_mean, ndn_mean = nup_mean_new, ndn_mean_new
            sgn = sgn_new

            #Zₙup_true = RecursionPF(lattice.nup, real(EigensUp.values))
            #Zₙdn_true = RecursionPF(lattice.ndn, real(EigensDn.values))
            #Zₙ_true = Zₙup_true[lattice.nup] * Zₙdn_true[lattice.ndn]
            #Zₙ = prod(Zₙup_ratio[1 : lattice.nup]) * prod(Zₙdn_ratio[1 : lattice.ndn])
            #if Zₙ / Zₙ_true > 1e4
                #println(EigensUp.values)
                #println(EigensDn.values)
                #break
            #end

            Uup_array, Udn_array = overlap_coefficient(lattice.nsites, EigensUp.vectors, EigensDn.vectors)
            Ek, Ep = measure_energy(lattice, Uup_array, Udn_array, nup_mean, ndn_mean)
            push!(sgn_array, sgn)
            push!(Ek_array, Ek * sgn)
            push!(Ep_array, Ep * sgn)
            push!(E_array, Ek * sgn + Ep * sgn)

        else    # reject
            #Zₙup_true = RecursionPF(lattice.nup, real(EigensUp.values))
            #Zₙdn_true = RecursionPF(lattice.ndn, real(EigensDn.values))
            #Zₙ_true = Zₙup_true[lattice.nup] * Zₙdn_true[lattice.ndn]
            #Zₙ = prod(Zₙup_ratio[1 : lattice.nup]) * prod(Zₙdn_ratio[1 : lattice.ndn])

            Uup_array, Udn_array = overlap_coefficient(lattice.nsites, EigensUp.vectors, EigensDn.vectors)
            Ek, Ep = measure_energy(lattice, Uup_array, Udn_array, nup_mean, ndn_mean)
            push!(sgn_array, sgn)
            push!(Ek_array, Ek * sgn)
            push!(Ep_array, Ep * sgn)
            push!(E_array, Ek * sgn + Ep * sgn)

        end
    end

    E_mean = mean(E_array[opt.nburn_in : qmc.samples])
    E_error = std(E_array[opt.nburn_in : qmc.samples]) / sqrt(qmc.samples - opt.nburn_in)
    Ek_mean = mean(Ek_array[opt.nburn_in : qmc.samples])
    Ek_error = std(Ek_array[opt.nburn_in : qmc.samples]) / sqrt(qmc.samples - opt.nburn_in)
    Ep_mean = mean(Ep_array[opt.nburn_in : qmc.samples])
    Ep_error = std(Ep_array[opt.nburn_in : qmc.samples]) / sqrt(qmc.samples - opt.nburn_in)
    sgn_avg = mean(sgn_array[opt.nburn_in : qmc.samples])

    return E_mean, E_error, Ek_mean, Ek_error, Ep_mean, Ep_error, sgn_avg
end
