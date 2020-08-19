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


function recursion(self::Lattice, expϵup::Array{Float64,1}, expϵdn::Array{Float64,1})
    nup_mean = zeros(Float64, self.nsites, self.nsites)
    ndn_mean = zeros(Float64, self.nsites, self.nsites)
    Znup_ratio = zeros(Float64, self.ntot)
    Zndn_ratio = zeros(Float64, self.ntot)
    expϵup = expϵup ./ expϵup[self.nsites]
    expϵdn = expϵdn ./ expϵdn[self.nsites]
    Znup_ratio[1], Zndn_ratio[1] = sum(expϵup), sum(expϵdn)
    nup_mean[:, 1], ndn_mean[:, 1] = expϵup / Znup_ratio[1], expϵdn / Zndn_ratio[1]
    for N = 2 : max(self.nup, self.ndn)
        Znup_ratio[N] = sum(expϵup .* (1 .- nup_mean[:, N-1])) / N
        Zndn_ratio[N] = sum(expϵdn .* (1 .- ndn_mean[:, N-1])) / N
        nup_mean[:, N] = expϵup .* (1 .- nup_mean[:, N - 1]) / Znup_ratio[N]
        ndn_mean[:, N] = expϵdn .* (1 .- ndn_mean[:, N - 1]) / Zndn_ratio[N]
    end
    return nup_mean, ndn_mean, Znup_ratio, Zndn_ratio
end

function measure_energy(
    self::Lattice, Uup_array::Array{Float64,3}, Udn_array::Array{Float64,3},
    nup_mean::Array{Float64,2}, ndn_mean::Array{Float64,2}
)
    Eonebody, Etwobody = 0, 0
    for i = 1 : self.nsites
        for j = 1 : self.nsites
            if abs(i - j) == 1 || abs(i - j) == self.nsites - 1
                Dij_up = sum(Uup_array[:, i, j] .* nup_mean[:, self.nup])
                Dij_dn = sum(Udn_array[:, i, j] .* ndn_mean[:, self.ndn])
                Eonebody += -self.t * ( Dij_up + Dij_dn )
            end
        end
    end
    for i = 1 : self.nsites
        Umatrix = Uup_array[:, i, i] * Udn_array[:, i, i]'
        nmatrix = nup_mean[:, self.nup] * ndn_mean[:, self.ndn]'
        Diiii = sum(Umatrix .* nmatrix)
        Etwobody += self.U * Diiii
    end
    return Eonebody, Etwobody
end

# rank-one perturbation calculator
function rankone_perturb!(
    γ::Float64, σil::Int64, Ns::Int64,
    MatProdUp_new::Array{Float64,2}, MatProdDn_new::Array{Float64,2}
)
    Δup, Δdn = Matrix{Float64}(I, size(MatProdUp_new)), Matrix{Float64}(I, size(MatProdDn_new))
    Δup[Ns, Ns] = exp(-4 * γ * σil)
    Δdn[Ns, Ns] = exp(4 * γ * σil)
    MatProdUp_new = Δup * MatProdUp_new
    MatProdDn_new = Δdn * MatProdDn_new
    return MatProdUp_new, MatProdDn_new
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
    NiMax_array = Vector{Float64}()
    # initialize a random configuration
    σlist = 2 * (rand(Float64, (lattice.nsites, qmc.L)) .< 0.5) .- 1
    MatProdUp, MatProdDn, EigensUp, EigensDn = matrix_product_QR(σlist, lattice, qmc, opt.stab_interval)
    nup_mean, ndn_mean, Znup_ratio, Zndn_ratio = recursion(lattice, real(EigensUp.values), real(EigensDn.values))
    sgn = Rsign(prod(Znup_ratio[1 : lattice.nup]) * prod(Zndn_ratio[1 : lattice.nup]))
    for n = 1 : qmc.samples
        σlist_new = copy(σlist)
        # random flip a spin
        #flipL, flipNs = rand(1 : qmc.L), rand(1 : lattice.nsites)
        flip = rand(1 : qmc.L * lattice.nsites)
        σlist_new[flip] *= -1
        #rankone_perturb!(qmc.γ, σlist[flipNs, flipL], flipNs, MatProdUp_new, MatProdDn_new)
        MatProdUp_new, MatProdDn_new, EigensUp_new, EigensDn_new = matrix_product_QR(σlist_new, lattice, qmc, opt.stab_interval)
        nup_mean_new, ndn_mean_new, Znup_ratio_new, Zndn_ratio_new = recursion(lattice, real(EigensUp_new.values), real(EigensDn_new.values))
        Rup = prod(Znup_ratio_new[1 : lattice.nup] ./ Znup_ratio[1 : lattice.nup])
        Rdn = prod(Zndn_ratio_new[1 : lattice.ndn] ./ Zndn_ratio[1 : lattice.ndn])
        sgn_new = Rsign(Rup * Rdn) * sgn
        R = abs(Rup * Rdn)
        accept_ratio = R / (1 + R)
        if rand() <= accept_ratio   # accept
            σlist = σlist_new
            Znup_ratio, Zndn_ratio = Znup_ratio_new, Zndn_ratio_new
            EigensUp, EigensDn = EigensUp_new, EigensDn_new
            nup_mean, ndn_mean = nup_mean_new, ndn_mean_new
            sgn = sgn_new
            Uup_array, Udn_array = overlap_coefficient(lattice.nsites, EigensUp.vectors, EigensDn.vectors)
            Eonebody, Etwobody = measure_energy(lattice, Uup_array, Udn_array, nup_mean, ndn_mean)
            push!(sgn_array, sgn)
            push!(Ek_array, Eonebody * sgn)
            push!(Ep_array, Etwobody * sgn)
            push!(E_array, Eonebody * sgn + Etwobody * sgn)
            push!(NiMax_array, maximum(nup_mean[:, lattice.nup]))
        else    # reject
            Uup_array, Udn_array = overlap_coefficient(lattice.nsites, EigensUp.vectors, EigensDn.vectors)
            Eonebody, Etwobody = measure_energy(lattice, Uup_array, Udn_array, nup_mean, ndn_mean)
            push!(sgn_array, sgn)
            push!(Ek_array, Eonebody * sgn)
            push!(Ep_array, Etwobody * sgn)
            push!(E_array, Eonebody * sgn + Etwobody * sgn)
            push!(NiMax_array, maximum(nup_mean[:, lattice.nup]))
        end
    end
    E_mean = mean(E_array[opt.nburn_in : qmc.samples])
    E_error = std(E_array[opt.nburn_in : qmc.samples]) / sqrt(qmc.samples - opt.nburn_in)
    Ek_mean = mean(Ek_array[opt.nburn_in : qmc.samples])
    Ek_error = std(Ek_array[opt.nburn_in : qmc.samples]) / sqrt(qmc.samples - opt.nburn_in)
    Ep_mean = mean(Ep_array[opt.nburn_in : qmc.samples])
    Ep_error = std(Ep_array[opt.nburn_in : qmc.samples]) / sqrt(qmc.samples - opt.nburn_in)
    sgn_avg = mean(sgn_array[opt.nburn_in : qmc.samples])
    NiMax_avg = mean(NiMax_array[opt.nburn_in : qmc.samples])
    return E_mean, E_error, Ek_mean, Ek_error, Ep_mean, Ep_error, sgn_avg, NiMax_avg
end
