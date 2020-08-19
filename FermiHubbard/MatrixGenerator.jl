using LinearAlgebra

# one-body matrix, nearest neighbour, PBC
function kinetic_matrix(self::Lattice)
    K_matrix = zeros(self.nsites, self.nsites)
    for i =  1 : self.nsites
        for j = 1 : self.nsites
            if abs(i - j) == 1 || abs(i - j) == self.nsites - 1
                K_matrix[i, j] = -self.t
            end
        end
    end
    return K_matrix
end

# σ-dependent generator
function auxiliary_field_matrix_stepwise(σ::Array{Int64,1}, γ::Float64, Δτ::Float64, U::Float64)
    AF_matrix_SpinUp = diagm(0 => exp.(2 * γ * σ .- Δτ * U / 2))
    AF_matrix_SpinDn = diagm(0 => exp.(-2 * γ * σ .- Δτ * U / 2))
    return AF_matrix_SpinUp, AF_matrix_SpinDn
end

function matrix_product_stepwise(σlist::Array{Int64,2}, lattice::Lattice, qmc::QMC)
    MatProdUp = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    MatProdDn = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    for i = 1 : qmc.L
        σ = σlist[:, i]
        AF_matrix_SpinUp, AF_matrix_SpinDn = auxiliary_field_matrix_stepwise(σ, qmc.γ, qmc.Δτ, lattice.U)
        MatProdUp = Symmetric(qmc.Kmatrix * AF_matrix_SpinUp * qmc.Kmatrix) * MatProdUp
        MatProdDn = Symmetric(qmc.Kmatrix * AF_matrix_SpinDn * qmc.Kmatrix) * MatProdDn
    end
    EigensUp = eigen(MatProdUp)
    EigensDn = eigen(MatProdDn)
    return MatProdUp, MatProdDn, EigensUp, EigensDn
end

# Use QR decomposition to stablize the propogation
function matrix_product_QR(σlist::Array{Int64,2}, lattice::Lattice, qmc::QMC, n_stab_interval::Int64)
    # Initialize QR matrices
    QRUp_Q = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    QRUp_R = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    QRDn_Q = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    QRDn_R = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    # Final Product Matrix
    MatProdUp = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    MatProdDn = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
    for i = 1 : div(qmc.L, n_stab_interval)
        matprodUp = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
        matprodDn = Matrix{Float64}(I, lattice.nsites, lattice.nsites)
        for j = 1 : n_stab_interval
            σ = σlist[:, (i - 1) * n_stab_interval + j]
            AF_matrix_SpinUp, AF_matrix_SpinDn = auxiliary_field_matrix_stepwise(σ, qmc.γ, qmc.Δτ, lattice.U)
            matprodUp = Symmetric(qmc.Kmatrix * AF_matrix_SpinUp * qmc.Kmatrix) * matprodUp
            matprodDn = Symmetric(qmc.Kmatrix * AF_matrix_SpinDn * qmc.Kmatrix) * matprodDn
        end
        # QR decomposition of the matrix product
        # leave here to see if can be replaced with rmul!
        QRUp_Q = matprodUp * QRUp_Q
        DecompUp = qr(QRUp_Q)
        QRDn_Q = matprodDn * QRDn_Q
        DecompDn = qr(QRDn_Q)
        # Remaining R_i goes to multiply R_1
        QRUp_Q = DecompUp.Q
        QRUp_R = DecompUp.R * QRUp_R
        QRDn_Q = DecompDn.Q
        QRDn_R = DecompDn.R * QRDn_R
    end
    mul!(MatProdUp, QRUp_Q, QRUp_R)
    mul!(MatProdDn, QRDn_Q, QRDn_R)
    EigensUp = eigen(MatProdUp)
    EigensDn = eigen(MatProdDn)
    return MatProdUp, MatProdDn, EigensUp, EigensDn
end

# generate overlap "tensor" / column major
function overlap_coefficient(nsites::Int64, Uup::Array{T1,2}, Udn::Array{T2,2}) where {T1<:FloatType, T2<:FloatType}
    Uup_inv, Udn_inv = inv(Uup), inv(Udn)
    Uup_array = complex(zeros(nsites, nsites, nsites))
    Udn_array = complex(zeros(nsites, nsites, nsites))
    for i = 1 : nsites
        for j = 1 : nsites
            for k = 1 : nsites
                Uup_array[k ,i, j] = Uup_inv[k, i] * Uup[j, k]
                Udn_array[k, i, j] = Udn_inv[k, i] * Udn[j, k]
            end
        end
    end
    real(Uup_array), real(Udn_array)
end
