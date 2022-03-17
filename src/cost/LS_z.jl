

######### z, compensation for α at the (i,k) level. (spin system index, partition element index).
# include singlets too.


function updatez!(  A::Matrix{T},
    b::Vector{T},
    z::Vector{Complex{T}},
    U_LS,
    Gs::Vector{NMRSpectraSimulator.zCompoundFIDType{T,SST}},
    w::Vector{T}) where {T <: Real,SST}

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_LS)

    evaldesignmatrixz!(A, U_LS, Gs, w)

    z[:] = A\b

    return nothing
end

function evaldesignmatrixz!(B::Matrix{T},
    U,
    Gs::Vector{NMRSpectraSimulator.zCompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType1{T}}},
    w::Vector{T}) where T <: Real

    #
    M = length(U)
    N = length(Gs)

    N_z, N_z_singlets = countz(Gs)
    @assert size(B) == (2*M, N_z + N_z_singlets)
    fill!(B, Inf) # debug.

    resetz!(Gs)
    j = 0

    # loop over each κ partition element in Gs.
    for n = 1:N
        A = Gs[n]
        @assert length(A.zs) == length(A.core.qs)

        # spin system.
        for i = 1:length(A.zs)

            # partition
            for k = 1:length(A.zs[i])

                j += 1

                # loop over each fit position.
                for m = 1:M

                    # taken from evalitproxysys()
                    r = 2*π*U[m] - A.core.ss_params.d[i]
                    out = w[n]*A.core.gs[i][k](r, A.core.ss_params.κs_λ[i])

                    B[m,j], B[m+M,j] = real(out), imag(out)
                end
            end
        end

        # singlets.
        for k = 1:length(A.zs_singlets)
            j += 1

            for m = 1:M

                tmp = w[n]*evalκsingletsz(U[m], A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.λ0, A.core.κs_λ_singlets)

                B[m,j] = real(tmp)
                B[m+M,j] = imag(tmp)
            end
        end

    end

    return nothing
end

# modified from NMRSpectraSimulator.evalκsinglets().
function evalκsingletsz(u::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    λ0::T, λ_multipliers::Vector{T})::Complex{T} where T <: Real

    u_rad = 2*π*u

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        #out += αs_singlets[i]*exp(im*βs_singlets[i])/(λ+im*(τ-Ω))
        out += αs_singlets[i]/(λ+im*(τ-Ω))
    end
    return out
end
