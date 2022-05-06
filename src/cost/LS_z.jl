#### z parsing.

function countz(Gs::Vector{NMRSpectraSimulator.zFIDModelType{T,SST}}) where {T,SST}

    N_z = 0
    N_z_singlets = 0
    for n = 1:length(Gs)
        # for i = 1:length(Gs[n].zs)
        #     #for l = 1:length(Gs[n].zs[i])
        #     N_z += length(Gs[n].zs[i])
        #
        #     #end
        # end
        N_z += sum(length.(Gs[n].zs))

        N_z_singlets += length(Gs[n].zs_singlets)
    end

    return N_z, N_z_singlets
end

function parsez!(Gs::Vector{NMRSpectraSimulator.zFIDModelType{T,SST}},
    z_vec::Vector{Complex{T}}) where {T,SST}

    @assert length(z_vec) == sum(countz(Gs))

    j = 0
    for n = 1:length(Gs)
        for i = 1:length(Gs[n].zs)
            for l = 1:length(Gs[n].zs[i])
                j += 1
                Gs[n].zs[i][l] = z_vec[j]
            end
        end

        for i = 1:length(Gs[n].zs_singlets)
            j += 1
            Gs[n].zs_singlets[i] = z_vec[j]
        end
    end

    return j
end

function resetz!(Gs::Vector{NMRSpectraSimulator.zFIDModelType{T,SST}}) where {T,SST}

    j = 0
    for n = 1:length(Gs)
        for i = 1:length(Gs[n].zs)
            for l = 1:length(Gs[n].zs[i])
                j += 1
                Gs[n].zs[i][l] = one(T)
            end
        end

        for i = 1:length(Gs[n].zs_singlets)
            j += 1
            Gs[n].zs_singlets[i] = one(T)
        end
    end

    return j
end


######### z, compensation for α at the (i,k) level. (spin system index, partition element index).
# include singlets too.


function updatez!(  A::Matrix{Complex{T}},
    b::Vector{T},
    z::Vector{Complex{T}},
    U_rad_LS,
    Gs::Vector{NMRSpectraSimulator.zFIDModelType{T,SST}},
    w::Vector{T}) where {T <: Real,SST}

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrixz!(A, U_rad_LS, Gs, w)

    W = reinterpret(T, A)
    z[:] = W\b

    return nothing
end

function evaldesignmatrixz!(B::Matrix{Complex{T}},
    U_rad,
    Es::Vector{NMRSpectraSimulator.zFIDModelType{T,NMRSpectraSimulator.SpinSysParamsType1{T}}},
    w::Vector{T}) where T <: Real

    #
    M = length(U_rad)
    N = length(Es)

    #N_κ, N_κ_singlets = countκ(Es)
    #@assert size(B) == (M, N_κ + N_κ_singlets)
    #fill!(B, Inf) # debug.

    #resetκ!(Es)
    j = 0

    # loop over each κ partition element in Es.
    for n = 1:N
        A = Es[n]
        @assert length(A.zs) == length(A.core.qs)

        #spin system.
        #for i = 1:length(A.zs)
        for i in eachindex(A.zs)

            # partition
            for k in eachindex(A.zs[i])
                j += 1

                B[:,j] = collect( w[n]*A.core.gs[i][k](U_rad[m] - A.core.ss_params.d[i], A.core.ss_params.κs_λ[i]) for m = 1:M )

            end
        end

        # singlets.
        for k = 1:length(A.zs_singlets)
            j += 1

            B[:,j] = collect( w[n]*evalκsingletsz(U_rad[m], A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets) for m = 1:M )
        end

    end

    return nothing
end

# modified from NMRSpectraSimulator.evalκsinglets().
function evalκsingletsz(u_rad::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    λ0::T, λ_multipliers::Vector{T})::Complex{T} where T <: Real

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
