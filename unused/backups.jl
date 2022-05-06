
function getmeanΔc(Δc_m_compound::Vector{Vector{Vector{T}}},
    part_inds_compound::Vector{Vector{Vector{Int}}}) where T

    @assert length(part_inds_compound) == length(Δc_m_compound)

    out = Vector{Vector{Vector{T}}}(undef, length(part_inds_compound))

    for i = 1:length(part_inds_compound)

        N_parition_elements = length(part_inds_compound[i])
        out[i] = Vector{Vector{T}}(undef, N_parition_elements)

        for k = 1:N_parition_elements
            inds = part_inds_compound[i][k]

            out[i][k] = Statistics.mean(Δc_m_compound[i][inds])
        end
    end

    return out
end

A = As[1]
Δc_m_compound_avg = getmeanΔc(A.Δc_m_compound, A.part_inds_compound)
#CNMRSpectraSimulator.combinevectors(Δc_m_compound_avg)

function estimateβ(inds_set::Vector{Vector{Int}},
    c_all::Vector{Vector{T}},
    phase_vec::Vector{T};
    st_ind::Int = 1) where T
    # loop through As, assemble. Δc matrix.
    #N_rows = sum( sum( length(As[n].part_inds_compound[i]) for i = 1:length(As[n].part_inds_compound) ) for n = 1:length(As) )

    @assert !isempty(c_all)
    N_partition_elements = length(inds_set)
    N_spins = length(c_all[1])

    C = Matrix{T}(undef, N_partition_elements, N_spins)
    y = Vector{T}(undef, N_partition_elements)

    for k = 1:N_partition_elements
        inds = inds_set[k]

        C[k,:] = Statistics.mean( c_all[inds] )
        y[k] = phase_vec[st_ind+k-1]
    end

    p_β = C\y

    fin_ind = st_ind + N_spins - 1

    return p_β, fin_ind, C, y
end

# z_flatten is the sum of all partition elements, of all spin groups in mixture.
# p_β is the sum of all β parameters, which is sum of all spins in each spin group in mixture.
function getinitialβ!(p_β::Vector{T},
    z_flatten::Vector{Complex{T}}) where T <: Real

    fill!(p_β, Inf)
    st_ind = 1

    for n = 1:length(As)
        A = As[n]

        for i = 1:length(A.part_inds_compound)
            p_β_i, fin_ind, C, y = estimateβ(A.part_inds_compound[i],
                A.Δc_m_compound[i], angle.(z_flatten); st_ind = st_ind)

            println("norm(C*p_β_i-y) = ", norm(C*p_β_i-y)) # debug.
            #@assert length(p_β_i) == length(st_ind:fin_ind)
            p_β[st_ind:fin_ind] = p_β_i
            st_ind = fin_ind + 1
        end
    end

    return nothing
end

# β_z = similar(β_initial)
# fill!(β_z, Inf)
# getinitialβ!(β_z, z_LS)
# #display([β_z; ])
# clamp!(β_z, -π+1e-5, π -1e-5)



# just use the real part for LS.
#yr, yi assumed to have LS_inds applied already.
# This does A\b but with bounds.
function solveBLScL( A::Matrix{T},
    b,
    x_lower::Vector{T},
    x_upper::Vector{T}) where T <: Real

    # This is the formula.
    # BoundedLeastSquares.Quadratic(A'*A, A'*b)

    AtA = forcesymmetric(A'*A)
    Aty = A'*b

    Q = BoundedLeastSquares.Quadratic(AtA, Aty)

    x_opt, out_flag = BoundedLeastSquares.min_bound_constrained_quadratic(Q,
                        x_lower, x_upper)

    return x_opt, out_flag
end

####### w.

function setupupdatew(N_positions::Int, N_compounds::Int)

    T = Float64

    # q_t.
    A = Matrix{T}(undef, 2*N_positions, N_compounds)
    w = zeros(T, N_compounds)

    return A, w
end

function updatew!(  A::Matrix{T},
    b::Vector{T},
    w::Vector{T},
    U_LS,
    Es::Vector{NMRSpectraSimulator.FIDModelType{T,SST}},
    w_lower::Vector{T},
    w_upper::Vector{T}) where {T <: Real, SST}

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_LS)

    evaldesignmatrixw!(A, U_LS, Es)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    w[:], status_flag = solveBLScL(A, b, w_lower, w_upper)

    return status_flag
end


function evaldesignmatrixw!(B::Matrix{T},
    U,
    Es::Vector{NMRSpectraSimulator.FIDModelType{T,SST}}) where {T <: Real, SST}

    #
    M = length(U)
    N = length(Es)

    @assert size(B) == (2*M,N)

    fill!(B, NaN) # debug.

    for n = 1:N
        for m = 1:M
            # speed up later.
            tmp = NMRSpectraSimulator.evalitpproxycompound(U[m], Es[n])
            B[m,n], B[m+M,n] = real(tmp), imag(tmp)
        end
    end

    return nothing
end


function updateκ!(  A::Matrix{T},
    #r_buffer::Vector{T},
    b::Vector{T},
    κ::Vector{T},
    U_LS,
    Es::Vector{NMRSpectraSimulator.καFIDModelType{T,SST}},
    w::Vector{T},
    κ_lower::Vector{T},
    κ_upper::Vector{T}) where {T <: Real,SST}

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_LS)

    r_buffer = Vector{T}(undef, 0)
    evaldesignmatrixκ!(A, r_buffer, U_LS, Es, w)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    # if !all(isfinite.(A))
    #     println("A is not finite")

    #     JLD.save("debug2.jld", "A", A, "b", b)

    #     # f = jldopen(filename, "r+")
    #     # write(f, "A", A, "b", b)
    #     # close(f)
    # end

    κ[:], status_flag = solveBLScL(A, b, κ_lower, κ_upper)
    # κ[:] = A\b
    # clamp!(κ, κ_lower[1], κ_upper[1])
    # status_flag = true

    return status_flag
end

function evaldesignmatrixκ!(B::Matrix{T},
    U,
    Es::Vector{NMRSpectraSimulator.καFIDModelType{T,NMRSpectraSimulator.SpinSysParamsType1{T}}},
    w::Vector{T}) where T <: Real

    #
    M = length(U)
    N = length(Es)

    N_κ, N_κ_singlets = countκ(Es)
    #println((N_κ, N_κ_singlets))
    #println(size(B))
    @assert size(B) == (2*M, N_κ + N_κ_singlets)
    fill!(B, Inf) # debug.

    resetκ!(Es)
    j = 0

    #U_rad = U .* (2*π)
    #resize!(r_buffer, M)

    # loop over each κ partition element in Es.
    for n = 1:N
        A = Es[n]
        @assert length(A.κs_α) == length(A.core.qs)

        # spin system.
        for i = 1:length(A.κs_α)

            # # setup.
            # for m = 1:M
            #     r_buffer[m] = 2*π*U[m] - A.core.ss_params.d[i]
            # end

            # partition
            for k = 1:length(A.κs_α[i])

                j += 1

                # loop over each fit position.
                for m = 1:M

                    # taken from evalitproxysys()
                    r = 2*π*U[m] - A.core.ss_params.d[i]
                    #r = U_rad[m] - A.core.ss_params.d[i]
                    #r = r_buffer[m]

                    out = w[n]*A.core.qs[i][k](r, A.core.ss_params.κs_λ[i])

                    B[m,j], B[m+M,j] = real(out), imag(out)
                end
            end
        end

        # singlets.
        for k = 1:length(A.κs_α_singlets)
            j += 1

            for m = 1:M

                tmp = w[n]*NMRSpectraSimulator.evalκsinglets(U[m], A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets)

                B[m,j] = real(tmp)
                B[m+M,j] = imag(tmp)
            end
        end

    end

    return nothing
end

function evaldesignmatrixκ!(B::Matrix{T}, r_buffer_unused::Vector{T},
    U,
    Es::Vector{NMRSpectraSimulator.καFIDModelType{T,NMRSpectraSimulator.SpinSysParamsType2{T}}},
    w::Vector{T}) where T <: Real

    #
    M = length(U)
    N = length(Es)

    N_κ, N_κ_singlets = countκ(Es)
    #println((N_κ, N_κ_singlets))
    #println(size(B))
    @assert size(B) == (2*M, N_κ + N_κ_singlets)
    fill!(B, Inf) # debug.

    resetκ!(Es)
    j = 0

    #resize!(r_buffer, M)

    # loop over each κ partition element in Es.
    for n = 1:N
        A = Es[n]
        @assert length(A.κs_α) == length(A.core.qs)

        # spin system.
        for i = 1:length(A.κs_α)

            # partition
            for k = 1:length(A.κs_α[i])

                j += 1

                # loop over each fit position.
                for m = 1:M

                    # taken from evalitproxysys()
                    r = 2*π*U[m] - A.core.ss_params.d[i][k]
                    out = w[n]*A.core.qs[i][k](r, A.core.ss_params.κs_λ[i][k])

                    B[m,j], B[m+M,j] = real(out), imag(out)
                end
            end
        end

        # singlets.
        for k = 1:length(A.κs_α_singlets)
            j += 1

            for m = 1:M

                tmp = w[n]*NMRSpectraSimulator.evalκsinglets(U[m], A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets)

                B[m,j] = real(tmp)
                B[m+M,j] = imag(tmp)
            end
        end

    end

    return nothing
end
