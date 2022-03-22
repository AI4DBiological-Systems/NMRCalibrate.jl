

# just use the real part for LS.
#yr, yi assumed to have LS_inds applied already.
# This does A\b but with bounds.
function solveBLScL( W::Matrix{Complex{T}},
    b,
    x_lower::Vector{T},
    x_upper::Vector{T}) where T <: Real

    # This is the formula.
    # BoundedLeastSquares.Quadratic(A'*A, A'*b)

    A = reinterpret(T, W)

    AtA = forcesymmetric(A'*A)
    Aty = A'*b

    Q = BoundedLeastSquares.Quadratic(AtA, Aty)

    x_opt, out_flag = BoundedLeastSquares.min_bound_constrained_quadratic(Q,
                        x_lower, x_upper)

    return x_opt, out_flag
end

####### w.

function setupupdatew(N_positions::Int, N_compounds::Int, y::Vector{Complex{T}}) where T <: Real

    A = Matrix{Complex{T}}(undef, N_positions, N_compounds)
    w = zeros(T, N_compounds)
    b = collect(reinterpret(T, y))

    return A, w, b
end

"""
A is an M x L complex-valued matrix, b is 2*M 1D array.
"""
function updatew!(  A::Matrix{Complex{T}},
    b,
    w::Vector{T},
    U_rad_LS,
    Es::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}},
    w_lower::Vector{T},
    w_upper::Vector{T}) where {T <: Real, SST}

    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrixw!(A, U_rad_LS, Es)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    w[:], status_flag = solveBLScL(A, b, w_lower, w_upper)

    return status_flag
end


function evaldesignmatrixw!(B::Matrix{Complex{T}},
    U_rad,
    Es::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}}) where {T <: Real, SST}

    #
    M = length(U_rad)
    N = length(Es)

    @assert size(B) == (M,N)

    fill!(B, NaN) # debug.

    for n = 1:N
        for m = 1:M
            # speed up later.
            B[m,n] = NMRSpectraSimulator.evalitpproxycompound(U_rad[m], Es[n])
        end
    end

    return nothing
end

######### κ, compensation for α.


function updateκ!(  A::Matrix{Complex{T}},
    #r_buffer::Vector{T},
    b,
    κ::Vector{T},
    U_rad_LS,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}},
    w::Vector{T},
    κ_lower::Vector{T},
    κ_upper::Vector{T}) where {T <: Real,SST}

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrixκ!(A, U_rad_LS, Es, w)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    κ[:], status_flag = solveBLScL(A, b, κ_lower, κ_upper)
    # κ[:] = A\b
    # clamp!(κ, κ_lower[1], κ_upper[1])
    # status_flag = true

    return status_flag
end





function evaldesignmatrixκ!(B::Matrix{Complex{T}},
    U_rad,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType2{T}}},
    w::Vector{T}) where T <: Real

    #
    M = length(U_rad)
    N = length(Es)

    N_κ, N_κ_singlets = countκ(Es)
    #println((N_κ, N_κ_singlets))
    #println(size(B))
    @assert size(B) == (M, N_κ + N_κ_singlets)
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
                    r = U_rad[m] - A.core.ss_params.d[i][k]
                    B[m,j] = w[n]*A.core.qs[i][k](r, A.core.ss_params.κs_λ[i][k])
                end
            end
        end

        # singlets.
        for k = 1:length(A.κs_α_singlets)
            j += 1

            for m = 1:M

                B[m,j] = w[n]*NMRSpectraSimulator.evalκsinglets(U_rad[m], A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets)
            end
        end

    end

    return nothing
end
