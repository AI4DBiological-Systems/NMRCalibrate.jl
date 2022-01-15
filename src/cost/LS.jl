

# just use the real part for LS.
#yr, yi assumed to have LS_inds applied already.
function solveBLScL( A,
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
    Es::Vector{NMRSpectraSimulator.CompoundFIDType{T}},
    w_lower::Vector{T},
    w_upper::Vector{T}) where T <: Real

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
    Es::Vector{NMRSpectraSimulator.CompoundFIDType{T}}) where T <: Real

    #
    M = length(U)
    N = length(Es)

    @assert size(B) == (2*M,N)

    for n = 1:N
        for m = 1:M
            # speed up later.
            tmp = NMRSpectraSimulator.evalitpproxycompound(U[m], Es[n])
            B[m,n], B[m+M,n] = real(tmp), imag(tmp)
        end
    end

    return nothing
end

######### κ, compensation for α.


function updateκ!(  A::Matrix{T},
    b::Vector{T},
    κ::Vector{T},
    U_LS,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T}},
    κ_lower::Vector{T},
    κ_upper::Vector{T}) where T <: Real

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_LS)

    evaldesignmatrixκ!(A, U_LS, Es)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    κ[:], status_flag = solveBLScL(A, b, κ_lower, κ_upper)

    return status_flag
end


function evaldesignmatrixκ!(B::Matrix{T},
    U,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T}}) where T <: Real

    #
    M = length(U)
    N = length(Es)

    N_κ, N_κ_singlets = NMRCalibrate.countκ(Es)
    #println((N_κ, N_κ_singlets))
    #println(size(B))
    @assert size(B) == (2*M, N_κ + N_κ_singlets)

    resetκ!(Es)
    j = 0

    # loop over each κ partition element in Es.
    for n = 1:N
        A = Es[n]
        @assert length(A.κ) == length(A.core.qs) == length(A.core.κs_λ) == length(A.core.κs_β) == length(A.core.d)
        
         # spin system.
        for i = 1:length(A.κ)
            
            # partition
            for k = 1:length(A.κ[i])
                
                j += 1
                
                # loop over each fit position.
                for m = 1:M

                    # taken from evalitproxysys()
                    r = 2*π*U[m] - A.core.d[i]
                    out = A.core.qs[i][k](r, A.core.κs_λ[i], A.core.κs_β[i])

                    #tmp = NMRSpectraSimulator.evalitpproxycompound(U[m], A)
                    # tmp = one κ partition.

                    B[m,j], B[m+M,j] = real(out), imag(out)
                end
            end
        end

        # singlets.
        for k = 1:length(A.κ_singlets)
            j += 1

            for m = 1:M

                tmp = NMRSpectraSimulator.evalκsinglets(U[m], A.κ_singlets, A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets)

                B[m,j] += real(tmp)
                B[m+M,j] += imag(tmp)
            end
        end

    end

    return nothing
end