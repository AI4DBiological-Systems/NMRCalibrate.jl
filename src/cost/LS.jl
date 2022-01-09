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

    evaldesignmatrix!(A, U_LS, Es)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    w[:], status_flag = solveBLScL(A, b, w_lower, w_upper)

    return status_flag
end

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

function evaldesignmatrix!(B::Matrix{T},
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