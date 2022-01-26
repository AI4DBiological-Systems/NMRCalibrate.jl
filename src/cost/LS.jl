

# just use the real part for LS.
#yr, yi assumed to have LS_inds applied already.
# This does A\b but with bounds.
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
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T}}) where T <: Real

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

                    # if !isfinite(out)
                    #     println("eval not finite!")
                    #     println("U[m] = ", U[m])

                    #     println("n,ik,m = ", (n,i,k,m))
                    #     println("Es[n].core.d = ", Es[n].core.d)
                    #     println("Es[n].core.κs_λ = ", Es[n].core.κs_λ)
                    #     println("Es[n].core.κs_β = ", Es[n].core.κs_β)

                    #     println("Es[n].κ = ", Es[n].κ)
                    #     println()

                    # end

                    B[m,j], B[m+M,j] = real(out), imag(out)
                end
            end
        end

        # singlets.
        for k = 1:length(A.κ_singlets)
            j += 1
            #println("singlet: j = ", j)
            for m = 1:M

                tmp = NMRSpectraSimulator.evalκsinglets(U[m], A.core.d_singlets,
                A.core.αs_singlets, A.core.Ωs_singlets,
                A.core.β_singlets, A.core.λ0, A.core.κs_λ_singlets)

                # if !isfinite(tmp)
                #     println("Es[n].core.d_singlets = ", Es[n].core.d_singlets)
                #     println("Es[n].core.κs_λ_singlets = ", Es[n].core.κs_λ_singlets)
                #     println("Es[n].core.β_singlets = ", Es[n].core.β_singlets)

                #     println("Es[n].κ_singlets = ", Es[n].κ_singlets)
                # end

                # B[m,j] += real(tmp)
                # B[m+M,j] += imag(tmp)
                B[m,j] = real(tmp)
                B[m+M,j] = imag(tmp)
            end
            #println("B[:,j] = ", B[:,j])
        end

    end

    ### debug.
    # #if !all(isfinite.(B))
    #     JLD.save("debug.jld", "B", B, "j", j)
    #     #println("B is not finite!")
    # #end

    return nothing
end
