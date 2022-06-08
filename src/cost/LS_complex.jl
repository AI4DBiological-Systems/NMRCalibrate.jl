

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

function setupupdateLS(N_positions::Int, N_compounds::Int, y::Vector{Complex{T}}) where T <: Real

    A = Matrix{Complex{T}}(undef, N_positions, N_compounds)
    w = zeros(T, N_compounds)
    b = collect(reinterpret(T, y))

    return A, w, b
end

"""
A is an M x L complex-valued matrix, b is 2*M 1D array.
Updates solution to w (mutates it), and mutates intermediate variables A, b.
"""
function updatew!(  A::Matrix{Complex{T}},
    b,
    w::Vector{T},
    U_rad_LS,
    Es,
    As,
    w_lower::Vector{T},
    w_upper::Vector{T}) where T <: Real

    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrixw2!(A, U_rad_LS, Es, As)

    # ## check.
    # A2 = copy(A)
    # evaldesignmatrixw!(A, U_rad_LS, Es, As)
    #
    # if norm(A-A2) > 1e-10
    #     println("problem)")
    # end

    ### LS solve.
    w[:], status_flag = solveBLScL(A, b, w_lower, w_upper)

    return status_flag
end


function evaldesignmatrixw!(B::Matrix{Complex{T}},
    U_rad, Es, As) where T <: Real

    M = length(U_rad)
    N = length(Es)

    @assert size(B) == (M,N)

    fill!(B, NaN) # debug.

    for n = 1:N
        for m = 1:M
            # speed up later.
            B[m,n] = NMRSpectraSimulator.evalitpproxycompound(U_rad[m], As[n], Es[n])
        end
    end

    return nothing
end

function evaldesignmatrixw2!(B::Matrix{Complex{T}},
    U_rad, Es, As) where T <: Real

    M = length(U_rad)
    N = length(Es)

    @assert size(B) == (M,N)

    fill!(B, NaN) # debug.

    for n = 1:N
        E = Es[n]

        for m = 1:M

            s = zero(Complex{T})
            for i in eachindex(E.κs_α)
                for k in eachindex(E.κs_α[i])
                    s += E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.κs_λ[i])
                end
            end
            B[m,n] = s

            #B[m,n] = sum( sum( E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.κs_λ[i]) for k in eachindex(E.κs_α[i])) for i in eachindex(E.κs_α))

            # singlets.
            for k = 1:length(E.κs_α_singlets)

                B[m,n] += NMRSpectraSimulator.evalsinglets(U_rad[m], E.core.d_singlets,
                A.αs_singlets, A.Ωs_singlets,
                E.core.β_singlets, E.core.λ0, E.core.κs_λ_singlets, E.κs_α_singlets)
            end

        end
    end

    # if !all(isfinite.(B))
    #     println("problem")
    # end

    return nothing
end


######### κ, compensation for α.


function updateκ2!(  A::Matrix{Complex{T}},
    #r_buffer::Vector{T},
    b,
    κ::Vector{T},
    U_rad_LS,
    Es,
    As,
    w::Vector{T},
    κ_lower::Vector{T},
    κ_upper::Vector{T}) where T <: Real

    #@assert size(A) == (length(b), length(αS))
    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrixκ!(A, U_rad_LS, Es, As, w)

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
    Es,
    As,
    w::Vector{T}) where T <: Real

    #
    M = length(U_rad)
    N = length(Es)

    N_κ, N_κ_singlets = countκ(Es)
    @assert size(B) == (M, N_κ + N_κ_singlets)
    fill!(B, Inf) # debug.

    resetκ!(Es)
    j = 0

    # loop over each κ partition element in Es.
    for n = 1:N
        E = Es[n]
        A = As[n]
        @assert length(E.κs_α) == length(E.core.qs)

        #spin system.
        #for i = 1:length(E.κs_α)
        for i in eachindex(E.κs_α)

            # partition
            #for k = 1:length(E.κs_α[i])
            for k in eachindex(E.κs_α[i])

                j += 1

                # # loop over each fit position.
                # for m = 1:M
                #
                #     # taken from evalitproxysys()
                #     r = U_rad[m] - E.core.ss_params.d[i]
                #
                #     #B[m,j] = w[n]*E.core.qs[i][k](r, E.core.ss_params.κs_λ[i])
                #     #w[n]*E.core.qs[i][k](r, E.core.ss_params.κs_λ[i])
                #     #B[m,j] = 1.9
                #
                #     out = w[n]*E.core.qs[i][k](r, E.core.ss_params.κs_λ[i])
                #
                #     #Q[j,m] = w[n]*E.core.qs[i][k](r, E.core.ss_params.κs_λ[i])
                # end

                #out = collect( w[n]*E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.κs_λ[i]) for m = 1:M )
                B[:,j] = collect( w[n]*E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.κs_λ[i]) for m = 1:M )
            end
        end

        # singlets.
        for k = 1:length(E.κs_α_singlets)
            j += 1

            # for m = 1:M
            #
            #     B[m,j] = w[n]*NMRSpectraSimulator.evalsinglets(U_rad[m], E.core.d_singlets,
            #     A.αs_singlets, A.Ωs_singlets,
            #     E.core.β_singlets, E.core.λ0, E.core.κs_λ_singlets, E.κs_α_singlets)
            # end

            B[:,j] = collect( w[n]*NMRSpectraSimulator.evalsinglets(U_rad[m], E.core.d_singlets,
            A.αs_singlets, A.Ωs_singlets,
            E.core.β_singlets, E.core.λ0, E.core.κs_λ_singlets, E.κs_α_singlets) for m = 1:M )


        end
    end

    #println("all(isfinite.(B)) = ", all(isfinite.(B)))
    @assert all(isfinite.(B))

    return nothing
end
