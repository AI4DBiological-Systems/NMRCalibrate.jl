

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
        A = As[n]

        for m = 1:M

            s = zero(Complex{T})
            for i in eachindex(E.??s_??)
                for k in eachindex(E.??s_??[i])
                    s += E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.??s_??[i])
                end
            end
            B[m,n] = s

            #B[m,n] = sum( sum( E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.??s_??[i]) for k in eachindex(E.??s_??[i])) for i in eachindex(E.??s_??))

            # singlets.
            for k = 1:length(E.??s_??_singlets)

                B[m,n] += NMRSpectraSimulator.evalsinglets(U_rad[m], E.core.d_singlets,
                A.??s_singlets, A.??s_singlets,
                E.core.??_singlets, E.core.??0, E.core.??s_??_singlets, E.??s_??_singlets)
            end

        end
    end

    # if !all(isfinite.(B))
    #     println("problem")
    # end

    return nothing
end


######### ??, compensation for ??.


function update??2!(  A::Matrix{Complex{T}},
    #r_buffer::Vector{T},
    b,
    ??::Vector{T},
    U_rad_LS,
    Es,
    As,
    w::Vector{T},
    ??_lower::Vector{T},
    ??_upper::Vector{T}) where T <: Real

    #@assert size(A) == (length(b), length(??S))
    @assert length(b) == 2*length(U_rad_LS)

    evaldesignmatrix??!(A, U_rad_LS, Es, As, w)

    ### LS solve.
    # w[:] = NonNegLeastSquares.nonneg_lsq(A, b; alg = :fnnls)
    # status_flag = true

    ??[:], status_flag = solveBLScL(A, b, ??_lower, ??_upper)
    # ??[:] = A\b
    # clamp!(??, ??_lower[1], ??_upper[1])
    # status_flag = true

    return status_flag
end

function evaldesignmatrix??!(B::Matrix{Complex{T}},
    U_rad,
    Es,
    As,
    w::Vector{T}) where T <: Real

    #
    M = length(U_rad)
    N = length(Es)

    N_??, N_??_singlets = count??(Es)
    @assert size(B) == (M, N_?? + N_??_singlets)
    fill!(B, Inf) # debug.

    reset??!(Es)
    j = 0

    # loop over each ?? partition element in Es.
    for n = 1:N
        E = Es[n]
        A = As[n]
        @assert length(E.??s_??) == length(E.core.qs)

        #spin system.
        #for i = 1:length(E.??s_??)
        for i in eachindex(E.??s_??)

            # partition
            #for k = 1:length(E.??s_??[i])
            for k in eachindex(E.??s_??[i])

                j += 1

                # # loop over each fit position.
                # for m = 1:M
                #
                #     # taken from evalitproxysys()
                #     r = U_rad[m] - E.core.ss_params.d[i]
                #
                #     #B[m,j] = w[n]*E.core.qs[i][k](r, E.core.ss_params.??s_??[i])
                #     #w[n]*E.core.qs[i][k](r, E.core.ss_params.??s_??[i])
                #     #B[m,j] = 1.9
                #
                #     out = w[n]*E.core.qs[i][k](r, E.core.ss_params.??s_??[i])
                #
                #     #Q[j,m] = w[n]*E.core.qs[i][k](r, E.core.ss_params.??s_??[i])
                # end

                #out = collect( w[n]*E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.??s_??[i]) for m = 1:M )
                B[:,j] = collect( w[n]*E.core.qs[i][k](U_rad[m] - E.core.ss_params.d[i], E.core.ss_params.??s_??[i]) for m = 1:M )
            end
        end

        # singlets.
        for k = 1:length(E.??s_??_singlets)
            j += 1

            # for m = 1:M
            #
            #     B[m,j] = w[n]*NMRSpectraSimulator.evalsinglets(U_rad[m], E.core.d_singlets,
            #     A.??s_singlets, A.??s_singlets,
            #     E.core.??_singlets, E.core.??0, E.core.??s_??_singlets, E.??s_??_singlets)
            # end

            B[:,j] = collect( w[n]*NMRSpectraSimulator.evalsinglets(U_rad[m], E.core.d_singlets,
            A.??s_singlets, A.??s_singlets,
            E.core.??_singlets, E.core.??0, E.core.??s_??_singlets, E.??s_??_singlets) for m = 1:M )


        end
    end

    #println("all(isfinite.(B)) = ", all(isfinite.(B)))
    @assert all(isfinite.(B))

    return nothing
end
