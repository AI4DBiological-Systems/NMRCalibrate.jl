function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

# average Δc vector for each partition element.
function viewaverageΔc(A::NMRSpectraSimulator.CompoundFIDType{T}) where T

    N_spin_groups = length(A.part_inds_compound)
    out = Vector{Vector{Vector{T}}}(undef, N_spin_groups)

    for i = 1:N_spin_groups

        part_size = length(A.part_inds_compound[i])

        out[i] = Vector{Vector{T}}(undef, part_size)


        for (k,inds) in enumerate(A.part_inds_compound[i])

            list_of_Δc_m = A.Δc_m_compound[i][inds]

            # insert error handing here later in case list_of_Δc_m is empty.
            N_spins = length(list_of_Δc_m[1])
            out[i][k] = zeros(T, N_spins)

            N_components = length(list_of_Δc_m)
            for j = 1:N_components
                #
                out[i][k] += list_of_Δc_m[j]
            end
            out[i][k] = out[i][k] ./ N_components

            #out[i][k] = Statistics.mean(tmp)
        end
    end

    return out
end

function findfreqrange(As::Vector{NMRSpectraSimulator.CompoundFIDType{T}}, hz2ppmfunc) where T

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end



function setupcostcLshiftLS(Es,
    Bs,
    fs::T,
    SW::T,
    LS_inds,
    U0,
    y0::Vector{Complex{T}},
    Δ_shifts::Vector{T};
    w = ones(T, length(Es)),
    κ_lb_default = 0.2,
    κ_ub_default = 5.0) where T <: Real

    #N = length(ΩS0)
    #@assert length(ΩS0) == length(αS) == length(w_lb) == length(w_ub)

    ## allocate buffers.
    #U_rad_LS = U0[LS_inds] .* (2*π)
    #A, w_BLS = setupLSwcL(U_rad_LS, αS)
    #
    #tmp = S_U0[LS_inds]
    #LS_b = [real.(tmp); imag(tmp)]

    # model.
    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)

    ##### update functions.
    N_d = sum( length(Bs[n].d) + length(Bs[n].d_singlets) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )
    N_λ = sum( getNλ(Bs[n]) for n = 1:length(Bs) )

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixtured!(Bs, pp, st_ind_d, fs, SW, Δ_shifts)

    st_ind_β = fin_ind_d + 1
    fin_ind_β = st_ind_β + N_β - 1
    updateβfunc = pp->updateβ!(Bs, pp, st_ind_β)
    #N_β = sum( sum(length(Bs[n].κs_β[l]) for l = 1:length(Bs[n].κs_β)) + length(Bs[n].β_singlets) for n = 1:length(Bs) )

    #λupdate.
    st_ind_λ = fin_ind_β + 1
    fin_ind_λ = st_ind_λ + N_λ -1
    updateλfunc = pp->updateλ!(Bs, pp, st_ind_λ)

    N_vars_set = [N_d; N_β; N_λ]

    ### LS κ.
    U_LS = U0[LS_inds]

    N_κ, N_κ_singlets = countκ(Es)
    N_κ_vars = N_κ + N_κ_singlets
    E_BLS, κ_BLS = setupupdatew(length(U_LS), N_κ_vars)

    κ_lb = ones(N_κ_vars) .* κ_lb_default
    κ_ub = ones(N_κ_vars) .* κ_ub_default


    b_BLS = [real.(y0[LS_inds]); imag.(y0[LS_inds])]

    updateκfunc = xx->updateκ!(E_BLS, b_BLS, κ_BLS,
    U_LS, Es, w, κ_lb, κ_ub)

    #### extract parameters from p.
    getshiftfunc = pp->pp[st_ind_d:fin_ind_d]
    getβfunc = pp->pp[st_ind_β:fin_ind_β]
    getλfunc = pp->pp[st_ind_λ:fin_ind_λ]

    return f, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
    κ_BLS, getshiftfunc, getβfunc, getλfunc, N_vars_set
end

function mergeinds(inds_set::Vector{Vector{Int}}, merge_set) where T

    # check.
    for l = 1:length(merge_set)
        @assert 0 < l <= length(inds_set)
    end

    # merge.
    inds = Vector{Int}(undef, 0)
    for i in merge_set
        push!(inds, inds_set[i]...)
    end

    unique!(inds)
    return inds
end


function costcLshift(U,
    S_U::Vector{Complex{T}},
    updatedfunc, updateβfunc, updateλfunc, updateκfunc,
    p::Vector{T}, Es, κ_BLS,
    f)::T where T <: Real

    #println("p = ", p)
    ### update main optim variables.
    updatedfunc(p)
    updateβfunc(p)
    updateλfunc(p)

    updateκfunc(p)
    parseκ!(Es, κ_BLS)
    # try
    #     updateκfunc(p)
    # catch err_var
    #     println("error: ", err_var)
    #     println("p = ", p)
    # end


    ## l-2 costfunc.
    cost = zero(T)
    for m = 1:length(S_U)

        f_u = f(U[m])

        cost += abs2( f_u - S_U[m] )
        #cost += (abs2(f_u) - abs2(S_U[m]))^2

        #cost += abs2( f_buffer[m] - DFT_s[m] ) + (abs2(f_buffer[m]) - abs2(DFT_s[m]))^2
    end

    return cost
end



function runNLopt!(  opt,
    p0::Vector{T},
    obj_func,
    grad_func,
    p_lbs,
    p_ubs;
    max_iters = 10000,
    xtol_rel = 1e-12,
    ftol_rel = 1e-12,
    maxtime = Inf) where T

    @assert length(p0) == length(p_lbs) == length(p_ubs)

    opt.maxeval = max_iters
    opt.lower_bounds = p_lbs
    opt.upper_bounds = p_ubs
    opt.xtol_rel = xtol_rel
    opt.ftol_rel = ftol_rel
    opt.maxtime = maxtime


    opt.min_objective = (xx, gg)->genericNLoptcostfunc!(xx, gg, obj_func, grad_func)

    # optimize.
    (minf, minx, ret) = NLopt.optimize(opt, p0)

    N_evals = opt.numevals

    return minf, minx, ret, N_evals
end

function genericNLoptcostfunc!(x::Vector{T}, df_x, f, df)::T where T <: Real

    #
    if length(df_x) > 0
        df_x[:] = df(x)
    end

    return f(x)
end
