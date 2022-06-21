
function setupβLSsolverw(optim_algorithm,
    Es,
    Bs,
    As,
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    κs_β_DOFs, κs_β_orderings;
    w_lb_default = 0.2,
    w_ub_default = 5.0,
    β_max_iters = 50,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T

    N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )

    β_lb = ones(T, N_β) .* (-π)
    β_ub = ones(T, N_β) .* (π)

    p_lb = β_lb
    p_ub = β_ub

    q, updateβfunc, updatewfunc, E_BLS, w_BLS, b_BLS,
    getβfunc = setupcostβLSw(Es, Bs, As, LS_inds, U_rad_cost, y_cost,
        κs_β_DOFs, κs_β_orderings;
        w_lb_default = w_lb_default,
        w_ub_default = w_ub_default)

    # q is simulated spectra. f is cost function.
    f = pp->costβLSw(U_rad_cost, y_cost, updateβfunc, updatewfunc, pp, Es,
    E_BLS, w_BLS, b_BLS, q)

    df = xx->FiniteDiff.finite_difference_gradient(f, xx)

    opt = NLopt.Opt(optim_algorithm, N_β)

    run_optim = pp->runNLopt!(opt,
        pp,
        f,
        df,
        p_lb,
        p_ub;
        max_iters = β_max_iters,
        xtol_rel = β_xtol_rel,
        ftol_rel = β_ftol_rel,
        maxtime = β_maxtime)

    return run_optim, f, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc, q
end

"""
w version.
"""
function setupcostβLSw(Es::Vector{NMRSpectraSimulator.καFIDModelType{T, SST}},
    Bs,
    As,
    LS_inds,
    U0_rad,
    y0::Vector{Complex{T}},
    κs_β_DOFs, κs_β_orderings;
    w_lb_default = 0.2,
    w_ub_default = 5.0) where {T <: Real, SST}

    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )


    st_ind_β = 1
    fin_ind_β = st_ind_β + N_β - 1
    #updateβfunc = pp->updateβ!(Bs, κs_β_orderings, κs_β_DOFs, pp, st_ind_β)
    updateβfunc = pp->updateβ!(Bs, pp, st_ind_β)

    ### LS w.
    U_rad_LS = U0_rad[LS_inds]

    N_w_vars = length(Es)
    E_BLS, w_BLS, b_BLS = setupupdateLS(length(U_rad_LS), N_w_vars, y0[LS_inds])

    w_lb = ones(N_w_vars) .* w_lb_default
    w_ub = ones(N_w_vars) .* w_ub_default

    updatewfunc = xx->updatew!(E_BLS, b_BLS, w_BLS,
    U_rad_LS, Es, As, w_lb, w_ub)

    #### extract parameters from p.
    getβfunc = pp->pp[st_ind_β:fin_ind_β]

    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w_BLS)

    return f, updateβfunc, updatewfunc, E_BLS, w_BLS, b_BLS, getβfunc
end

function costβLSw(U_rad,
    S_U::Vector{Complex{T}},
    updateβfunc, updatewfunc,
    p::Vector{T}, Es,
    E_BLS::Matrix{Complex{T}}, w_BLS::Vector{T}, b_BLS,
    f)::T where T <: Real

    updateβfunc(p)
    updatewfunc(1.0)

    # ## l-2 costfunc.
    # cost = zero(T)
    # for m = 1:length(S_U)
    #
    #     f_u = f(U_rad[m])
    #
    #     cost += abs2( f_u - S_U[m] )
    # end
    #
    # # faster l-2 cost compute.
    # B = reinterpret(T, E_BLS)
    # tmp = B*w_BLS - b_BLS
    # cost = dot(tmp, tmp)
    cost = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2 # compact version.

    return cost
end
