
"""
an example of optim_algorithm is NLopt.LN_BOBYQA
"""
function setupβLSsolverMultistartoptim(optim_algorithm,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
    As::Vector{NMRSpectraSimulator.CompoundFIDType{T, SST}},
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    κs_β_DOFs, κs_β_orderings;
    κ_lb_default = 0.2,
    κ_ub_default = 5.0,
    max_iters = 50,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf,
    N_starts = 100,
    inner_xtol_rel = 1e-12,
    inner_maxeval = 100,
    inner_maxtime = Inf) where {T,SST}

    #
    N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    β_lb = ones(T, N_β) .* (-π)
    β_ub = ones(T, N_β) .* (π)
    #β_initial = zeros(T, N_β)

    p_lb = β_lb
    p_ub = β_ub

    q, updateβfunc, updateκfunc, E_BLS, κ_BLS, b_BLS,
    getβfunc = setupcostβLS(Es, Bs, As, LS_inds, U_rad_cost, y_cost,
        κs_β_DOFs,
        κs_β_orderings;
        κ_lb_default = κ_lb_default,
        κ_ub_default = κ_ub_default)

    f = pp->costβLS(U_rad_cost, y_cost, updateβfunc, updateκfunc, pp, Es,
    E_BLS, κ_BLS, b_BLS, q)

    P = MultistartOptimization.MinimizationProblem(f,
        p_lb, p_ub)
    #
    local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = optim_algorithm,
    xtol_rel = inner_xtol_rel,
    maxeval = inner_maxeval,
    maxtime = inner_maxtime)

    multistart_method = MultistartOptimization.TikTak(N_starts)
    run_optim = pp->MultistartOptimization.multistart_minimization(multistart_method,
        local_method, P)

    return run_optim, f, E_BLS, κ_BLS, b_BLS, updateβfunc, q
end

function setupβLSsolver(optim_algorithm,
    Es,
    Bs,
    As,
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    κs_β_DOFs::Vector{Vector{Int}},
    κs_β_orderings::Vector{Vector{Vector{Int}}};
    κ_lb_default = 0.2,
    κ_ub_default = 5.0,
    max_iters = 50,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf) where T

    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )

    β_lb = ones(T, N_β) .* (-π)
    β_ub = ones(T, N_β) .* (π)
    #β_initial = zeros(T, N_β)

    p_lb = β_lb
    p_ub = β_ub

    q, updateβfunc, updateκfunc, E_BLS, κ_BLS, b_BLS,
    getβfunc = setupcostβLS(Es, Bs, As, LS_inds, U_rad_cost, y_cost,
    κs_β_DOFs, κs_β_orderings;
        κ_lb_default = κ_lb_default,
        κ_ub_default = κ_ub_default)

    # q is simulated spectra. f is cost function.
    f = pp->costβLS(U_rad_cost, y_cost, updateβfunc, updateκfunc, pp, Es,
    E_BLS, κ_BLS, b_BLS, q)

    df = xx->FiniteDiff.finite_difference_gradient(f, xx)

    opt = NLopt.Opt(optim_algorithm, N_β)

    run_optim = pp->runNLopt!(opt,
        pp,
        f,
        df,
        p_lb,
        p_ub;
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        maxtime = maxtime)

    return run_optim, f, E_BLS, κ_BLS, b_BLS, updateβfunc, updateκfunc, q
end

"""
κ version.
"""
function setupcostβLS(Es::Vector{NMRSpectraSimulator.καFIDModelType{T, SST}},
    Bs,
    As,
    LS_inds,
    U0_rad,
    y0::Vector{Complex{T}},
    κs_β_DOFs,
    κs_β_orderings;
    κ_lb_default = 0.2,
    κ_ub_default = 5.0) where {T <: Real, SST}

    w = ones(T, length(Es))

    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )


    st_ind_β = 1
    fin_ind_β = st_ind_β + N_β - 1
    #updateβfunc = pp->updateβ!(Bs, κs_β_orderings, κs_β_DOFs, pp, st_ind_β)
    updateβfunc = pp->updateβ!(Bs, pp, st_ind_β)


    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Bs, Es; w = w)

    ### LS κ.
    U_rad_LS = U0_rad[LS_inds]

    N_κ, N_κ_singlets = countκ(Es)
    N_κ_vars = N_κ + N_κ_singlets
    E_BLS, κ_BLS, b_BLS = setupupdateLS(length(U_rad_LS), N_κ_vars, y0[LS_inds])

    κ_lb = ones(N_κ_vars) .* κ_lb_default
    κ_ub = ones(N_κ_vars) .* κ_ub_default

    updateκfunc = xx->updateκ2!(E_BLS, b_BLS, κ_BLS,
    U_rad_LS, Es, As, w, κ_lb, κ_ub)

    #### extract parameters from p.
    getβfunc = pp->pp[st_ind_β:fin_ind_β]

    return f, updateβfunc, updateκfunc, E_BLS, κ_BLS, b_BLS, getβfunc
end

function costβLS(U_rad,
    S_U::Vector{Complex{T}},
    updateβfunc, updateκfunc,
    p::Vector{T}, Es,
    E_BLS::Matrix{Complex{T}}, κ_BLS::Vector{T}, b_BLS,
    f)::T where T <: Real

    updateβfunc(p)

    # try updateκfunc(p)
    # catch e
    #     println("error in updateκfunc(p)")
    # 	bt = catch_backtrace()
    # 	showerror(stdout, e, bt)
    # 	rethrow(e)
    # end

    updateκfunc(p)
    parseκ!(Es, κ_BLS)
    #println("all(isfinite.(Es)) = ", all(isfinite.(Es)))
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
    # tmp = B*κ_BLS - b_BLS
    # cost = dot(tmp, tmp)
    cost = norm(reinterpret(T, E_BLS)*κ_BLS - b_BLS)^2 # compact version.

    return cost
end
