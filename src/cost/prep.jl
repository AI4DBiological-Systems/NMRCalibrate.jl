"""
region_min_dist is the minimum horizontal distance between regions, in ppm.
"""
function prepareoptim(cs_config_path,
    molecule_names,
    hz2ppmfunc,
    U_cost0, y_cost0::Vector{Complex{T}},
    As;
    region_min_dist = 0.1) where T <: Real

    cs_delta_group = NMRSpecifyRegions.extractinfofromconfig( cs_config_path, molecule_names)
    Δsys_cs = NMRSpecifyRegions.condenseΔcsconfig(cs_delta_group)

    ΩS0 = NMRSpecifyRegions.getΩS(As)
    ΩS0_ppm = NMRSpecifyRegions.getPs(ΩS0, hz2ppmfunc)

    exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_names, ΩS0_ppm, Δsys_cs;
    min_dist = region_min_dist)

    P_cost0 = hz2ppmfunc.(U_cost0)
    cost_inds, cost_inds_set = NMRSpecifyRegions.getcostinds(exp_info, P_cost0)

    U_cost = U_cost0[cost_inds]
    P_cost = P_cost0[cost_inds]
    y_cost = y_cost0[cost_inds]

    return Δsys_cs, y_cost, U_cost, P_cost, exp_info, cost_inds, cost_inds_set
end



function runalignment(Δ_shifts::Vector{T},
    U_cost,
    y_cost::Vector{Complex{T}},
    LS_inds,
    Es,
    As,
    fs::T,
    SW::T;
    max_iters = 5000,
    xtol_rel = 1e-7,
    ftol_rel = 1e-12,
    maxtime = Inf,
    w = ones(T, length(Es)),
    κ_lb_default = 0.2,
    κ_ub_default = 50.0,
    λ_each_lb = 0.7,
    λ_each_ub = 5.0) where T <: Real

    @assert length(U_cost) == length(y_cost) == length(LS_inds)

    #
    N_d = sum( length(As[n].d) + length(As[n].d_singlets) for n = 1:length(As) )
    N_β = sum( getNβ(As[n]) for n = 1:length(As) )
    N_λ = sum( getNλ(As[n]) for n = 1:length(As) )
    #N_vars = N_d + N_β + N_λ

    #### initial values.
    shift_lb = -copy(Δ_shifts)
    shift_ub = copy(Δ_shifts)
    shift_initial = zeros(T, N_d)

    β_lb = ones(T, N_β) .* (-π)
    β_ub = ones(T, N_β) .* (π)
    β_initial = zeros(T, N_β)

    λ_lb = λ_each_lb .* ones(T, N_λ)
    λ_ub = λ_each_ub .* ones(T, N_λ)
    λ_initial = ones(T, N_λ)

    ### initial phase estimate. Constant phase model.
    #β_est = estimateconstantphase(U_cost, y_cost)
    #βS_base = initializeFIDparameter(αS, β_est)
    #fill!(β_initial, β_est)
    #fill!(β_initial, 0.0)

    ### set up.
    p_lb = [ shift_lb; β_lb; λ_lb ]
    p_ub = [ shift_ub; β_ub; λ_ub ]
    p_initial = [shift_initial; β_initial; λ_initial]
    println("p_lb = ", p_lb)
    println("p_ub = ", p_ub)

    q, updatedfunc, updateβfunc, updateλfunc, updateκfunc, #updatewfunc,
    κ_BLS, getshiftfunc, getβfunc, getλfunc,
    N_vars_set = setupcostcLshiftLS(Es, As, fs, SW,
        LS_inds, U_cost, y_cost, Δ_shifts;
        w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)
    #updatewfunc(p_initial)
    updateκfunc(p_initial)
    parseκ!(Es, κ_BLS)

    # TODO: I am here.
    #q_star = uu->q_star_rad(uu*2*π) # input: Hz.

    ### cost function.
    #U_cost_rad = U_cost .* (2*π)

    obj_func = pp->costcLshift(U_cost, y_cost,
    updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es, κ_BLS, q)

    grad_func = xx->FiniteDiff.finite_difference_gradient(obj_func, xx)

    # force eval for debugging.
    println("Timing: obj_func")
    @time cost_initial = obj_func(p_initial)
    cost_initial = obj_func(p_initial)
    println("obj_func(p_initial) = ", cost_initial)
    println()

    opt = NLopt.Opt(:GN_ESCH, length(p_initial)) # global evolutionary.

    minf, minx, ret, numevals = runNLopt!(  opt,
        p_initial,
        obj_func,
        grad_func,
        p_lb,
        p_ub;
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        maxtime = maxtime)

    println("numevals = ", numevals)
    println("norm(p_initial-minx) = ", norm(p_initial-minx))

    println("Timing: obj_func")
    @time cost_final = obj_func(minx)
    println("obj_func(minx) = ", cost_final)
    println()

    return minx, q, κ_BLS, getshiftfunc, getβfunc, getλfunc,
        obj_func, N_vars_set
end
