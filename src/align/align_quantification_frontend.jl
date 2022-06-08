# based on the contents of align.jl
# consider changing terminology from align to calibrate.
# assumes κs have been updated into As, Bs, Es.

function alignquantificationregion(y_cost::Vector{Complex{T}},
    U_cost,
    P_cost,
    As,
    Bs,
    Es,
    fs,
    SW,
    Δsys_cs,
    a_setp, b_setp, #κs_β_DOFs, κs_β_orderings,
    shift_lb::Vector{T},
    shift_ub::Vector{T};
    LS_inds = 1:length(U_cost),
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 1e-3,
    w_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    # prepare.
    N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
    @assert length(shift_ub) == length(shift_lb) == N_d

    U_rad_cost = U_cost .* (2*π)

    # setup inner optim over β.
    q, updatedfunc, getshiftfunc, N_vars_set,
    run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc,
    q_β = setupcostnesteddwarpw(Es, Bs, As, fs, SW, LS_inds, U_rad_cost,
        y_cost, Δsys_cs, a_setp, b_setp;
        optim_algorithm = β_optim_algorithm,
        w_lb_default = w_lb_default,
        w_ub_default = w_ub_default,
        max_iters = β_max_iters,
        xtol_rel = β_xtol_rel,
        ftol_rel = β_ftol_rel,
        maxtime = β_maxtime)

    # set up outer optim over shifts.
    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )
    p_β = zeros(T, N_β) # persistant buffer.

    obj_func = pp->costnesteddw(U_rad_cost, y_cost, updatedfunc, pp,
    Es, Bs, run_optim, E_BLS, w_BLS, b_BLS, p_β)

    # optim.
    prob = MultistartOptimization.MinimizationProblem(obj_func, shift_lb, shift_ub)

    local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = local_optim_algorithm,
    xtol_rel = xtol_rel,
    maxeval = maxeval,
    maxtime = maxtime)

    multistart_method = MultistartOptimization.TikTak(N_starts)
    ret_mo = MultistartOptimization.multistart_minimization(multistart_method,
        local_method, prob)
    #
    println("ret_mo = ", ret_mo)
    minf = ret_mo.value
    minx = ret_mo.location
    ret = :None
    if hasproperty(ret_mo, :ret)
        ret = ret_mo.ret
    end

    # force w_BLS to update.
    obj_func(minx)

    return obj_func, minf, minx, ret, w_BLS
end

"""
Requires user to supply a surrogate.
"""
function alignquantificationcompound(y::Vector{Complex{T}}, U_y, P_y, As, Bs, Es, fs, SW,
    Δsys_cs, a_setp, b_setp, #κs_β_DOFs, κs_β_orderings::Vector{Vector{Vector{Int}}},
    shift_lb::Vector{T},
    shift_ub::Vector{T},
    cost_inds_set::Vector{Vector{Int}};
    loop_range = 1:length(cost_inds_set),
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 1e-2,
    w_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    N_regions = length(loop_range)
    minxs = Vector{Vector{T}}(undef, N_regions)
    rets = Vector{Symbol}(undef, N_regions)
    minfs = Vector{T}(undef, N_regions)
    obj_funcs = Vector{Function}(undef, N_regions)
    ws = Vector{Vector{T}}(undef, N_regions)

    for r in loop_range
        println("Working on region $(r)")
        obj_funcs[r], minfs[r], minxs[r], rets[r],
        ws[r] = alignquantificationregion(y[cost_inds_set[r]],
            U_y[cost_inds_set[r]],
            P_y[cost_inds_set[r]],
            As,
            Bs,
            Es,
            fs,
            SW,
            Δsys_cs,
            a_setp, b_setp, #κs_β_DOFs, κs_β_orderings,
            shift_lb,
            shift_ub;
            N_starts = N_starts,
            local_optim_algorithm = local_optim_algorithm,
            xtol_rel = xtol_rel,
            maxeval = maxeval,
            maxtime = maxtime,
            β_optim_algorithm = β_optim_algorithm,
            w_lb_default = w_lb_default,
            w_ub_default = w_ub_default,
            β_max_iters = β_max_iters,
            β_xtol_rel = β_xtol_rel,
            β_ftol_rel = β_ftol_rel,
            β_maxtime = β_maxtime)
    end

    return obj_funcs, minfs, minxs, rets, ws
end
