# based on the contents of align.jl
# consider changing terminology from align to calibrate.

function alignregion(y_cost::Vector{Complex{T}},
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
    w = ones(T, length(As)),
    LS_inds = 1:length(U_cost),
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = 1e-2,
    κ_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    # prepare.
    N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )
    @assert length(shift_ub) == length(shift_lb) == N_d

    U_rad_cost = U_cost .* (2*π)

    # setup inner optim over β.
    q, updatedfunc, getshiftfunc, N_vars_set,
    run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc,
    q_β = setupcostnesteddwarp(Es, Bs, As, fs, SW, LS_inds, U_rad_cost,
        y_cost, Δsys_cs, a_setp, b_setp; #κs_β_DOFs, κs_β_orderings;
        w = w,
        optim_algorithm = β_optim_algorithm,
        κ_lb_default = κ_lb_default,
        κ_ub_default = κ_ub_default,
        max_iters = β_max_iters,
        xtol_rel = β_xtol_rel,
        ftol_rel = β_ftol_rel,
        maxtime = β_maxtime)

    # set up outer optim over shifts.
    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )
    p_β = zeros(T, N_β) # persistant buffer.

    obj_func = pp->costnestedd(U_rad_cost, y_cost, updatedfunc, pp,
    #Es, Bs, κs_β_orderings, κs_β_DOFs, q, run_optim, E_BLS, κ_BLS, b_BLS, p_β)
    Es, Bs, run_optim, E_BLS, κ_BLS, b_BLS, p_β)

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

    return obj_func, minf, minx, ret
end

"""
Requires user to supply a surrogate.
"""
function aligncompound(y::Vector{Complex{T}}, U_y, P_y, As, Bs, Es, fs, SW,
    Δsys_cs, a_setp, b_setp, #κs_β_DOFs, κs_β_orderings::Vector{Vector{Vector{Int}}},
    shift_lb::Vector{T},
    shift_ub::Vector{T},
    cost_inds_set::Vector{Vector{Int}};
    loop_range = 1:length(cost_inds_set),
    w = ones(T, length(As)),
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = 1e-2,
    κ_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    N_regions = length(loop_range)
    minxs = Vector{Vector{T}}(undef, N_regions)
    rets = Vector{Symbol}(undef, N_regions)
    minfs = Vector{T}(undef, N_regions)
    obj_funcs = Vector{Function}(undef, N_regions)

    for r in loop_range
        println("Working on region $(r)")
        obj_funcs[r], minfs[r], minxs[r], rets[r] = alignregion(y[cost_inds_set[r]],
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
            w = w,
            N_starts = N_starts,
            local_optim_algorithm = local_optim_algorithm,
            xtol_rel = xtol_rel,
            maxeval = maxeval,
            maxtime = maxtime,
            β_optim_algorithm = β_optim_algorithm,
            κ_lb_default = κ_lb_default,
            κ_ub_default = κ_ub_default,
            β_max_iters = β_max_iters,
            β_xtol_rel = β_xtol_rel,
            β_ftol_rel = β_ftol_rel,
            β_maxtime = β_maxtime)
    end

    return obj_funcs, minfs, minxs, rets
end


function aligncompoundsingleregion(y::Vector{Complex{T}}, U_y, P_y, As, Bs, Es, fs, SW,
    Δsys_cs, a_setp, b_setp, #κs_β_DOFs, κs_β_orderings::Vector{Vector{Vector{Int}}},
    shift_lb::Vector{T},
    shift_ub::Vector{T},
    cost_inds_set::Vector{Vector{Int}};
    loop_range = 1:length(cost_inds_set),
    w = ones(T, length(As)),
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = 1e-2,
    κ_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where T <: Real

    tmp = collect( cost_inds_set[r] for r = 1:length(cost_inds_set) )
    inds_redundant = NMRSpectraSimulator.combinevectors(tmp)
    inds = inds_redundant[indexin( unique(inds_redundant), inds_redundant)]

    U_cost = U_y[inds]
    P_cost = P_y[inds]
    y_cost = y[inds]

    obj_func, minf, minx, ret = alignregion(y_cost,
        U_cost,
        P_cost,
        As,
        Bs,
        Es,
        fs,
        SW,
        Δsys_cs,
        a_setp, b_setp, #κs_β_DOFs, κs_β_orderings,
        shift_lb,
        shift_ub;
        w = w,
        N_starts = 100,
        local_optim_algorithm = NLopt.LN_BOBYQA,
        xtol_rel = 1e-3,
        maxeval = 50,
        maxtime = Inf,
        β_optim_algorithm = :GN_DIRECT_L,
        κ_lb_default = 1e-2,
        κ_ub_default = 1e2,
        β_max_iters = 500,
        β_xtol_rel = 1e-9,
        β_ftol_rel = 1e-9,
        β_maxtime = Inf)

    return obj_func, minf, minx, ret, P_cost, y_cost
end


"""
Includes surrogate construction. Taken from prep_save.jl
"""
function alignproject(save_path::String,
    SH_config_path,
    surrogate_config_path,
    fit_config_path,
    H_params_path,
    dict_compound_to_filename,
    molecule_names::Vector{String},
    w::Vector{T},
    dummy_SSFID::SST,
    y::Vector{Complex{T}},
    U_y,
    SW,
    fs,
    ν_0ppm,
    λ0::T;
    u_offset = 0.5,
    unique_cs_atol = 1e-6,
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05,
    Δc_partition_radius = 0.3,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05,
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.7,
    κ_λ_ub_default = 1.5,
    offset_ppm::T = 0.3,
    region_min_dist = 0.1,
    sigmoid_lb = 0.1,
    sigmoid_ub = 0.7,
    N_sigmoid_samples = 10,
    sigmoid_optim_algorithm = :LN_BOBYOA,
    save_flag = false,
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf) where {T, SST}

    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    ### surrogate.
    # get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
    MEs = NMRSpectraSimulator.getmageqinfomixture(molecule_names,
        H_params_path,
        dict_compound_to_filename;
        unique_cs_atol = unique_cs_atol)

    mixture_params = NMRSpectraSimulator.setupmixtureSH(molecule_names,
        H_params_path, dict_compound_to_filename, fs, SW,
        ν_0ppm;
        config_path = SH_config_path,
        tol_coherence = tol_coherence,
        α_relative_threshold = α_relative_threshold,
        Δc_partition_radius = Δc_partition_radius)
    As = mixture_params

    ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
    ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

    u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
    u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)


    Bs = NMRSpectraSimulator.fitproxies(As, dummy_SSFID, λ0;
        names = molecule_names,
        config_path = surrogate_config_path,
        Δcs_max_scalar_default = Δcs_max_scalar_default,
        κ_λ_lb_default = κ_λ_lb_default,
        κ_λ_ub_default = κ_λ_ub_default,
        u_min = u_min,
        u_max = u_max,
        Δr_default = Δr_default,
        Δκ_λ_default = Δκ_λ_default)

    ### prepare positions.
    P_y = hz2ppmfunc.(U_y)

    ΩS0 = NMRSpectraSimulator.getΩS(As)
    ΩS0_ppm = NMRSpectraSimulator.getPs(ΩS0, hz2ppmfunc)

    Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds, cost_inds_set,
        λ_lbs, λ_ubs, κs_β_orderings,
        κs_β_DOFs = prepareoptim(fit_config_path, molecule_names,
        hz2ppmfunc, U_y, y, As;
        region_min_dist = region_min_dist)

    ### align.
    N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )
    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )

    shift_lb = -ones(T, N_d)
    shift_ub = ones(T, N_d)
    w = ones(T, length(As))

    Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )

    a_setp, b_setp, minxs,
        rets = setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
    #

    obj_funcs, minfs, minxs, rets = aligncompound(y,
        U_y,
        P_y,
        As,
        Bs,
        Es,
        fs,
        SW,
        Δsys_cs,
        a_setp, b_setp, #κs_β_DOFs, κs_β_orderings,
        shift_lb,
        shift_ub,
        cost_inds_set;
        w = w,
        N_starts = N_starts,
        local_optim_algorithm = local_optim_algorithm,
        xtol_rel = xtol_rel,
        maxeval = maxeval,
        maxtime = maxtime,
        β_optim_algorithm = :GN_DIRECT_L,
        κ_lb_default = κ_λ_lb_default,
        κ_ub_default = κ_λ_ub_default,
        β_max_iters = β_max_iters,
        β_xtol_rel = β_xtol_rel,
        β_ftol_rel = β_ftol_rel,
        β_maxtime = β_maxtime)

    if save_flag
        BSON.bson(save_path, region_min_dist = region_min_dist,
        minfs = minfs,
        minxs = minxs,
        rets = rets,
        cost_inds_set = cost_inds_set,
        Es = Es,
        As = As,
        w = w,
        molecule_names = molecule_names)
    end

    return obj_funcs, minfs, minxs, rets, As, Bs, Es, cost_inds_set, u_min, u_max
end
