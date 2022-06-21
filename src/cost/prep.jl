"""
region_min_dist is the minimum horizontal distance between regions, in ppm.
"""
function prepareoptim(config_path::String,
    molecule_names,
    hz2ppmfunc,
    U_cost0, y_cost0::Vector{Complex{T}},
    As;
    region_min_dist = 0.1) where T <: Real

    ## parse from config file.
    config_dict = Dict()
    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        config_dict = JSON.parsefile(config_path)
    end

    N_compounds = length(molecule_names)

    λ_lbs = Vector{Vector{T}}(undef, N_compounds)
    λ_ubs = Vector{Vector{T}}(undef, N_compounds)
    Δsys_cs = Vector{Vector{T}}(undef, N_compounds)
    κs_β_orderings = Vector{Vector{Vector{Int}}}(undef, N_compounds)
    κs_β_DOFs = Vector{Vector{Int}}(undef, N_compounds)

    for n = 1:N_compounds

        dict = config_dict[molecule_names[n]] # TODO graceful error-handle.

        λ_lbs[n] = convert(Vector{T}, dict["λ_lb"])
        λ_ubs[n] = convert(Vector{T}, dict["λ_ub"])
        Δsys_cs[n] = convert(Vector{T}, dict["maximum chemical shift"])
        κs_β_orderings[n] = convert(Vector{Vector{Int}}, dict["κs_β ordering"])
        κs_β_DOFs[n] = convert(Vector{Int}, dict["κs_β degrees of freedom"])

    end


    # cs_delta_group = NMRSpecifyRegions.extractinfofromconfig( cs_config_path, molecule_names)
    # Δsys_cs = NMRSpecifyRegions.condenseΔcsconfig(cs_delta_group)

    ## get regions.
    ΩS0 = NMRSpectraSimulator.getΩS(As)
    ΩS0_ppm = NMRSpectraSimulator.getPs(ΩS0, hz2ppmfunc)

    exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_names, ΩS0_ppm, Δsys_cs;
    min_dist = region_min_dist)

    P_cost0 = hz2ppmfunc.(U_cost0)
    cost_inds, cost_inds_set = NMRSpecifyRegions.getcostinds(exp_info, P_cost0)

    U_cost = U_cost0[cost_inds]
    P_cost = P_cost0[cost_inds]
    y_cost = y_cost0[cost_inds]

    return Δsys_cs, y_cost, U_cost, P_cost, exp_info, cost_inds, cost_inds_set,
    λ_lbs,
    λ_ubs,
    κs_β_orderings,
    κs_β_DOFs
end


function fitβLSκ(U_rad_cost,
    y_cost::Vector{Complex{T}},
    LS_inds,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
    As::Vector{NMRSpectraSimulator.CompoundFIDType{T, SST}},
    β_initial::Vector{T},
    β_lb,
    β_ub;
    optim_algorithm::Symbol = :GN_ESCH,
    max_iters = 5000,
    xtol_rel = 1e-7,
    ftol_rel = 1e-12,
    maxtime = Inf,
    w = ones(T, length(Es)),
    κ_lb_default = 0.2,
    κ_ub_default = 50.0) where {T <: Real, SST}

    @assert length(U_rad_cost) == length(y_cost) == length(LS_inds)

    q, updateβfunc, updateκfunc,
    E_BLS, κ_BLS, b_BLS, getβfunc = setupcostβLS(Es, As, LS_inds, U_rad_cost, y_cost)

    # updateκfunc(0.0)
    # NMRCalibrate.parseκ!(Es, κ_BLS)
    obj_func = pp->costβLS(U_rad_cost, y_cost, updateβfunc, updateκfunc, pp, Es, E_BLS, κ_BLS, b_BLS, q)
    grad_func = xx->FiniteDiff.finite_difference_gradient(obj_func, xx)

    # force eval for debugging.
    println("Timing: obj_func")
    @time cost_initial = obj_func(β_initial)
    cost_initial = obj_func(β_initial)
    println("obj_func(β_initial) = ", cost_initial)
    println()

    opt = NLopt.Opt(optim_algorithm, length(β_initial)) # global evolutionary.
println(β_initial)
    minf, minx, ret, numevals = runNLopt!(  opt,
        β_initial,
        obj_func,
        grad_func,
        β_lb,
        β_ub;
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        maxtime = maxtime)

    println("Timing: obj_func")
    @time cost_final = obj_func(minx)
    println("obj_func(minx) = ", cost_final)
    println()

    return minx, q, κ_BLS, getβfunc, obj_func
end
