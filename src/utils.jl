

function extractinfofromconfig( config_path::String,
    molecule_names::Vector{String}) where T <: Real

    N_compounds = length(molecule_names)

    # read.
    file_strings = readlines(config_path)

    cs_delta_group = Vector{Vector{Vector{Float64}}}(undef, N_compounds)
    #λ_group_labels = Vector{Vector{Vector{Int}}}(undef, N_compounds)

    for n = 1:N_compounds

        cs_delta_group[n] = extractmoleculeinfofromconfig(file_strings,
                                    molecule_names[n], config_path)
    end

    return cs_delta_group#, λ_group_labels
end

function condenseΔcsconfig(cs_delta_group::Vector{Vector{Vector{T}}}) where T
    #
    out = Vector{Vector{T}}(undef, length(cs_delta_group))

    for n = 1:length(cs_delta_group)
        out[n] = Vector{T}(undef, length(cs_delta_group[n]))

        for i = 1:length(cs_delta_group[n])
            out[n][i] = cs_delta_group[n][i][1]
        end
    end

    return out
end

function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function displaypartitionsizes(A::NMRSpectraSimulator.FIDModelType{T,SST}) where {T,SST}

    return collect( length(A.part_inds_compound[k]) for k = 1:length(A.part_inds_compound) )
end

function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

# average Δc vector for each partition element.
function viewaverageΔc(A::NMRSpectraSimulator.FIDModelType{T,SST}) where {T,SST}

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


###### tmp.
function runalignment2(shift_constants,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    LS_inds,
    Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
    As::Vector{NMRSpectraSimulator.CompoundFIDType{T, SST}},
    fs::T,
    SW::T,
    λ_lbs, λ_ubs, κs_β_DOFs, κs_β_orderings;
    optim_algorithm::Symbol = :GN_ESCH,
    max_iters = 5000,
    xtol_rel = 1e-7,
    ftol_rel = 1e-12,
    maxtime = Inf,
    w = ones(T, length(Es)),
    κ_lb_default = 0.2,
    κ_ub_default = 50.0,
    λ_each_lb = 0.7,
    λ_each_ub = 5.0) where {T <: Real, SST}

    @assert length(U_rad_cost) == length(y_cost) == length(LS_inds)

    #
    N_d = sum( getNd(As[n]) for n = 1:length(As) )
    N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_λ = sum( getNλ(As[n]) for n = 1:length(As) )
    # N_vars = N_d + N_β + N_λ

    #### initial values.
    shift_lb = -ones(T, N_d)
    shift_ub = ones(T, N_d)
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
    println("p_initial = ", p_initial)
    println("length(p_initial) = ", length(p_initial))
    println()

    q, updatedfunc, updateβfunc, updateλfunc, updateκfunc, #updatewfunc,
    E_BLS, κ_BLS, b_BLS, getshiftfunc, getβfunc, getλfunc,
    N_vars_set = setupcostcLshiftLS(Es, As, fs, SW,
        LS_inds, U_rad_cost, y_cost, shift_constants;
        w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)
    #updatewfunc(p_initial)
    updateκfunc(p_initial)
    parseκ!(Es, κ_BLS)

    ### cost function.
    obj_func = pp->costcLshift(U_rad_cost, y_cost,
    updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es, E_BLS, κ_BLS, b_BLS, q)

    grad_func = xx->FiniteDiff.finite_difference_gradient(obj_func, xx)

    # force eval for debugging.
    println("Timing: obj_func")
    @time cost_initial = obj_func(p_initial)
    cost_initial = obj_func(p_initial)
    println("obj_func(p_initial) = ", cost_initial)
    println()

    opt = NLopt.Opt(optim_algorithm, length(p_initial)) # global evolutionary.

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
        obj_func, N_vars_set, p_initial
end
