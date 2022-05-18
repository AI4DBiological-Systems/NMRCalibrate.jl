# organize/sort this later.

"""
Info on run_optim usage: minf, minx, ret, N_evals = run_optim(p_β)
"""
function setupcostnestedλd(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
    As::Vector{NMRSpectraSimulator.CompoundFIDType{T, SST}},
    fs::T,
    SW::T,
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    shift_constants;
    w = ones(T, length(Es)),
    optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = 0.2,
    κ_ub_default = 5.0,
    max_iters = 500,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf) where {T <: Real, SST}

    # model.
    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)

    ##### update functions.
    N_d = sum( getNd(As[n]) for n = 1:length(As) )
    N_λ = sum( getNλ(As[n]) for n = 1:length(As) )

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixtured!(As, pp, st_ind_d, fs, SW, shift_constants)

    #λupdate.
    st_ind_λ = fin_ind_d + 1
    fin_ind_λ = st_ind_λ + N_λ -1
    updateλfunc = pp->updateλ!(As, pp, st_ind_λ)

    N_vars_set = [N_d; N_λ]

    # β, κ update.
    run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc,
    q_β = setupβLSsolver(optim_algorithm,
        Es, As, LS_inds, U_rad_cost, y_cost;
        κ_lb_default = κ_lb_default,
        κ_ub_default = κ_ub_default,
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        maxtime = maxtime)

    #### extract parameters from p.
    getshiftfunc = pp->pp[st_ind_d:fin_ind_d]
    getλfunc = pp->pp[st_ind_λ:fin_ind_λ]

    return f, updatedfunc, updateλfunc, getshiftfunc, getλfunc, N_vars_set,
    run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc, q_β
end

# function setupcostnestedλdwarp(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
#     Bs,
#     As,
#     fs::T,
#     SW::T,
#     LS_inds,
#     U_rad_cost,
#     y_cost::Vector{Complex{T}},
#     Δsys_cs::Vector{Vector{T}},
#     itp_a,
#     itp_b;
#     w = ones(T, length(Es)),
#     optim_algorithm = :GN_DIRECT_L,
#     κ_lb_default = 0.2,
#     κ_ub_default = 5.0,
#     max_iters = 500,
#     xtol_rel = 1e-9,
#     ftol_rel = 1e-9,
#     maxtime = Inf) where {T <: Real, SST}
#
#     # model.
#     f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)
#
#     ##### update functions.
#     N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )
#     N_λ = sum( getNλ(Bs[n]) for n = 1:length(Bs) )
#
#     st_ind_d = 1
#     fin_ind_d = st_ind_d + N_d - 1
#     updatedfunc = pp->updatemixturedwarp!(Bs, pp, st_ind_d, fs, SW,
#         Δsys_cs, itp_a, itp_b)
#
#     #λupdate.
#     st_ind_λ = fin_ind_d + 1
#     fin_ind_λ = st_ind_λ + N_λ -1
#     updateλfunc = pp->updateλ!(Bs, pp, st_ind_λ)
#
#     N_vars_set = [N_d; N_λ]
#
#     # β, κ update.
#     run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc,
#     q_β = setupβLSsolver(optim_algorithm,
#         Es, Bs, As, LS_inds, U_rad_cost, y_cost,
#         κs_β_DOFs, κs_β_orderings;
#         κ_lb_default = κ_lb_default,
#         κ_ub_default = κ_ub_default,
#         max_iters = max_iters,
#         xtol_rel = xtol_rel,
#         ftol_rel = ftol_rel,
#         maxtime = maxtime)
#
#     #### extract parameters from p.
#     getshiftfunc = pp->pp[st_ind_d:fin_ind_d]
#     getλfunc = pp->pp[st_ind_λ:fin_ind_λ]
#
#     return f, updatedfunc, updateλfunc, getshiftfunc, getλfunc, N_vars_set,
#     run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc, q_β
# end

function costnestedλd(U,
    S_U::Vector{Complex{T}},
    updatedfunc, updateλfunc,
    p::Vector{T}, Es, As, f,
    run_optim_β_κ::Function,
    E_BLS::Matrix{Complex{T}}, κ_BLS::Vector{T}, b_BLS,
    p_β::Vector{T})::T where T <: Real

    updatedfunc(p)
    updateλfunc(p)

    ### minimize inner problem.
    fill!(p_β, zero(T)) # always start from 0-phase?
    minf, minx, ret, N_evals = run_optim_β_κ(p_β)
    p_β[:] = minx

    # ensure Es/As is updated with the latest β and κ.
    updateβ!(As, minx, 1)
    parseκ!(Es, κ_BLS)

    # evaluate cost.
    cost = norm(reinterpret(T, E_BLS)*κ_BLS - b_BLS)^2

    return cost
end




##### only d.


function setupcostnesteddwarp(Es,
    Bs,
    As,
    fs::T,
    SW::T,
    LS_inds,
    U_rad_cost,
    y_cost::Vector{Complex{T}},
    Δsys_cs::Vector{Vector{T}},
    itp_a,
    itp_b;
    #κs_β_DOFs,
    #κs_β_orderings;
    w = ones(T, length(Es)),
    optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = 0.2,
    κ_ub_default = 5.0,
    max_iters = 500,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf) where T <: Real

    # model.
    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)

    ##### update functions.
    #N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
    N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixturedwarp!(Bs, pp, st_ind_d, fs, SW,
        Δsys_cs, itp_a, itp_b)

    N_vars_set = [N_d; ]

    # β, κ update.
    run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc,
    q_β = setupβLSsolver(optim_algorithm,
        Es, Bs, As, LS_inds, U_rad_cost, y_cost;
        #κs_β_DOFs, κs_β_orderings;
        κ_lb_default = κ_lb_default,
        κ_ub_default = κ_ub_default,
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        maxtime = maxtime)

    #### extract parameters from p.
    getshiftfunc = pp->pp[st_ind_d:fin_ind_d]

    return f, updatedfunc, getshiftfunc, N_vars_set,
    run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc, q_β
end

##### cost.
"""
Note that NLopt handle exceptions gracefully. Check the termination status to see if an exception occurs (would return FORCED_STOP)
"""
function costnestedd(U,
    S_U::Vector{Complex{T}},
    updatedfunc,
    p::Vector{T}, Es, Bs, #κs_β_orderings, κs_β_DOFs, f,
    run_optim_β_κ::Function,
    E_BLS::Matrix{Complex{T}}, κ_BLS::Vector{T}, b_BLS,
    p_β::Vector{T})::T where T <: Real

    updatedfunc(p)

    ### minimize inner problem.
    fill!(p_β, zero(T)) # always start from 0-phase?
    minf, minx, ret, N_evals = run_optim_β_κ(p_β)
    p_β[:] = minx

    # ensure Es/As is updated with the latest β and κ.
    #updateβ!(Bs, κs_β_orderings, κs_β_DOFs, minx, 1)
    updateβ!(Bs, minx, 1)
    parseκ!(Es, κ_BLS)

    # evaluate cost.
    cost = norm(reinterpret(T, E_BLS)*κ_BLS - b_BLS)^2

    return cost
end
