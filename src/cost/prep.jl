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
    cost_inds = NMRSpecifyRegions.getcostinds(exp_info, P_cost0)

    U_cost = U_cost0[cost_inds]
    P_cost = P_cost0[cost_inds]
    y_cost = y_cost0[cost_inds]

    return Δsys_cs, y_cost, U_cost, P_cost, exp_info, cost_inds
end



# ## TODO I am here
# function runalignment(Δ_shifts::Vector{T},
#     U_cost,
#     y_cost,
#     Es,
#     Bs,
#     mixture_LUT::Vector{Vector{Vector{Int}}},
#     fs::T, SW::T, LS_inds::Vector{Int},
#     w_lb, w_ub;
#     max_iters = 5000,
#     shift_initial = zeros(T, 0),
#     xtol_rel = 1e-7,
#     ftol_rel = 1e-12,
#     maxtime = Inf,
#     warp_param_set::Vector{Piecewise2DLineType{T}} = Vector{Piecewise2DLineType{T}}(undef, 0)) where T <: Real

#     #
#     N_d = sum( length(Bs[n].d) + length(Bs[n].d_singlets) for n = 1:length(Bs) )

#     st_ind = 1
#     updatedfunc = pp->NMRCalibrate.updatemixtured!(Bs, pp, st_ind, fs, SW, Δ_shifts)

#     st_ind_β = N_d + 1
#     updateβfunc = pp->NMRCalibrate.updateβ!(Bs, pp, st_ind_β)
#     #N_β = sum( sum(length(Bs[n].κs_β[l]) for l = 1:length(Bs[n].κs_β)) + length(Bs[n].β_singlets) for n = 1:length(Bs) )
#     N_β = sum( NMRCalibrate.getNβ(Bs[n]) for n = 1:length(Bs) )

#     # I am here. do λupdate.
#     st_ind_λ = st_ind_β + N_β
#     updateλfunc = pp->NMRCalibrate.updateλ!(Bs, pp, st_ind_λ)
#     N_λ = sum( NMRCalibrate.getNλ(Bs[n]) for n = 1:length(Bs) )

#     N_vars = N_d + N_β + N_λ

#     #### initial values.
#     shift_lb = zeros(T, N_d)
#     shift_ub = zeros(T, N_d)
#     shift_initial = zeros(T, N_d)

#     β_lb = ones(T, N_β) .* (-π)
#     β_ub = ones(T, N_β) .* (π)
#     β_initial = zeros(T, N_β)

#     λ_lb = zeros(T, N_λ)
#     λ_ub = zeros(T, N_λ)
#     λ_initial = zeros(T, N_λ)

#     ### initial phase estimate. Constant phase model.
#     #β_est = estimateconstantphase(U_cost, y_cost)
#     #βS_base = initializeFIDparameter(αS, β_est)
#     #fill!(β_initial, β_est)
#     #fill!(β_initial, 0.0)

#     ### set up.
#     p_lb = [ -ones(T, N_spin_groups); β_lb ]
#     p_ub = [ ones(T, N_spin_groups); β_ub ]
#     p_initial = [shift_initial; β_initial]

#     q_star_rad, updatecsfunc, updateβfunc, updatewfunc,
#         w_star, ΩS_star, βS_star,
#         getshiftfunc, getβfunc = setupcostcLshiftLS(αS, ΩS,
#         λS, βS_base, fs, SW, U_cost, LS_inds, y_cost, w_lb, w_ub, Δ_shifts, warp_param_set)
#     updatewfunc(p_initial)

#     q_star = uu->q_star_rad(uu*2*π) # input: Hz.

#     ### cost function.
#     U_cost_rad = U_cost .* (2*π)

#     obj_func = pp->costcLshift(U_cost_rad, y_cost,
#         updatecsfunc, updateβfunc, updatewfunc, pp, q_star_rad)

#     grad_func = xx->FiniteDiff.finite_difference_gradient(obj_func, xx)

#     # force eval for debugging.
#     println("Timing: obj_func")
#     @time cost_initial = obj_func(p_initial)
#     cost_initial = obj_func(p_initial)
#     println("obj_func(p_initial) = ", cost_initial)
#     println()

#     opt = NLopt.Opt(:GN_ESCH, length(p_initial)) # global evolutionary.

#     minf, minx, ret, numevals = runNLopt!(  opt,
#         p_initial,
#         obj_func,
#         grad_func,
#         p_lb,
#         p_ub;
#         max_iters = max_iters,
#         xtol_rel = xtol_rel,
#         ftol_rel = ftol_rel,
#         maxtime = maxtime)

#     println("numevals = ", numevals)
#     println("norm(p_initial-minx) = ", norm(p_initial-minx))

#     println("Timing: obj_func")
#     @time cost_final = obj_func(minx)
#     println("obj_func(minx) = ", cost_final)
#     println()

#     return minx, q_star, w_star, ΩS_star, βS_star, obj_func, p_initial,
#     getshiftfunc, getβfunc
# end