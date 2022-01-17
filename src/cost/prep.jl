function prepareoptim(regions::Vector{RegionInfoType{T}},
    region_select::Int,
    hz2ppmfunc, ppm2hzfunc,
    cost_border_ppm::T,
    αS0::Vector{Vector{Vector{T}}},
    βS0, ΩS0, λS0,
    Δsys_shift_lb, Δsys_shift_ub, sys_w_lb, sys_w_ub,
    Δsys_shift_initial,
    U_cost0, y_cost0) where T <: Real

    ###
    B = regions[region_select]

    cost_start = B.st - cost_border_ppm
    cost_fin = B.fin + cost_border_ppm

    P_cost0 = hz2ppmfunc.(U_cost0)
    band_inds = filterfreqpositions(P_cost0, [cost_start;], [cost_fin;])
    U_cost = U_cost0[band_inds]
    y_cost = y_cost0[band_inds]
    P_cost = P_cost0[band_inds]

    ppm_lower = P_cost[1]
    ppm_upper = P_cost[end]

    ###
    αS = assembleFIDbounds(B.spin_group_indices, αS0)
    ΩS = assembleFIDbounds(B.spin_group_indices, ΩS0)
    λS = assembleFIDbounds(B.spin_group_indices, λS0)
    βS = assembleFIDbounds(B.spin_group_indices, βS0)

    shift_lb, shift_ub, w_lb, w_ub,
        shift_initial = assembleshiftwbounds(B.spin_group_indices,
        Δsys_shift_lb, Δsys_shift_ub, sys_w_lb, sys_w_ub,
        Δsys_shift_initial)

    return shift_lb, shift_ub, shift_initial, w_lb, w_ub, αS, ΩS, λS, βS,
    y_cost, U_cost
end