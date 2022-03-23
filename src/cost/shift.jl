function costcLshift(U_rad,
    S_U::Vector{Complex{T}},
    updatedfunc, updateβfunc, updateλfunc, updateκfunc,
    p::Vector{T},
    Es,
    E_BLS::Matrix{Complex{T}}, κ_BLS::Vector{T}, b_BLS,
    f)::T where T <: Real

    #println("p = ", p)
    ### update main optim variables.
    updatedfunc(p)
    updateβfunc(p)
    updateλfunc(p)

    updateκfunc(p)
    parseκ!(Es, κ_BLS)

    # ## l-2 costfunc.
    # cost = zero(T)
    # for m = 1:length(S_U)
    #
    #     f_u = f(U_rad[m])
    #
    #     cost += abs2( f_u - S_U[m] )
    #     #cost += (abs2(f_u) - abs2(S_U[m]))^2
    #     #cost += abs2( f_buffer[m] - DFT_s[m] ) + (abs2(f_buffer[m]) - abs2(DFT_s[m]))^2
    # end
    B = reinterpret(T, E_BLS)
    tmp = B*κ_BLS - b_BLS
    cost = dot(tmp, tmp)

    return cost
end

function setupcostcLshiftLS(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
    As::Vector{NMRSpectraSimulator.CompoundFIDType{T, SST}},
    fs::T,
    SW::T,
    LS_inds,
    U0_rad,
    y0::Vector{Complex{T}},
    shift_constants;
    w = ones(T, length(Es)),
    κ_lb_default = 0.2,
    κ_ub_default = 5.0) where {T <: Real, SST}

    #N = length(ΩS0)
    #@assert length(ΩS0) == length(αS) == length(w_lb) == length(w_ub)

    ## allocate buffers.
    #U_rad_LS = U0[LS_inds] .* (2*π)
    #A, w_BLS = setupLSwcL(U_rad_LS, αS)
    #
    #tmp = S_U0[LS_inds]
    #LS_b = [real.(tmp); imag(tmp)]

    # model.
    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)

    ##### update functions.
    N_d = sum( getNd(As[n]) for n = 1:length(As) )
    N_β = sum( getNβ(As[n]) for n = 1:length(As) )
    N_λ = sum( getNλ(As[n]) for n = 1:length(As) )

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixtured!(As, pp, st_ind_d, fs, SW, shift_constants)

    st_ind_β = fin_ind_d + 1
    fin_ind_β = st_ind_β + N_β - 1
    updateβfunc = pp->updateβ!(As, pp, st_ind_β)

    #λupdate.
    st_ind_λ = fin_ind_β + 1
    fin_ind_λ = st_ind_λ + N_λ -1
    updateλfunc = pp->updateλ!(As, pp, st_ind_λ)

    N_vars_set = [N_d; N_β; N_λ]

    ### LS κ.
    U_rad_LS = U0_rad[LS_inds]

    N_κ, N_κ_singlets = countκ(Es)
    N_κ_vars = N_κ + N_κ_singlets
    E_BLS, κ_BLS, b_BLS = setupupdateLS(length(U_rad_LS), N_κ_vars, y0[LS_inds])

    κ_lb = ones(N_κ_vars) .* κ_lb_default
    κ_ub = ones(N_κ_vars) .* κ_ub_default

    updateκfunc = xx->updateκ!(E_BLS, b_BLS, κ_BLS,
    U_rad_LS, Es, w, κ_lb, κ_ub)

    #### extract parameters from p.
    getshiftfunc = pp->pp[st_ind_d:fin_ind_d]
    getβfunc = pp->pp[st_ind_β:fin_ind_β]
    getλfunc = pp->pp[st_ind_λ:fin_ind_λ]

    return f, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
    E_BLS, κ_BLS, b_BLS, getshiftfunc, getβfunc, getλfunc, N_vars_set
end
