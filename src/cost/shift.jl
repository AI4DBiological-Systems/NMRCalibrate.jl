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

# WIP.
function setupcostalign0(Es::Vector{NMRSpectraSimulator.καFIDModelType{T, SST}},
    Bs::Vector{NMRSpectraSimulator.FIDModelType{T, SST}},
    As::Vector{NMRSpectraSimulator.SHType{T}},
    κs_β_orderings::Vector{Vector{Vector{Int}}},
    κs_β_DOFs::Vector{Vector{Int}},
    fs::T,
    SW::T,
    LS_inds,
    U0_rad,
    y0::Vector{Complex{T}},
    Δsys_cs::Vector{Vector{T}};
    w = ones(T, length(Es))) where {T <: Real, SST}

    # model.
    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Bs, Es; w = w)

    ##### update functions.
    N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )
    N_β = sum( getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixtured!(Bs,
        p, st_ind, fs, SW, Δsys_cs)

    st_ind_β = fin_ind_d + 1
    fin_ind_β = st_ind_β + N_β - 1
    updateβfunc = pp->updateβ!(Bs, pp, st_ind_β)

    N_vars_set = [N_d; N_β]

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
