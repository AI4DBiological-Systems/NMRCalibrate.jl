# organize/sort this later.

function setupcostnestedλd(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T, SST}},
    As::Vector{NMRSpectraSimulator.CompoundFIDType{T, SST}},
    fs::T,
    SW::T,
    LS_inds,
    U0,
    y0::Vector{Complex{T}},
    shift_constants;
    w = ones(T, length(Es))) where {T <: Real, SST}

    # model.
    f = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)

    ##### update functions.
    N_d = sum( getNd(As[n]) for n = 1:length(As) )
    N_λ = sum( getNλ(As[n]) for n = 1:length(As) )

    st_ind_d = 1
    fin_ind_d = st_ind_d + N_d - 1
    updatedfunc = pp->updatemixtured!(As, pp, st_ind_d, fs, SW, shift_constants)

    #λupdate.
    st_ind_λ = fin_ind_β + 1
    fin_ind_λ = st_ind_λ + N_λ -1
    updateλfunc = pp->updateλ!(As, pp, st_ind_λ)

    N_vars_set = [N_d; N_λ]

    #### extract parameters from p.
    getshiftfunc = pp->pp[st_ind_d:fin_ind_d]
    getλfunc = pp->pp[st_ind_λ:fin_ind_λ]

    return f, updatedfunc, updateλfunc, getshiftfunc, getλfunc, N_vars_set
end

function costnestedλd(U,
    S_U::Vector{Complex{T}},
    updatedfunc, updateλfunc,
    p::Vector{T}, Es, q,
    run_inner_optim::Function,
    κ_BLS::Vector{T},
    p_β::Vector{T})::T where T <: Real

    updatedfunc(p)
    updateλfunc(p)

    ### minimize inner problem.
    minf, minx, ret, numevals = run_inner_optim(p_β)
    p_β[:] = minx

    # ensure Es is updated with the latest.
    updateβfunc(minx)
    updateκfunc(minx)
    parseκ!(Es, κ_BLS)

    ### l-2 costfunc for outer problem.
    cost = zero(T)
    for m = 1:length(S_U)

        ## eval.
        q_u = q(U[m])
        cost += abs2( q_u - S_U[m] )
    end

    return cost
end
