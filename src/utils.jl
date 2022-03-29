
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function displaypartitionsizes(A::NMRSpectraSimulator.CompoundFIDType{T,SST}) where {T,SST}

    return collect( length(A.part_inds_compound[k]) for k = 1:length(A.part_inds_compound) )
end

function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

# average Δc vector for each partition element.
function viewaverageΔc(A::NMRSpectraSimulator.CompoundFIDType{T,SST}) where {T,SST}

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

function findfreqrange(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}}, hz2ppmfunc) where {T,SST}

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
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
