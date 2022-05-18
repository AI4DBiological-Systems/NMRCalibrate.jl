
function applywarptoshifts(p::Vector{T},
    Bs,
    st_ind::Int,
    itp_a,
    itp_b)::Vector{T} where T <: Real

    x = similar(p)
    applywarptoshifts!(x, Bs, p, st_ind, itp_a, itp_b)

    return x
end

function applywarptoshifts!(x::Vector{T},
    Bs,
    p::Vector{T},
    st_ind::Int,
    itp_a,
    itp_b) where T <: Real

    j = st_ind - 1

    for n = 1:length(Bs)

        # first.
        i = 1
        j += 1
        x[j] = p[j]

        # itp.
        target = convertcompactdomain(p[j], -one(T), one(T), zero(T), one(T))
        a = itp_a(target)
        b = itp_b(target)

        for i = 2:length(Bs[n].ss_params.d)
            j += 1

            x[j] = MonotoneMaps.evalcompositelogisticprobit(p[j], a, b, -one(T), one(T))
        end

        for i = 1:length(Bs[n].d_singlets)
            j += 1

            x[j] = p[j]
        end
    end

    return j
end

# assume shifts always between [-1,1]
function updatemixturedwarp!(Bs::Vector{NMRSpectraSimulator.FIDModelType{T,SST}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    #warp_param_set::Vector{Piecewise2DLineType{T}},
    Δsys_cs::Vector{Vector{T}},
    itp_a,
    itp_b)::Int where {T <: Real, SST}

    #@assert length(warp_param_set) == length(Δ_shifts)

    j1 = applywarptoshifts!(p, Bs, p, st_ind, itp_a, itp_b)
    j = updatemixtured!(Bs, p, st_ind, fs, SW, Δsys_cs)

    @assert j1 == j # debug.

    return j
end

function updatemixtured!(As::Vector{NMRSpectraSimulator.FIDModelType{T,NMRSpectraSimulator.SpinSysParamsType1{T}}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    Δsys_cs::Vector{Vector{T}})::Int where T <: Real

    #@assert length(warp_param_set) == length(Δ_shifts)

    j = st_ind - 1

    for n = 1:length(As)

        N_spins_sys = length(As[n].ss_params.d)
        @assert length(Δsys_cs[n]) == N_spins_sys + length(As[n].d_singlets)

        for i = 1:N_spins_sys
            j += 1

            #p2 = p[j]*0.05 #*Δ_shifts[j]
            p2 = p[j]*Δsys_cs[n][i]
            As[n].ss_params.d[i] = convertΔcstoΔω0(p2, fs, SW)
        end

        for i = 1:length(As[n].d_singlets)
            j += 1

            p2 = p[j]*Δsys_cs[n][i+N_spins_sys]
            As[n].d_singlets[i] = convertΔcstoΔω0(p2, fs, SW)
        end
    end

    return j
end

function convertΔcstoΔω0(x::T, fs::T, SW::T)::T where T
    return x*2*π*fs/SW
end

function updatemixtured!(As::Vector{NMRSpectraSimulator.FIDModelType{T,NMRSpectraSimulator.SpinSysParamsType2{T}}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    #warp_param_set::Vector{Piecewise2DLineType{T}},
    shift_constants::Tuple{Vector{Vector{T}},T})::Int where T <: Real

    Δsys_cs, γ = shift_constants
    @assert length(Δsys_cs) == length(As)

    j = st_ind - 1

    d_i::T = NaN

    for n = 1:length(As)

        for i = 1:length(As[n].ss_params.d)

            # first.
            j += 1
            p2 = p[j]*Δsys_cs[n][i]
            d_i = convertΔcstoΔω0(p2, fs, SW)

            As[n].ss_params.d[i][1] = d_i

            # rest.
            for k = 2:length(As[n].ss_params.d[i])

                j += 1
                As[n].ss_params.d[i][k] = p[j]*γ + d_i
            end
        end

        for i = 1:length(As[n].d_singlets)
            j += 1

            p2 = p[j]*Δ_shifts[j]
            As[n].d_singlets[i] = convertΔcstoΔω0(p2, fs, SW)
        end
    end

    return j
end


function getNd(A)::Int

    counter_sys = 0
    for i = 1:length(A.ss_params.d)
        counter_sys += length(A.ss_params.d[i])
    end

    return counter_sys + length(A.d_singlets)
end

# function getNd(A::NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType2{T}}) where T
#
#     counter_sys = 0
#     for i = 1:length(A.ss_params.d)
#         for k = 1:length(A.ss_params.d[i])
#             counter_sys += length(A.ss_params.d[i][k])
#         end
#     end
#
#     return counter_sys + length(A.d_singlets)
# end

function updateβ!(Bs::Vector{NMRSpectraSimulator.FIDModelType{T,SST}},
    κs_β_orderings::Vector{Vector{Vector{Int}}},
    κs_β_DOFs::Vector{Vector{Int}},
    p::Vector{T},
    st_ind::Int)::Int where {T <: Real, SST}

    j = st_ind - 1

    for n = 1:length(Bs)
        for i = 1:length(Bs[n].ss_params.κs_β)

            @assert length(κs_β_orderings[n][i]) == length(Bs[n].ss_params.κs_β[i])

            for l = 1:length(Bs[n].ss_params.κs_β[i])

                idx = κs_β_orderings[n][i][l] + j # use the idx-th unique equivalent value.

                Bs[n].ss_params.κs_β[i][l] = p[idx]
            end
            j += κs_β_DOFs[n][i] # add the number of unique values for this spin system we just processed.

        end

        for i = 1:length(Bs[n].β_singlets)
            j += 1

            Bs[n].β_singlets[i] = p[j]
        end
    end

    return j
end

function updateβ!(Bs::Vector{NMRSpectraSimulator.FIDModelType{T,SST}},
    p::Vector{T},
    st_ind::Int)::Int where {T <: Real, SST}

    j = st_ind - 1

    for n = 1:length(Bs)
        for i = 1:length(Bs[n].ss_params.κs_β)
            for l = 1:length(Bs[n].ss_params.κs_β[i])

                j += 1
                Bs[n].ss_params.κs_β[i][l] = p[j]
            end
        end

        for i = 1:length(Bs[n].β_singlets)
            j += 1
            Bs[n].β_singlets[i] = p[j]
        end
    end

    return j
end

function getNβ(A::NMRSpectraSimulator.FIDModelType{T,SST}) where {T,SST}

    counter_sys = 0
    for i = 1:length(A.ss_params.κs_β)
        counter_sys += length(A.ss_params.κs_β[i])
    end

    return counter_sys + length(A.β_singlets)
end

function getNβ(κs_β_DOF::Vector{Int}, B::NMRSpectraSimulator.FIDModelType{T,SST}) where {T,SST}

    counter_sys = sum( κs_β_DOF[i] for i = 1:length(κs_β_DOF) )

    return counter_sys + length(B.β_singlets)
end

function updateλ!(Bs::Vector{NMRSpectraSimulator.FIDModelType{T,NMRSpectraSimulator.SpinSysParamsType1{T}}},
    p::Vector{T},
    st_ind::Int)::Int where T <: Real

    j = st_ind - 1

    for n = 1:length(Bs)
        for i = 1:length(Bs[n].ss_params.κs_λ)

            j += 1
            Bs[n].ss_params.κs_λ[i] = p[j]
        end

        for i = 1:length(Bs[n].κs_λ_singlets)

            j += 1
            Bs[n].κs_λ_singlets[i] = p[j]
        end
    end

    return j
end

"""
use common λ for type2, just to get code going.
"""
function updateλ!(As::Vector{NMRSpectraSimulator.FIDModelType{T,NMRSpectraSimulator.SpinSysParamsType2{T}}},
    p::Vector{T},
    st_ind::Int)::Int where T <: Real

    j = st_ind - 1

    κ_λ_i::T = NaN

    for n = 1:length(As)

        for i = 1:length(As[n].ss_params.κs_λ)

            # first.
            j += 1
            κ_λ_i = p[j]

            As[n].ss_params.κs_λ[i][1] = κ_λ_i

            # rest.
            for k = 2:length(As[n].ss_params.κs_λ[i])

                ### uncommon λ.
                #j += 1
                #As[n].ss_params.κs_λ[i][k] = κ_λ_i + p[j]*γ_λ
                ### end uncommon λ.

                ### common λ
                As[n].ss_params.κs_λ[i][k] = κ_λ_i
                ### end common λ
            end
        end

        for i = 1:length(As[n].κs_λ_singlets)

            j += 1
            As[n].κs_λ_singlets[i] = p[j]
        end
    end

    return j
end

function getNλ(A::NMRSpectraSimulator.FIDModelType{T,NMRSpectraSimulator.SpinSysParamsType1{T}}) where T

    counter_sys = 0
    for i = 1:length(A.ss_params.κs_λ)
        counter_sys += length(A.ss_params.κs_λ[i])
    end

    return counter_sys + length(A.κs_λ_singlets)
end

function getNλ(A::NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType2{T}}) where T

    counter_sys = 0

    ### uncommon λ
    # for i = 1:length(A.ss_params.κs_λ)
    #     for k = 1:length(A.ss_params.κs_λ[i])
    #         counter_sys += length(A.ss_params.κs_λ[i])
    #     end
    # end

    ### common λ.
    for i = 1:length(A.ss_params.κs_λ)
        counter_sys += length(A.ss_params.κs_λ[i])
    end

    return counter_sys + length(A.κs_λ_singlets)
end
