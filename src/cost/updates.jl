
# assume shifts always between [-1,1]
function updatemixturedwarp!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType1{T}}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    #warp_param_set::Vector{Piecewise2DLineType{T}},
    Δ_shifts::Vector{T},
    itp_a,
    itp_b)::Int where T <: Real

    #@assert length(warp_param_set) == length(Δ_shifts)

    j = st_ind - 1

    for n = 1:length(As)

        # first.
        i = 1
        j += 1
        p2 = p[j]*Δ_shifts[j]
        As[n].ss_params.d[i] = convertΔcstoΔω0(p2, fs, SW)

        # itp.
        target = convertcompactdomain(p[j], -one(T), one(T), zero(T), one(T))
        a = itp_a(target)
        b = itp_b(target)

        for i = 2:length(As[n].ss_params.d)
            j += 1

            p_j = MonotoneMaps.evalcompositelogisticprobit(p[j], a, b)
            p2 = p[j]*Δ_shifts[j]
            As[n].ss_params.d[i] = convertΔcstoΔω0(p2, fs, SW)
        end

        for i = 1:length(As[n].d_singlets)
            j += 1

            p2 = p[j]*Δ_shifts[j]
            As[n].d_singlets[i] = convertΔcstoΔω0(p2, fs, SW)
        end
    end

    return j
end

function updatemixtured!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType1{T}}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    #warp_param_set::Vector{Piecewise2DLineType{T}},
    Δ_shifts::Vector{T})::Int where T <: Real

    #@assert length(warp_param_set) == length(Δ_shifts)

    j = st_ind - 1

    for n = 1:length(As)
        for i = 1:length(As[n].ss_params.d)
            j += 1

            #p2 = p[j]*0.05 #*Δ_shifts[j]
            p2 = p[j]*Δ_shifts[j]
            As[n].ss_params.d[i] = convertΔcstoΔω0(p2, fs, SW)
        end

        for i = 1:length(As[n].d_singlets)
            j += 1

            p2 = p[j]*Δ_shifts[j]
            As[n].d_singlets[i] = convertΔcstoΔω0(p2, fs, SW)
        end
    end

    return j
end

function convertΔcstoΔω0(x::T, fs::T, SW::T)::T where T
    return x*2*π*fs/SW
end

function updatemixtured!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType2{T}}},
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

function updateβ!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}},
    p::Vector{T},
    st_ind::Int)::Int where {T <: Real, SST}

    j = st_ind - 1

    for n = 1:length(As)
        for i = 1:length(As[n].ss_params.κs_β)
            for l = 1:length(As[n].ss_params.κs_β[i])
                j += 1

                As[n].ss_params.κs_β[i][l] = p[j]
            end
        end

        for i = 1:length(As[n].β_singlets)
            j += 1

            As[n].β_singlets[i] = p[j]
        end
    end

    return j
end

function getNβ(A::NMRSpectraSimulator.CompoundFIDType{T,SST}) where {T,SST}

    counter_sys = 0
    for i = 1:length(A.ss_params.κs_β)
        counter_sys += length(A.ss_params.κs_β[i])
    end

    return counter_sys + length(A.β_singlets)
end

function updateλ!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType1{T}}},
    p::Vector{T},
    st_ind::Int)::Int where T <: Real

    j = st_ind - 1

    for n = 1:length(As)
        for i = 1:length(As[n].ss_params.κs_λ)

            j += 1
            As[n].ss_params.κs_λ[i] = p[j]
        end

        for i = 1:length(As[n].κs_λ_singlets)

            j += 1
            As[n].κs_λ_singlets[i] = p[j]
        end
    end

    return j
end

"""
use common λ for type2, just to get code going.
"""
function updateλ!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType2{T}}},
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

function getNλ(A::NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType1{T}}) where T

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
