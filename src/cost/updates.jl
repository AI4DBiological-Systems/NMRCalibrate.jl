function updatemixtured!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    #warp_param_set::Vector{Piecewise2DLineType{T}},
    Δ_shifts::Vector{T})::Int where T <: Real

    #@assert length(warp_param_set) == length(Δ_shifts)

    j = st_ind - 1

    for n = 1:length(As)
        for i = 1:length(As[n].d)
            j += 1

            #p2 = p[j]*0.05 #*Δ_shifts[j]
            p2 = p[j]*Δ_shifts[j]
            As[n].d[i] = convertΔcstoΔω0(p2, fs, SW)
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


function updateβ!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T}},
    p::Vector{T},
    st_ind::Int)::Int where T <: Real

    j = st_ind - 1

    for n = 1:length(As)
        for i = 1:length(As[n].κs_β)
            for l = 1:length(As[n].κs_β[i])
                j += 1

                As[n].κs_β[i][l] = p[j]
            end
        end

        for i = 1:length(As[n].β_singlets)
            j += 1

            As[n].β_singlets[i] = p[j]
        end
    end

    return j
end

function getNβ(A::NMRSpectraSimulator.CompoundFIDType{T}) where T

    counter_sys = 0
    for i = 1:length(A.κs_β)
        counter_sys += length(A.κs_β[i])
    end

    N_β = counter_sys + length(A.β_singlets)

    return N_β
end