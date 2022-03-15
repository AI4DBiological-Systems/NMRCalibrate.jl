function updatemixtured!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}},
    p::Vector{T},
    st_ind::Int,
    fs::T,
    SW::T,
    #warp_param_set::Vector{Piecewise2DLineType{T}},
    Δ_shifts::Vector{T})::Int where {T <: Real, SST}

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

function updateλ!(As::Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}},
    p::Vector{T},
    st_ind::Int)::Int where {T <: Real, SST}

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


function getNλ(A::NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType1{T}}) where T

    counter_sys = 0
    for i = 1:length(A.ss_params.κs_λ)
        counter_sys += length(A.ss_params.κs_λ[i])
    end

    return counter_sys + length(A.κs_λ_singlets)
end

function getNλ(A::NMRSpectraSimulator.CompoundFIDType{T,NMRSpectraSimulator.SpinSysFIDType2{T}}) where T

    counter_sys = 0
    for i = 1:length(A.ss_params.κs_λ)
        for k = 1:length(A.ss_params.κs_λ[i])
            counter_sys += length(A.ss_params.κs_λ[i])
        end
    end

    return counter_sys + length(A.κs_λ_singlets)
end

#### κ parsing.


function countκ(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}}) where {T,SST}

    N_κ = 0
    N_κ_singlets = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κ)
            #for l = 1:length(Es[n].κ[i])
            N_κ += length(Es[n].κ[i])

            #end
        end

        N_κ_singlets += length(Es[n].κ_singlets)
    end

    return N_κ, N_κ_singlets
end

function parseκ!(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}},
    κ_vec::Vector{T}) where {T,SST}

    @assert length(κ_vec) == sum(countκ(Es))

    j = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κ)
            for l = 1:length(Es[n].κ[i])
                j += 1
                Es[n].κ[i][l] = κ_vec[j]
            end
        end

        for i = 1:length(Es[n].κ_singlets)
            j += 1
            Es[n].κ_singlets[i] = κ_vec[j]
        end
    end

    return j
end

function resetκ!(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}}) where {T,SST}

    j = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κ)
            for l = 1:length(Es[n].κ[i])
                j += 1
                Es[n].κ[i][l] = one(T)
            end
        end

        for i = 1:length(Es[n].κ_singlets)
            j += 1
            Es[n].κ_singlets[i] = one(T)
        end
    end

    return j
end
