
#### κ parsing.


function countκ(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}}) where {T,SST}

    N_κ = 0
    N_κ_singlets = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κs_α)
            #for l = 1:length(Es[n].κs_α[i])
            N_κ += length(Es[n].κs_α[i])

            #end
        end

        N_κ_singlets += length(Es[n].κs_α_singlets)
    end

    return N_κ, N_κ_singlets
end

function parseκ!(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}},
    κ_vec::Vector{T}) where {T,SST}

    @assert length(κ_vec) == sum(countκ(Es))

    j = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κs_α)
            for l = 1:length(Es[n].κs_α[i])
                j += 1
                Es[n].κs_α[i][l] = κ_vec[j]
            end
        end

        for i = 1:length(Es[n].κs_α_singlets)
            j += 1
            Es[n].κs_α_singlets[i] = κ_vec[j]
        end
    end

    return j
end

function resetκ!(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T,SST}}) where {T,SST}

    j = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κs_α)
            for l = 1:length(Es[n].κs_α[i])
                j += 1
                Es[n].κs_α[i][l] = one(T)
            end
        end

        for i = 1:length(Es[n].κs_α_singlets)
            j += 1
            Es[n].κs_α_singlets[i] = one(T)
        end
    end

    return j
end

#### z parsing.

function countz(Gs::Vector{NMRSpectraSimulator.zCompoundFIDType{T,SST}}) where {T,SST}

    N_z = 0
    N_z_singlets = 0
    for n = 1:length(Gs)
        for i = 1:length(Gs[n].zs)
            #for l = 1:length(Gs[n].zs[i])
            N_z += length(Gs[n].zs[i])

            #end
        end

        N_z_singlets += length(Gs[n].zs_singlets)
    end

    return N_z, N_z_singlets
end

function parsez!(Gs::Vector{NMRSpectraSimulator.zCompoundFIDType{T,SST}},
    z_vec::Vector{Complex{T}}) where {T,SST}

    @assert length(z_vec) == sum(countz(Gs))

    j = 0
    for n = 1:length(Gs)
        for i = 1:length(Gs[n].zs)
            for l = 1:length(Gs[n].zs[i])
                j += 1
                Gs[n].zs[i][l] = z_vec[j]
            end
        end

        for i = 1:length(Gs[n].zs_singlets)
            j += 1
            Gs[n].zs_singlets[i] = z_vec[j]
        end
    end

    return j
end

function resetz!(Gs::Vector{NMRSpectraSimulator.zCompoundFIDType{T,SST}}) where {T,SST}

    j = 0
    for n = 1:length(Gs)
        for i = 1:length(Gs[n].zs)
            for l = 1:length(Gs[n].zs[i])
                j += 1
                Gs[n].zs[i][l] = one(T)
            end
        end

        for i = 1:length(Gs[n].zs_singlets)
            j += 1
            Gs[n].zs_singlets[i] = one(T)
        end
    end

    return j
end
