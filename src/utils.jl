function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

# average Δc vector for each partition element.
function viewaverageΔc(A::NMRSpectraSimulator.CompoundFIDType{T}) where T
    
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

function countκ(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T}}) where T
    
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

function parseκ!(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T}},
    κ_vec::Vector{T}) where T

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

    return nothing
end

function resetκ!(Es::Vector{NMRSpectraSimulator.κCompoundFIDType{T}}) where T

    j = 0
    for n = 1:length(Es)
        for i = 1:length(Es[n].κ)
            for l = 1:length(Es[n].κ[i])
                j += 1
                Es[n].κ[i][l] = one(T)
            end
        end
    end

    return nothing
end