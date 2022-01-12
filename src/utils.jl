function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end


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