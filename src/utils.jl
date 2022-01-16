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

function findfreqrange(As::Vector{NMRSpectraSimulator.CompoundFIDType{T}}, hz2ppmfunc) where T

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))
    
    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) )
        
        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end