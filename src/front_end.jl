function trydiffΔcradius(Δc_partition_radius_candidates::Vector{T},
    molecule_names, base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc,
    fs, SW, λ0, ν_0ppm, early_exit_part_size, Δcs_max, tol_coherence,
    α_relative_threshold,
    dummy_SSFID::SST)::Tuple{Vector{NMRSpectraSimulator.CompoundFIDType{T,SST}}, T} where {T <: Real, SST}

    @assert early_exit_part_size > 0
    @assert all(Δc_partition_radius_candidates .> zero(T))

    Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

    for Δc_partition_radius in Δc_partition_radius_candidates[1:end-1]

        mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
            base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW,
            λ0, ν_0ppm, dummy_SSFID;
            tol_coherence = tol_coherence,
            α_relative_threshold = α_relative_threshold,
            Δc_partition_radius = Δc_partition_radius)
        As = mixture_params

        if all( all(NMRCalibrate.displaypartitionsizes(As[n]) .<= early_exit_part_size) for n = 1:length(As) )
            return mixture_params, Δc_partition_radius
        end
    end

    mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm, dummy_SSFID;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius_candidates[end])

    return mixture_params, Δc_partition_radius_candidates[end]
end
