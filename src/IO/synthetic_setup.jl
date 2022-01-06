
Base.@kwdef mutable struct SyntheticCalibrationDataType
    BMRB_setup_label = "BMRB-700" # for use with the NMRData package. See that package's readme for details.

    save_dir = "/home/roy/MEGAsync/outputs/NMR/combined/"
    project_name = "cal01"

    target_names = Vector{String}(undef, 0)

    # mixing weights for each experiment.
    weights = ones(Float64, 0)

    solvent_ppm_guess = 4.7 # guess where the solvent ppm is for each BMRB experiment.
    solvent_window_ppm = 0.1 # allowed solvent shift, in ppm.

end

# taken from an example script in the NMRDataSetup repo.
function setupsyntheticcalibrationdata(config::SyntheticCalibrationDataType)

    ## get BMRB experiment paths from the NMRData package. You can specify your own array of experiments to mix.
    paths = Vector{String}(undef, 0)
    for i = 1:length(config.target_names)
        push!( paths, NMRData.getBMRBexperimentpath(config.BMRB_setup_label, config.target_names[i])[1] )
    end
    @assert all(isdir.(paths)) # check for bad paths.

    ## mix.
    s_t, ν_0ppm, fs, SW,
    s_t_set, ν_0ppm_set, fs_set, SW_set = NMRDataSetup.assemblemixture(paths; weights = config.weights)

    ## store.
    save_folder_path = joinpath(config.save_dir, config.project_name)
    isdir(save_folder_path) || mkdir(save_folder_path); # make save folder if it doesn't exist.

    save_path = joinpath(save_folder_path, "$(config.project_name).bson")
    BSON.bson(save_path,
    s_t = s_t,
    molecule_names = config.target_names,
    fs = fs,
    SW = SW,
    ν_0ppm = ν_0ppm,
    s_t_set = s_t_set,
    ν_0ppm_set = ν_0ppm_set,
    fs_set = fs_set,
    SW_set = SW_set)

    return s_t, ν_0ppm, fs, SW,
        s_t_set, ν_0ppm_set, fs_set, SW_set
end
