
import NMRDataSetup

using FFTW
import PyPlot
import BSON

import Base.@kwdef

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

@kwdef mutable struct StoreConfigType
    save_dir::String
    project_name::String
    experiment_full_path::String
    solvent_ppm_guess::Float64 = 4.7
    solvent_window_ppm::Float64 = 0.1
end
# function StoreConfigType(x,y,z)
#     return StoreConfigType(x,y,z,4.7,0.1)
# end

metadata_set = Vector{StoreConfigType}(undef, 0)

### user inputs.
save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final"

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Alanine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Alanine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Arginine-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Arginine"))

# missing aspartic.

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Cysteine-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Cysteine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Glutamic acid-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Glutamic acid"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Glutamine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Glutamine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "Glycine-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/Glycine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Histidine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Histidine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Isoleucine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Isoleucine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Leucine-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Leucine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Lysine-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Lysine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Methionine-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Methionine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Phenylalanine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Phenylalanine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Proline-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Proline"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Serine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Serine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Threonine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Threonine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Tryptophan-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Tryptophan"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Tyrosine-500-2mM",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/L-Tyrosine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "L-Valine-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Valine"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "D-(+)-Glucose-700",
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/D-(+)-Glucose"))

push!(metadata_set, StoreConfigType( save_dir = save_dir,
project_name = "D-(+)-Glucose-NRC-600",
experiment_full_path = "/home/roy/Documents/repo/NMRData/src/experiments/NRC/misc/glucose/Sep-25-2018"))

### end inputs.

## load.
function batchstoreexp(metadata_set)

    for i = 1:length(metadata_set)
        save_dir = metadata_set[i].save_dir
        project_name = metadata_set[i].project_name
        experiment_full_path = metadata_set[i].experiment_full_path
        solvent_ppm_guess = metadata_set[i].solvent_ppm_guess
        solvent_window_ppm = metadata_set[i].solvent_window_ppm

        s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW, α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
            results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
            results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
            solvent_ppm = solvent_ppm_guess,
            solvent_window_ppm = solvent_window_ppm)

        ## store.
        save_folder_path = joinpath(save_dir, project_name)
        isdir(save_folder_path) || mkpath(save_folder_path); # make save folder if it doesn't exist.

        save_path = joinpath(save_folder_path, "$(project_name).bson")
        BSON.bson(save_path,
        s_t = s_t,
        fs = fs,
        SW = SW,
        ν_0ppm = ν_0ppm)
    end
end

batchstoreexp(metadata_set)
