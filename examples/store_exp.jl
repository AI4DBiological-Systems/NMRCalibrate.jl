
import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


### user inputs.
save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
project_name = "phenylalanine-700"

experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Phenylalanine"
solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

### end inputs.

## load.
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

