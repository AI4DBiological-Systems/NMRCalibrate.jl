
import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


### user inputs.
# save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
# project_name = "phenylalanine-700"
#
# experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Phenylalanine"
# solvent_ppm_guess = 4.7
# solvent_window_ppm = 0.1

# save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
# project_name = "isoleucine-700"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/L-Isoleucine"
# solvent_ppm_guess = 4.7
# solvent_window_ppm = 0.1



save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
project_name = "glucose-700"
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/D-(+)-Glucose"
solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

# save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
# project_name = "glucose-500-2mM"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-2mM/D-(+)-Glucose"
# solvent_ppm_guess = 4.7
# solvent_window_ppm = 0.1

# save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
# project_name = "glucose-500-0.5mM"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-500-0.5mM/D-(+)-Glucose"
# solvent_ppm_guess = 4.7
# solvent_window_ppm = 0.1

# save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
# project_name = "glucose-nrc-2018"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/src/experiments/NRC/misc/glucose/Sep-25-2018"
# solvent_ppm_guess = 4.7
# solvent_window_ppm = 0.1

# save_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"
# project_name = "Nam-Jan2022"
# solvent_ppm_guess = 4.7
# solvent_window_ppm = 0.1
# experiment_full_path = "/home/roy/MEGAsync/data/NMR/NRC/Nam_4_amino_acid_mixture_Jan_2022/NHK220113_4AA_glu_PBS_d2o/1"

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
ν_0ppm = ν_0ppm,
α_0ppm = α_0ppm,
β_0ppm = β_0ppm,
λ_0ppm = λ_0ppm)


## visualize.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(N, fs, offset_Hz)

S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)



q = uu->NMRDataSetup.evalcomplexLorentzian(uu, α_0ppm, β_0ppm, λ_0ppm, 2*π*ν_0ppm)*fs
q_U = q.(U_y)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(S_U), label = "data")
PyPlot.plot(P_y, real.(q_U), "--", label = "estimated 0ppm resonance component")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("real part of data spectrum (shifted and scaled Discrete Fourier Transform)")



###normalized.

## normalize by max magnitude in the 0 ppm neighbourhood.
# inds = findall(xx->(-0.3<xx<0.3), P_y)
# Z = maximum(abs.(S_U[inds]))/9 # α at 0ppm for DSS is 9 because there are 9 protons.

# normalize by the 0 ppm estimate, if it had an α of 9, which corresponds to the 9 protons for DSS' 0ppm component.
c = NMRDataSetup.evalcomplexLorentzian(ν_0ppm, 9.0, β_0ppm, λ_0ppm, 2*π*ν_0ppm)
Z = q(ν_0ppm)/c
y = S_U ./ Z

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y))

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("real part of normalized data spectrum s.t. 0ppm peak has magnitude 1")
