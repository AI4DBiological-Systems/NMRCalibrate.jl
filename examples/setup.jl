# load an experiment, save as JLD.
# load a proxy, save as JLD.import NMRData
# mixture of synthetics.


#include("../src/NMRCalibrate.jl")
using NMRCalibrate

#import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

## user inputs.
BMRB_setup_label = "BMRB-700" # for use with the NMRData package. See that package's readme for details.

save_dir = "/home/roy/MEGAsync/outputs/NMR/combined/"
project_name = "cal01"

target_names = ["D-(+)-Glucose"; # compounds you wish to add to the mixture.
"L-Alanine"]


weights = ones(Float64, length(target_names)) # mixing weights for each experiment.

solvent_ppm_guess = 4.7 # guess where the solvent ppm is for each BMRB experiment.
solvent_window_ppm = 0.1 # allowed solvent shift, in ppm.
## end inputs.

config = NMRCalibrate.SyntheticCalibrationDataType(
BMRB_setup_label = BMRB_setup_label,
save_dir = save_dir,
project_name = project_name,
target_names = target_names,
weights = weights,
solvent_ppm_guess = solvent_ppm_guess,
solvent_window_ppm = solvent_window_ppm)

### end inputs.

s_t, ν_0ppm, fs, SW, s_t_set, ν_0ppm_set, fs_set,
SW_set = NMRCalibrate.setupsyntheticcalibrationdata(config)

#### plot.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U, U_inds = getwraparoundDFTfreqs(N, fs, offset_Hz)
y = DFT_s[U_inds]

##############################

S_U = DFT_s[U_inds]
P = hz2ppmfunc.(U)

S_set_U = collect( (fft(s_t_set[n]) ./fs)[U_inds] for n = 1:length(s_t_set))
#P_cost = hz2ppmfunc.(U_cost)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(S_U), label = "combined spectrum")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("simulated spectra")
