

import NMRDataSetup
import NMRSpectraSimulator
import NMRSpecifyRegions

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JSON

import Statistics
import Random

# include("./helpers/final_helpers.jl")
# include("./helpers/solute_calibrate.jl")
# include("./helpers/loop_entries1.jl")


##### global constants.
#cs_config_path = "/home/roy/Documents/repo/NMRData/src/input/reduced_cs_config.txt"
SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs.json"
surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"
fit_config_path = "/home/roy/Documents/repo/NMRData/input/fit_configs/calibrate_700MHz_type1_select_compounds.json"

# get mapping from molecule names to their spin system info json files.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

# where the experiment bson file is located.
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
experiment_path = joinpath(root_folder, "NRC-glucose-2018") # second input is subfolder name.
load_path = joinpath(experiment_path, "experiment.bson")

# where the saved calibration project folder is for this experiment.
project_name = "NRC-glucose-2018"

save_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
save_folder = joinpath(save_base_folder, project_name)
save_path = joinpath(save_folder, "proxy.bson")
isdir(save_path) || mkpath(save_folder)
#####

##### user input for loading data.

molecule_names = ["D-(+)-Glucose";]
w = [1.0;]

# molecule_names = ["D-(+)-Glucose"; "DSS"; "Ethanol"]
# w = [1.0; 1.0; 1.0]

# molecule_names = ["D-(+)-Glucose"; "DSS"; ]
# w = [1.0; 1.0;]

# note that 0.1% DSS is 0.0046 M = 4.6 mM.

##### end inputs for loading data.

##### user input for surrogate.
dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0) # level 2 model.

##### end inputs for surrogate.

println("Now on $(project_name)")

##### load data that was saved in store_exp.jl. See that script for details of the following commands.
load_path = joinpath(load_path)
dic = BSON.load(load_path)
s_t = dic[:s_t]
fs = dic[:fs]
SW = dic[:SW]
??_0ppm = dic[:??_0ppm]
??_0ppm = dic[:??_0ppm]
??_0ppm = dic[:??_0ppm]
??_0ppm = dic[:??_0ppm]

hz2ppmfunc = uu->(uu - ??_0ppm)*SW/fs
ppm2hzfunc = pp->(??_0ppm + pp*fs/SW)

# start the first entry at the frequency corresponding to offset_Hz, which we set to 0.3 ppm.
offset_ppm = 0.3
offset_Hz = ??_0ppm - (ppm2hzfunc(offset_ppm)-ppm2hzfunc(0.0))
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(length(s_t), fs, offset_Hz)

S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

val, ind = findmin( abs.(U_y .- ??_0ppm) )
Z = abs(S_U[ind])

y = S_U ./ Z

### end loading spectrum.
#@assert 1==2

####### compute surrogate.

# TODO document how to use the configuration files. then delete these default values from tutorial.
tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
??_relative_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
??c_partition_radius = 0.3 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
??0 = ??_0ppm

??r_default = 1.0 # the samples used to build the surrogate is taken every `??r` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
????_??_default = 0.05 # the samples used to build thes urrogate for ??_?? are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
??cs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
??_??_lb_default = 0.7 # interpolation lower limit for ??_??.
??_??_ub_default = 1.5 # interpolation upper limit for ??_??.


# get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
println("Timing: mag equivalence")
@time MEs = NMRSpectraSimulator.getmageqinfomixture(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

println("Timing: setupmixtureproxies()")
@time mixture_params = NMRSpectraSimulator.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ??_0ppm;
    config_path = SH_config_path,
    tol_coherence = tol_coherence,
    ??_relative_threshold = ??_relative_threshold,
    ??c_partition_radius = ??c_partition_radius)
As = mixture_params


??S_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
??S_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(??S_ppm))

u_offset = 0.5 # in ppm.
u_min = ppm2hzfunc(??S_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(??S_ppm_sorted[end] + u_offset)



println("Timing: fitproxies!()")
dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0)
Bs = NMRSpectraSimulator.fitproxies(As, dummy_SSFID, ??0;
    names = molecule_names,
    config_path = surrogate_config_path,
    ??cs_max_scalar_default = ??cs_max_scalar_default,
    ??_??_lb_default = ??_??_lb_default,
    ??_??_ub_default = ??_??_ub_default,
    u_min = u_min,
    u_max = u_max,
    ??r_default = ??r_default,
    ????_??_default = ????_??_default)

####### end mixture proxy.


### save.
BSON.bson(save_path,
As = As,
Bs = Bs,
surrogate_config_path = surrogate_config_path,
SH_config_path = SH_config_path,
u_min = u_min,
u_max = u_max,
y = y,
U_y = U_y,
SW = SW,
fs = fs,
??_0ppm  = ??_0ppm,
??0 = ??0,
molecule_names = molecule_names,
w = w)
