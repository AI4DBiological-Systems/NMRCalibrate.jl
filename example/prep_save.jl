

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
ν_0ppm = dic[:ν_0ppm]
α_0ppm = dic[:α_0ppm]
β_0ppm = dic[:β_0ppm]
λ_0ppm = dic[:λ_0ppm]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# start the first entry at the frequency corresponding to offset_Hz, which we set to 0.3 ppm.
offset_ppm = 0.3
offset_Hz = ν_0ppm - (ppm2hzfunc(offset_ppm)-ppm2hzfunc(0.0))
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(length(s_t), fs, offset_Hz)

S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

val, ind = findmin( abs.(U_y .- ν_0ppm) )
Z = abs(S_U[ind])

y = S_U ./ Z

### end loading spectrum.
#@assert 1==2

####### compute surrogate.

# TODO document how to use the configuration files. then delete these default values from tutorial.
tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
α_relative_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
Δc_partition_radius = 0.3 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
λ0 = λ_0ppm

Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
κ_λ_lb_default = 0.7 # interpolation lower limit for κ_λ.
κ_λ_ub_default = 1.5 # interpolation upper limit for κ_λ.


# get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
println("Timing: mag equivalence")
@time MEs = NMRSpectraSimulator.getmageqinfomixture(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

println("Timing: setupmixtureproxies()")
@time mixture_params = NMRSpectraSimulator.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ν_0ppm;
    config_path = SH_config_path,
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius)
As = mixture_params


ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

u_offset = 0.5 # in ppm.
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)



println("Timing: fitproxies!()")
dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0)
Bs = NMRSpectraSimulator.fitproxies(As, dummy_SSFID, λ0;
    names = molecule_names,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    κ_λ_lb_default = κ_λ_lb_default,
    κ_λ_ub_default = κ_λ_ub_default,
    u_min = u_min,
    u_max = u_max,
    Δr_default = Δr_default,
    Δκ_λ_default = Δκ_λ_default)

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
ν_0ppm  = ν_0ppm,
λ0 = λ0,
molecule_names = molecule_names,
w = w)
