

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
mkpath(save_folder)
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
offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))
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

tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
α_relative_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
Δc_partition_radius = 0.3 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
λ0 = λ_0ppm

Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
Δc_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
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

# We only work with a single compound in this tutorial. Assign a new object for this compound to reduce clutter.
A = As[1];

ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)



println("Timing: fitproxies!()")
dummy_SSFID = NMRSpectraSimulator.SpinSysParamsType1(0.0)
Bs = NMRSpectraSimulator.fitproxies(As, dummy_SSFID, λ0;
    names = molecule_names,
    config_path = surrogate_config_path,
    Δc_max_scalar_default = Δc_max_scalar_default,
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


@assert 1==209999

####### cost func.

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)

##### set up updates.

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)
#ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )



## parameters that affect qs.
# A.ss_params.d, A.ss_params.κs_λ, A.ss_params.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
# purposely perturb κ.

Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )


κ_lb_default = 0.2
κ_ub_default = 50.0

# TODO here, change shift delta cs to each group.


# default, type2.
γ = 0.001 # in ppm.
shift_constants = (Δsys_cs, γ)

# type1.
if typeof(As) == Vector{NMRSpectraSimulator.CompoundFIDType{Float64, NMRSpectraSimulator.SpinSysFIDType1{Float64}}}
    shift_constants = Δ_shifts
end


r = 3 # for glucose nrc.

if project_name == "test_serine1"
    r = 1 # for serine 700.
end
y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

U_rad_cost = U_cost .* (2*π)

###
LS_inds = 1:length(U_cost)

println("Timing: runalignment(), r = $(r)")
@time p_star, q, κ_BLS, getshiftfunc, getβfunc, getλfunc,
obj_func, N_vars_set, p_initial = NMRCalibrate.runalignment(shift_constants,
U_rad_cost, y_cost, LS_inds, Es, As, fs, SW;
max_iters = max_iters,
xtol_rel = 1e-5,
ftol_rel = 1e-5,
w = w,
κ_lb_default = κ_lb_default,
κ_ub_default = κ_ub_default,
λ_each_lb = 0.9,
λ_each_ub = 1.2)

p_star
κ_BLS
q

d_star = getshiftfunc(p_star)
β_star = getβfunc(p_star)
λ_star = getλfunc(p_star)



###### visualize.
T = Float64

import PlotlyJS
import Plots
Plots.plotly()

projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/$(project_name)"
plots_save_folder = joinpath(projects_dir, "plots")
isdir(plots_save_folder) || mkdir(plots_save_folder)
project_title = "test" # TODO change later.

# prep initial.
initial_cost = obj_func(p_initial)
NMRCalibrate.parseκ!(Es, ones(T, length(κ_BLS)))
fill!(w, 1.0)
q_initial_U = q.(U_rad)
#println("norm(q_initial_U) = ", norm(q_initial_U))

final_cost = obj_func(p_star)
q_final_U = q.(U_rad)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). cost vs. data (y)")



PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_initial_U), label = "q initial")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_final_U), "--", label = "q final")
#PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). data (y) vs. fit")

plots_save_path = joinpath(plots_save_folder, "region_$(r)_real.html")
title_string = "$(project_title): region $(r), real"
savefigfitresult(plots_save_path, title_string,
 real.(q_final_U), P, P_cost, real.(y_cost);
initial_fit = real.(q_initial_U))

plots_save_path = joinpath(plots_save_folder, "region_$(r)_imaginary.html")
title_string = "$(project_title): region $(r), imaginary"
savefigfitresult(plots_save_path, title_string,
imag.(q_final_U), P, P_cost, imag.(y_cost);
initial_fit = imag.(q_initial_U))

plots_save_path = joinpath(plots_save_folder, "region_$(r)_modulus.html")
title_string = "$(project_title): region $(r), modulus"
savefigfitresult(plots_save_path, title_string,
abs.(q_final_U), P, P_cost, abs.(y_cost);
initial_fit = abs.(q_initial_U))



### save.
save_folder_path = joinpath(projects_dir, project_name)
save_path = joinpath(save_folder_path, "results_$(r).bson")
BSON.bson(save_path,
p_star = p_star,
κ_lb_default = κ_lb_default,
κ_ub_default = κ_ub_default,
κ_star = κ_BLS,
d_star = d_star,
β_star = β_star,
λ_star = λ_star,
w = w,
# proxy setup-related below.
Δc_partition_radius = Δc_partition_radius,
tol_coherence = tol_coherence,
α_relative_threshold = α_relative_threshold,
λ0 = λ0,
Δcs_max = Δcs_max,
κ_λ_lb = κ_λ_lb,
κ_λ_ub = κ_λ_ub,
# experiement/cost-related below.
cost_inds_set = cost_inds_set,
Δsys_cs = Δsys_cs,
y_cost = y_cost,
U_cost = U_cost,
P_cost = P_cost,
As = As,
Es = Es,
y = y,
U_y = U_y,
fs = fs,
SW = SW,
ν_0ppm = ν_0ppm)
