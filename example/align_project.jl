

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

import NLopt
import Statistics
import Random

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


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

# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
# experiment_path = joinpath(root_folder, "NRC-glucose-2018") # second input is subfolder name.
#
# molecule_names = ["D-(+)-Glucose";]
# w = [1.0;]
#
# project_name = "NRC-glucose-2018"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
# ###

# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-700-20mM"
# experiment_path = joinpath(root_folder, "L-Serine") # second input is subfolder name.
#
# molecule_names = ["L-Serine";]
# w = [1.0;]
# save_flag = true
#
# project_name = "L-Serine"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/BMRB-700-20mM"
# ###

# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
# experiment_path = joinpath(root_folder, "L-Serine") # second input is subfolder name.
#
# molecule_names = ["L-Serine";]
# w = [1.0;]
# save_flag = true
#
# project_name = "L-Serine"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/BMRB-500-0.5mM"
# ###

# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-2mM"
# experiment_path = joinpath(root_folder, "L-Serine") # second input is subfolder name.
#
# molecule_names = ["L-Serine";]
# w = [1.0;]
# save_flag = true
#
# project_name = "L-Serine"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/BMRB-500-2mM"
# ###
#
# load_path = joinpath(experiment_path, "experiment.bson")
#
# project_folder = joinpath(project_base_folder, project_name)
# save_path = joinpath(project_folder, "alignment_results.bson")
# mkpath(project_folder)
# #####

# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
# experiment_path = joinpath(root_folder, "NRC-4_amino_acid-Jan2022-1") # second input is subfolder name.
#
# molecule_names = ["L-Serine";]
# w = [1.0;]
# save_flag = true
#
# project_name = "Jan2022-Serine"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
# ###
#
# load_path = joinpath(experiment_path, "experiment.bson")
#
# project_folder = joinpath(project_base_folder, project_name)
# save_path = joinpath(project_folder, "alignment_results.bson")
# mkpath(project_folder)
# #####


# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
# experiment_path = joinpath(root_folder, "NRC-4_amino_acid-Jan2022-1") # second input is subfolder name.
#
# molecule_names = ["L-Isoleucine";]
# w = [1.0;]
# save_flag = true
#
# project_name = "Jan2022-Isoleucine"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
# ###
#
# load_path = joinpath(experiment_path, "experiment.bson")
#
# project_folder = joinpath(project_base_folder, project_name)
# save_path = joinpath(project_folder, "alignment_results.bson")
# mkpath(project_folder)
# #####

# # where the experiment bson file is located.
# root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-0.5mM"
# experiment_path = joinpath(root_folder, "L-Isoleucine") # second input is subfolder name.
#
# molecule_names = ["L-Isoleucine";]
# w = [1.0;]
# save_flag = true
#
# project_name = "L-Isoleucine"
# project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/BMRB-500-0.5mM"
# ###
#
# load_path = joinpath(experiment_path, "experiment.bson")
#
# project_folder = joinpath(project_base_folder, project_name)
# save_path = joinpath(project_folder, "alignment_results.bson")
# mkpath(project_folder)
# #####

# where the experiment bson file is located.
root_folder = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
experiment_path = joinpath(root_folder, "NRC-4_amino_acid-Jan2022-1") # second input is subfolder name.

molecule_names = ["L-Leucine";]
w = [1.0;]
save_flag = true

project_name = "Jan2022-Leucine"
project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
###

load_path = joinpath(experiment_path, "experiment.bson")

project_folder = joinpath(project_base_folder, project_name)
save_path = joinpath(project_folder, "alignment_results.bson")
mkpath(project_folder)
#####

##### user input for loading data.

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


@time obj_funcs, minfs, minxs, rets, As, Bs, Es, cost_inds_set, u_min,
u_max = NMRCalibrate.alignproject(save_path,
    SH_config_path,
    surrogate_config_path,
    fit_config_path,
    H_params_path,
    dict_compound_to_filename,
    molecule_names,
    w,
    dummy_SSFID,
    y,
    U_y,
    SW,
    fs,
    ??_0ppm,
    ??_0ppm;
    u_offset = 0.5,
    unique_cs_atol = 1e-6,
    tol_coherence = 1e-2,
    ??_relative_threshold = 0.05,
    ??c_partition_radius = 0.3,
    ??r_default = 1.0,
    ????_??_default = 0.05,
    ??cs_max_scalar_default = 0.2,
    ??_??_lb_default = 0.7,
    ??_??_ub_default = 1.5,
    offset_ppm = 0.3,
    region_min_dist = 0.1,
    sigmoid_lb = 0.1,
    sigmoid_ub = 0.7,
    N_sigmoid_samples = 10,
    sigmoid_optim_algorithm = :LN_BOBYOA,
    save_flag = save_flag,
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 500,
    maxtime = Inf,
    ??_max_iters = 500,
    ??_xtol_rel = 1e-9,
    ??_ftol_rel = 1e-9,
    ??_maxtime = Inf);

####
# if load_dict_flag
#     load_path = joinpath(project_folder, "alignment_results.bson")
#     dict = BSON.bson(load_path)
#     region_min_dist = dict[:region_min_dist]
#     minfs = dict[:minfs]
#     minxs = dict[:minxs]
#     rets = dict[:rets]
# end

#### visualize.



function visualizereresults(U, P, As, Es, w, y, U_y, P_y, cost_inds_set, minxs, fig_num)

    @assert length(cost_inds_set) == length(minxs)

    U_rad = U .* (2*??)

    for r = 1:length(minxs)

        y_cost = y[cost_inds_set[r]]
        U_cost = U_y[cost_inds_set[r]]
        P_cost = P_y[cost_inds_set[r]]



        q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)

        obj_funcs[r](minxs[r])

        #NMRCalibrate.reset??!(Es)
        q_U = q2.(U_rad)

        PyPlot.figure(fig_num)
        fig_num += 1

        PyPlot.plot(P_y, real.(y), label = "y")
        PyPlot.plot(P, real.(q_U), label = "q")
        PyPlot.plot(P_cost, real.(y_cost), "x")

        PyPlot.legend()
        PyPlot.xlabel("ppm")
        PyPlot.ylabel("real")
        PyPlot.title("Alignment result: region $(r)")

    end

    return fig_num
end

U = LinRange(u_min, u_max, 50000)
P = hz2ppmfunc.(U)
fig_num = visualizereresults(U, P, As, Es, w, y, U_y, P_y, cost_inds_set, minxs, fig_num)
