
# TODO: move plot resonance to NMRSpectraSimulator

import NMRDataSetup
import NMRSpectraSimulator
import NMRSpecifyRegions

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JLD

#import Clustering
import Statistics
import Random

include("./helpers/final_helpers.jl")
include("./helpers/solute_calibrate.jl")
include("./helpers/loop_entries1.jl")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


max_iters = 50000
#max_iters = 5


# dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType1(0.0)
# project_name = "test_serine1"
# molecule_names = ["L-Serine";]
# w = [1.0; ]

dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType1(0.0)
project_name = "test_glucose1"
molecule_names = ["D-(+)-Glucose";]
w = [1.0; ]

# dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType1(0.0)
# project_name = "D-(+)-Glucose-700-r3"
# molecule_names = ["D-(+)-Glucose";]
# w = [1.0; ]

# dummy_SSFID = NMRSpectraSimulator.SpinSysFIDType2(0.0)
# project_name = "test_glucose2"
# molecule_names = ["D-(+)-Glucose";]
# w = [1.0; ]

# project_name = "D-(+)-Glucose-NRC-600"
# molecule_names = ["D-(+)-Glucose";]
# w = [1.0; ]
#
# project_name = "L-Serine-700"
# molecule_names = ["L-Serine"; "DSS"]
# w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
#
# project_name = "D-(+)-Glucose-700"
# molecule_names = ["D-(+)-Glucose"; "DSS"]
# w = [20.0/4.6; 1.0] # BMRB: DSS is 0.1 % => 4.6 mM
#
# project_name = "L-Isoleucine-700"
# molecule_names = ["L-Isoleucine"; "DSS"]
# w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
#
# project_name = "L-Leucine-500-2mM"
# molecule_names = ["L-Leucine";]
# w = [1.0]
#
# project_name = "L-Glutamine-700"
# molecule_names = ["L-Glutamine"; "DSS"]
# w = [20.0/0.5; 1.0] # BMRB: DSS is 500 uM => 0.5 mM

projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final"
base_path_JLD = "/home/roy/Documents/repo/NMRData//src/input/molecules"
cs_config_path = "/home/roy/Documents/repo/NMRData/src/input/reduced_cs_config.txt"




println("Now on $(project_name)")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.
# 0.1% DSS is 0.0046 M = 4.6 mM.
save_folder_path = joinpath(projects_dir, project_name)
#@assert 1==2

w = w ./ norm(w) # since the fit data, y, is normalized.

# proxy-related.
tol_coherence = 1e-2
??_relative_threshold = 0.05
#??c_partition_radius = 1e-1 # 17 element partition for glucose spin group 2.
??c_partition_radius_candidates = [1e-1; 0.3; 0.5; 0.7; 0.8; 0.9]
#??c_partition_radius_candidates = [0.9;]
??0 = 3.4
??cs_max = 0.2 # for proxy.
??_??_lb = 0.5
??_??_ub = 2.5

# if a ??c_partition_radius candidate gives max partition sizes less than of
#   equal to this value, then use that candidate.
early_exit_part_size = 7

### end inputs.



## load data.
load_folder_path = joinpath(projects_dir, project_name)
load_path = joinpath(load_folder_path, "$(project_name).bson")
dic = BSON.load(load_path)
s_t = dic[:s_t]
fs = dic[:fs]
SW = dic[:SW]
??_0ppm = dic[:??_0ppm]


hz2ppmfunc = uu->(uu - ??_0ppm)*SW/fs
ppm2hzfunc = pp->(??_0ppm + pp*fs/SW)

offset_Hz = ??_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(N, fs, offset_Hz)

Z = maximum(abs.(DFT_s))
y = DFT_s[U_inds] ./ Z


#
S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(S_U), label = "data spectrum")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data spectra")




####### mixture proxy.
??cs_max_mixture = collect( ??cs_max for i = 1:length(molecule_names))

println("Timing: trydiff??cradius()")
@time mixture_params, ??c_partition_radius = NMRCalibrate.trydiff??cradius(??c_partition_radius_candidates,
    molecule_names, base_path_JLD, ??cs_max_mixture, hz2ppmfunc, ppm2hzfunc,
    fs, SW, ??0, ??_0ppm, early_exit_part_size, ??cs_max, tol_coherence,
    ??_relative_threshold, dummy_SSFID)
As = mixture_params




??S_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
??S_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(??S_ppm))

println("$(project_name): Partition sizes:")
display(NMRCalibrate.displaypartitionsizes(As[1]))
println("??c_partition_radius = ", ??c_partition_radius)
println()

#@assert 1==2

# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(3.0)

u_offset = 0.5
u_min = ppm2hzfunc(??S_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(??S_ppm_sorted[end] + u_offset)

println("Timing: fitproxies!()")
@time NMRSpectraSimulator.fitproxies!(As;
??_??_lb = ??_??_lb,
??_??_ub = ??_??_ub,
u_min = u_min,
u_max = u_max,
??r = 1.0,
????_?? = 0.05)

######## end mixture proxy.

N_d = sum( NMRCalibrate.getNd(As[n]) for n = 1:length(As) )
#N_?? = sum( NMRCalibrate.getN??(As[n]) for n = 1:length(As) )
N_?? = sum( NMRCalibrate.getN??(As[n]) for n = 1:length(As) )

####### cost func.
combinevectors = NMRSpectraSimulator.combinevectors

??sys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds,
cost_inds_set = NMRCalibrate.prepareoptim(cs_config_path, molecule_names, hz2ppmfunc,
U_y, y, As; region_min_dist = 0.1)


# visualize cost regions.
# sort(P_cost_set[1]) # etc..
P_cost_set = collect( P_y[cost_inds_set[r]] for r = 1:length(cost_inds_set) )


### optim all regions. # 900 secs.
y_cost = y_cost_all
U_cost = U_cost_all
P_cost = P_cost_all

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "data spectrum")
PyPlot.plot(P_cost, real.(y_cost), "^", label = "positions")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("positions against data spectrum, real part")

#@assert 1==2

??_shifts = NMRSpectraSimulator.combinevectors(??sys_cs)

##### set up updates.

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*??)
#??S_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.??s) ./ (2*??) ) for A in mixture_params )



## parameters that affect qs.
# A.ss_params.d, A.ss_params.??s_??, A.ss_params.??s_??
# A.d_singlets, A.??s_singlets, A.??s_singlets, A.??_singlets, A.??0, A.??s_??_singlets
# purposely perturb ??.

Es = collect( NMRSpectraSimulator.??CompoundFIDType(As[i]) for i = 1:length(As) )


??_lb_default = 0.2
??_ub_default = 50.0

# TODO here, change shift delta cs to each group.


# default, type2.
?? = 0.001 # in ppm.
shift_constants = (??sys_cs, ??)

# type1.
if typeof(As) == Vector{NMRSpectraSimulator.CompoundFIDType{Float64, NMRSpectraSimulator.SpinSysFIDType1{Float64}}}
    shift_constants = ??_shifts
end


r = 3 # for glucose nrc.

if project_name == "test_serine1"
    r = 1 # for serine 700.
end
y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

U_rad_cost = U_cost .* (2*??)

###
LS_inds = 1:length(U_cost)

println("Timing: runalignment(), r = $(r)")
@time p_star, q, ??_BLS, getshiftfunc, get??func, get??func,
obj_func, N_vars_set, p_initial = NMRCalibrate.runalignment(shift_constants,
U_rad_cost, y_cost, LS_inds, Es, As, fs, SW;
max_iters = max_iters,
xtol_rel = 1e-5,
ftol_rel = 1e-5,
w = w,
??_lb_default = ??_lb_default,
??_ub_default = ??_ub_default,
??_each_lb = 0.9,
??_each_ub = 1.2)

p_star
??_BLS
q

d_star = getshiftfunc(p_star)
??_star = get??func(p_star)
??_star = get??func(p_star)



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
NMRCalibrate.parse??!(Es, ones(T, length(??_BLS)))
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
save_path = joinpath(save_folder_path, "results_$(r).bson")
BSON.bson(save_path,
p_star = p_star,
??_lb_default = ??_lb_default,
??_ub_default = ??_ub_default,
??_star = ??_BLS,
d_star = d_star,
??_star = ??_star,
??_star = ??_star,
w = w,
# proxy setup-related below.
??c_partition_radius = ??c_partition_radius,
tol_coherence = tol_coherence,
??_relative_threshold = ??_relative_threshold,
??0 = ??0,
??cs_max = ??cs_max,
??_??_lb = ??_??_lb,
??_??_ub = ??_??_ub,
# experiement/cost-related below.
cost_inds_set = cost_inds_set,
??sys_cs = ??sys_cs,
y_cost = y_cost,
U_cost = U_cost,
P_cost = P_cost,
As = As,
Es = Es,
y = y,
U_y = U_y,
fs = fs,
SW = SW,
??_0ppm = ??_0ppm)
