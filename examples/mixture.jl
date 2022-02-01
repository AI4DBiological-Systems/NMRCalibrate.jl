

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

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"

# project_name = "glucose-700"
# molecule_names = ["D-(+)-Glucose"; "DSS"]
# w = [20.0/4.6; 1.0] # BMRB-700 glucose: DSS is 0.0046 M = 4.6 mM.

project_name = "Nam-Jan2022"

# molecule_names = ["L-Serine"; "L-Glutamine"]
# w = [1.0; 1.0] # BMRB-700 phenylalanine: DSS is 500 micro M.
# #w = w ./ norm(w)

molecule_names = ["L-Serine"; "D-(+)-Glucose"]
w = [1.0; 1.0] # BMRB-700 phenylalanine: DSS is 500 micro M.


# path to the GISSMO Julia storage folder.
base_path_JLD = "/home/roy/Documents/repo/NMRData//src/input/molecules"

# proxy-related.
tol_coherence = 1e-2
α_relative_threshold = 0.05
λ0 = 3.4
Δcs_max = 0.2
κ_λ_lb = 0.5
κ_λ_ub = 2.5

### end inputs.



## load data.
load_folder_path = joinpath(projects_dir, project_name)
load_path = joinpath(load_folder_path, "$(project_name).bson")
dic = BSON.load(load_path)
s_t = dic[:s_t]
fs = dic[:fs]
SW = dic[:SW]
ν_0ppm = dic[:ν_0ppm]

# normalize.
s_t = s_t

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRCalibrate.getwraparoundDFTfreqs(N, fs, offset_Hz)

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


Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold)
As = mixture_params

ΩS_ppm = NMRCalibrate.findfreqrange(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(3.0)

u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

NMRSpectraSimulator.fitproxies!(As;
#NMRSpectraSimulator.fitproxiessimple!(As;
    κ_λ_lb = κ_λ_lb,
    κ_λ_ub = κ_λ_ub,
    u_min = u_min,
    u_max = u_max,
    Δr = 1.0,
    Δκ_λ = 0.05)

### cost func.
combinevectors = NMRSpectraSimulator.combinevectors

cs_config_path = "/home/roy/MEGAsync/inputs/NMR/configs/reduced_cs_config.txt"

Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds,
cost_inds_set = NMRCalibrate.prepareoptim(cs_config_path, molecule_names, hz2ppmfunc,
U_y, y, As; region_min_dist = 0.1)


# visualize cost regions.
# sort(P_cost_set[1]) # etc..
P_cost_set = collect( P_y[cost_inds_set[r]] for r = 1:length(cost_inds_set) )


# ### just optim region r.
r = 3
y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

# ### just optim a few regions.
# R = [1; 2; 3]
# m_inds = mergeinds(cost_inds_set, R)
# y_cost = y[m_inds...]
# U_cost = U_y[m_inds...]
# P_cost = P_y[m_inds...]

# ### optim all regions. # 900 secs.
# y_cost = y_cost_all
# U_cost = U_cost_all
# P_cost = P_cost_all

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "data spectrum")
PyPlot.plot(P_cost, real.(y_cost), "^", label = "positions")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("positions against data spectrum, real part")

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)

##### set up updates.

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
#ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )


# lorentzian (oracle/ground truth)
f = uu->NMRSpectraSimulator.evalmixture(uu, mixture_params)

U_LS = U

## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
# purposely perturb κ.

As2 = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:length(As) )

# assemble proxy and compare against oracle at locations in U.
#w = ones(length(As))
#q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As2; w = w)
Es = As2
#w = ones(length(Es))
LS_inds = 1:length(U_cost)
max_iters = 5000
#max_iters = 20000

q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
κ_BLS, getshiftfunc, getβfunc, getλfunc,
N_vars_set = NMRCalibrate.setupcostcLshiftLS(Es, As, fs, SW,
    LS_inds, U_cost, y_cost, Δ_shifts;
    w = w, κ_lb_default = 0.2, κ_ub_default = 50.0)

N_d, N_β, N_λ = N_vars_set
p_initial = [zeros(N_d); zeros(N_β); ones(N_λ)]

#updateκfunc(p_initial)
fill!(κ_BLS, 1.0)

### plot.

#f_U = f.(U)
q_U = q.(U)
### for now.
# discrepancy = abs.(f_U-q_U)
# max_val, ind = findmax(discrepancy)
# println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
# println("max discrepancy: ", max_val)
# println()
#
# PyPlot.figure(fig_num)
# fig_num += 1
#
# PyPlot.plot(P, real.(f_U), label = "f")
# PyPlot.plot(P, real.(q_U), "--", label = "q")
#
# PyPlot.legend()
# PyPlot.xlabel("ppm")
# PyPlot.ylabel("real")
# PyPlot.title("f vs q")



# @assert 1==2
println("Timing: runalignment()")
@time p_star, q, κ_BLS, getshiftfunc, getβfunc, getλfunc,
obj_func, N_vars_set = NMRCalibrate.runalignment(Δ_shifts,
U_cost, y_cost, LS_inds, Es, As, fs, SW;
max_iters = 50000,
xtol_rel = 1e-7,
ftol_rel = 1e-12,
w = w,
κ_lb_default = 0.2,
κ_ub_default = 50.0,
λ_each_lb = 0.9,
λ_each_ub = 1.1)

d_star = getshiftfunc(p_star)
β_star = getβfunc(p_star)
λ_star = getλfunc(p_star)
println("d_star = ", d_star)
println("β_star = ", β_star)
println("λ_star = ", λ_star)
println("κ_BLS = ", κ_BLS)
println()

### reference.

# reference, zero shift, phase.
p_test = zeros(sum(N_vars_set))
initial_cost = obj_func(p_test)
q_initial_U = q.(U)

p_test = copy(p_star)
final_cost = obj_func(p_test)
q_final_U = q.(U)


PyPlot.figure(fig_num)
fig_num += 1


PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("cost vs. data (y)")


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_initial_U), label = "q initial")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_final_U), "--", label = "q final")
#PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data (y) vs. fit")



####################


#include("del.jl") # explore here.
