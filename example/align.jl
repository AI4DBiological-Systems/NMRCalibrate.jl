# packaged version of nested_optim_d_only.jl

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

import OffsetArrays
import Interpolations

import Statistics
import Random

import NLopt
import MultistartOptimization

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# for convinence.
T = Float64
combinevectors = NMRSpectraSimulator.combinevectors

### user inputs.
max_iters = 50000
#max_iters = 5

fit_config_path = "/home/roy/Documents/repo/NMRData/input/fit_configs/calibrate_700MHz_type1_select_compounds.json"


project_name = "NRC-glucose-2018"

project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
project_folder = joinpath(project_base_folder, project_name)
#mkpath(project_folder)
### end inputs.


### load block.
load_path = joinpath(project_folder, "proxy.bson")
dict = BSON.load(load_path)


u_min = dict[:u_min]
u_max = dict[:u_max]
y = convert(Vector{Complex{Float64}}, dict[:y])
U_y = convert(Vector{Float64}, dict[:U_y])
SW = dict[:SW]
fs = dict[:fs]
ν_0ppm  = dict[:ν_0ppm]
λ0 = dict[:λ0]
molecule_names = dict[:molecule_names]
#w = dict[:w]
w = ones(length(dict[:As]))

As = collect( dict[:As][i] for i = 1:length(dict[:As]) )
Bs = collect( dict[:Bs][i] for i = 1:length(dict[:Bs]) )

### end load.

### prepare positions.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
P_y = hz2ppmfunc.(U_y)

ΩS0 = NMRSpectraSimulator.getΩS(As)
ΩS0_ppm = NMRSpectraSimulator.getPs(ΩS0, hz2ppmfunc)

Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds, cost_inds_set,
    λ_lbs, λ_ubs, κs_β_orderings,
    κs_β_DOFs = NMRCalibrate.prepareoptim(fit_config_path, molecule_names,
    hz2ppmfunc, U_y, y, As; region_min_dist = 0.1)


# visualize cost regions.
# sort(P_cost_set[1]) # etc..
P_cost_set = collect( P_y[cost_inds_set[r]] for r = 1:length(cost_inds_set) )


### optim all regions. # 900 secs.
region_select = 3
y_cost = y[cost_inds_set[region_select]]
U_cost = U_y[cost_inds_set[region_select]]
P_cost = P_y[cost_inds_set[region_select]]

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
PyPlot.title("positions against data spectrum, r = $(region_select), real part")



### optim.
N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
N_β = sum( NMRCalibrate.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
LS_inds = 1:length(U_cost)

# inputs.
κ_lb_default = 1e-2
κ_ub_default = 1e2
p_initial = zeros(N_d)
p_β = zeros(N_β)

# prepare.
U_rad_cost = U_cost .* (2*π)
Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )

a_setp, b_setp, minxs,
    rets = NMRCalibrate.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#
q, updatedfunc, getshiftfunc, N_vars_set,
run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc,
q_β = NMRCalibrate.setupcostnesteddwarp(Es, Bs, As, fs, SW, LS_inds, U_rad_cost,
    y_cost, Δsys_cs, a_setp, b_setp, κs_β_DOFs, κs_β_orderings;
    w = w,
    optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = κ_lb_default,
    κ_ub_default = κ_ub_default,
    max_iters = 500,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf)
#
obj_func = pp->NMRCalibrate.costnestedd(U_rad_cost, y_cost, updatedfunc, pp,
Es, Bs, κs_β_orderings, κs_β_DOFs, q, run_optim, E_BLS, κ_BLS, b_BLS, p_β)

grad_func = xx->FiniteDiff.finite_difference_gradient(f, xx)

#NMRCalibrate.evaldesignmatrixκ!(E_BLS, U_rad_cost, Es, As, w)

#optim_algorithm = :LN_BOBYQA
optim_algorithm = :GN_ESCH
#optim_algorithm = :GN_ISRES
#optim_algorithm = :GN_DIRECT_L # larger objective!
opt = NLopt.Opt(optim_algorithm, length(p_initial))


println("obj_func(p_initial) = ", obj_func(p_initial))
println()


# for display.
U = LinRange(u_min, u_max, 50000)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

q_U = q.(U_rad)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_U), label = "q")
#PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(region_select). initial vs data")

# ## test updates.
# N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
# N_β = sum( NMRCalibrate.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
#
# #p = ones(N_d)
# p = -ones(N_d) .* 0.8
#
# st_ind = 1
# NMRCalibrate.updatemixtured!(Bs, p, st_ind, fs, SW, Δsys_cs)
#
# p2 = similar(p)
# NMRCalibrate.applywarptoshifts!(p2, Bs, p, st_ind, a_setp, b_setp)
#
# p = ones(N_β)
# p[end-1] = 0.23
# p[end] = 1.23
# st_ind = 1
# NMRCalibrate.updateβ!(Bs, κs_β_orderings, κs_β_DOFs, p, st_ind)


println("Starting multi-start optim")
N_starts = 100
optim_algorithm = NLopt.LN_BOBYQA
xtol_rel = 1e-3
maxeval = 100
maxtime = Inf
# 900 seconds. obj_func(minx) = 0.01690795832619695

shift_lb = -ones(T, N_d)
shift_ub = ones(T, N_d)
shift_initial = zeros(T, N_d)

p_lb = shift_lb
p_ub = shift_ub
p_initial = shift_initial

prob = MultistartOptimization.MinimizationProblem(obj_func, p_lb, p_ub)

local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = optim_algorithm,
xtol_rel = xtol_rel,
maxeval = maxeval,
maxtime = maxtime)

multistart_method = MultistartOptimization.TikTak(N_starts)
@time ret_mo = MultistartOptimization.multistart_minimization(multistart_method,
    local_method, prob)
#
minf = ret_mo.value
minx = ret_mo.location
ret = ret_mo.ret
# 1249 secs.
#
#minx = [-0.09; -0.1; 1.1; 1.1] # a good answer.

minf2 = obj_func(minx)
println("obj_func(minx) = ", minf2)
println()

q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)


q_U = q2.(U_rad)
# q0_U = q.(U_rad)
# norm(q_U-q0_U)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_U), label = "q")
#PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(region_select). nested optim result vs data")

@assert 1==2

### cost func
r = 3

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)

y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

U_rad_cost = U_cost .* (2*π)

Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )


q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
E_BLS, κ_BLS, b_BLS, getshiftfunc, getβfunc, getλfunc,
N_vars_set = NMRCalibrate.setupcostalign(Es, Bs, As, κs_β_orderings,
κs_β_DOFs, fs, SW, LS_inds, U_rad_cost, y_cost, Δ_shifts; w = w)
