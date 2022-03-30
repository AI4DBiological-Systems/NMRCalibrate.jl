

using LinearAlgebra, FFTW
import BSON, Statistics, PyPlot, Random

import NMRSpectraSimulator

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

import MultistartOptimization, NLopt
import MonotoneMaps
import Interpolations

include("./helpers/warp_helpers.jl")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])



projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/test_glucose1/"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-700-r3/"
r = 3

T = Float64

### load block.
#load_path = joinpath(joinpath(projects_dir, project_name), "results_full.bson")
load_path = joinpath(projects_dir, "results_$(r).bson")
dict = BSON.load(load_path)

Es = collect( dict[:Es][i] for i = 1:length(dict[:Es]) )
As = collect( Es[n].core for n = 1:length(Es))

Δsys_cs = convert(Vector{Vector{Float64}}, dict[:Δsys_cs])
y = convert(Vector{Complex{Float64}}, dict[:y])
U_y = convert(Vector{Float64}, dict[:U_y])
SW = dict[:SW]
fs = dict[:fs]
ν_0ppm = dict[:ν_0ppm]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
P_y = hz2ppmfunc.(U_y)


ΩS_ppm = NMRCalibrate.findfreqrange(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))
u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)


P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

p_star = dict[:p_star]
κ_lb_default = dict[:κ_lb_default]
κ_ub_default = dict[:κ_ub_default]
κ_star = dict[:κ_star]
d_star = dict[:d_star]
β_star = dict[:β_star]
λ_star = dict[:λ_star]
cost_inds_set = dict[:cost_inds_set]
w = dict[:w]


### end load block.

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)

y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

U_rad_cost = U_cost .* (2*π)

LS_inds = 1:length(U_cost)

q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
E_BLS, κ_BLS, b_BLS, getshiftfunc, getβfunc, getλfunc,
N_vars_set = NMRCalibrate.setupcostcLshiftLS(Es, As, fs, SW,
LS_inds, U_rad_cost, y_cost, Δ_shifts;
w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)

obj_func = pp->NMRCalibrate.costcLshift(U_rad_cost, y_cost,
updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es,
E_BLS, κ_BLS, b_BLS, q)





#q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)

q_U = q.(U_rad)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_U), label = "q")
#PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). load result vs data")

#@assert 1==2



### nested optim.

N_d = sum( NMRCalibrate.getNd(As[n]) for n = 1:length(As) )
N_β = sum( NMRCalibrate.getNβ(As[n]) for n = 1:length(As) )
N_λ = sum( NMRCalibrate.getNλ(As[n]) for n = 1:length(As) )
p_β = zeros(N_β)

#### initial values.
shift_lb = -ones(T, N_d)
shift_ub = ones(T, N_d)
shift_initial = zeros(T, N_d)
λ_each_lb = 0.2
λ_each_ub = 50.0

λ_lb = λ_each_lb .* ones(T, N_λ)
λ_ub = λ_each_ub .* ones(T, N_λ)
λ_initial = ones(T, N_λ)

### set up constraints and initial guess.
p_lb = [ shift_lb; λ_lb ]
p_ub = [ shift_ub; λ_ub ]
#p_initial = [shift_initial; λ_initial]

# loaded answer.
p_initial = [ -0.11701290729148561;
 -0.9360412587205644;
 1.1992580872763856
 1.1984875571682057]

## want to get to this answer.
# p_initial = [ -0.11701290729148561;
# -0.1;
# 1.1992580872763856
# 1.1984875571682057]

##########
# st_ind_d = 1
# fin_ind_d = st_ind_d + N_d - 1
# updatedfunc = pp->NMRCalibrate.updatemixtured!(As, pp, st_ind_d, fs, SW, shift_constants)
#
# #λupdate.
# st_ind_λ = fin_ind_d + 1
# fin_ind_λ = st_ind_λ + N_λ -1
# updateλfunc = pp->NMRCalibrate.updateλ!(As, pp, st_ind_λ)

a_setp, b_setp,
    minxs, rets = setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)


q, updatedfunc, updateλfunc, getshiftfunc, getλfunc, N_vars_set,
run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc,
q_β = NMRCalibrate.setupcostnestedλdwarp(Es, As, fs, SW, LS_inds, U_rad_cost,
    y_cost, Δ_shifts, a_setp, b_setp;
    w = w,
    optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = κ_lb_default,
    κ_ub_default = κ_ub_default,
    max_iters = 500,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf)

obj_func = pp->NMRCalibrate.costnestedλd(U_rad_cost, y_cost, updatedfunc, updateλfunc, pp,
Es, As, q, run_optim, E_BLS, κ_BLS, b_BLS, p_β)

grad_func = xx->FiniteDiff.finite_difference_gradient(f, xx)

#optim_algorithm = :LN_BOBYQA
optim_algorithm = :GN_ESCH
#optim_algorithm = :GN_ISRES
#optim_algorithm = :GN_DIRECT_L # larger objective!
opt = NLopt.Opt(optim_algorithm, length(p_initial))

max_iters = 1000
xtol_rel = 1e-3
ftol_rel = 1e-3
maxtime = Inf

println("obj_func(p_initial) = ", obj_func(p_initial))
println()

q_U = q.(U_rad)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_U), label = "q")
#PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). initial vs data")

# println("Starting optim.")
# @time minf, minx, ret, N_evals = NMRCalibrate.runNLopt!(opt,
#     p_initial,
#     obj_func,
#     grad_func,
#     p_lb,
#     p_ub;
#     max_iters = max_iters,
#     xtol_rel = xtol_rel,
#     ftol_rel = ftol_rel,
#     maxtime = maxtime)

#
println("Starting multi-start optim")
N_starts = 100
optim_algorithm = NLopt.LN_BOBYQA
xtol_rel = 1e-3
maxeval = 100
maxtime = Inf

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

q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)


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
PyPlot.title("r = $(r). nested optim result vs data")
