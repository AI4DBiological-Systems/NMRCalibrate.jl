
# packaged version of align.jl. Run align_prep.jl or prep_full.jl first.

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

import NLopt
import PlotlyJS
using Plots; plotly()

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# based on one region from pkg_calibrate.jl

r = 1
y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

local_optim_algorithm = NLopt.LN_BOBYQA
xtol_rel = 1e-9
maxeval = 50
maxtime = Inf
β_optim_algorithm = :GN_DIRECT_L
κ_lb_default = 1e-2
κ_ub_default = 1e2
β_max_iters = 500
β_xtol_rel = 1e-9
β_ftol_rel = 1e-9
β_maxtime = Inf

N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )

U_rad_cost = U_cost .* (2*π)

q, updatedfunc, getshiftfunc, N_vars_set,
run_optim, obj_func_β, E_BLS, κ_BLS, b_BLS, updateβfunc, updateκfunc,
q_β = NMRCalibrate.setupcostnesteddwarp(Es, Bs, As, fs, SW, LS_inds, U_rad_cost,
    y_cost, Δsys_cs, a_setp, b_setp; #κs_β_DOFs, κs_β_orderings;
    w = w,
    optim_algorithm = β_optim_algorithm,
    κ_lb_default = κ_lb_default,
    κ_ub_default = κ_ub_default,
    max_iters = β_max_iters,
    xtol_rel = β_xtol_rel,
    ftol_rel = β_ftol_rel,
    maxtime = β_maxtime)

#
#N_β = sum( NMRCalibrate.getNβ(Bs[n]) for n = 1:length(Bs) )
p_β = zeros(T, N_β) # persistant buffer.

obj_func = pp->NMRCalibrate.costnestedd(U_rad_cost, y_cost, updatedfunc, updateκfunc, pp,
Es, Bs, run_optim, E_BLS, κ_BLS, b_BLS, p_β)

### test.
p_test = rand(N_d)
obj_func(p_test)

#
# optim.
N_starts = 100
prob = MultistartOptimization.MinimizationProblem(obj_func, shift_lb, shift_ub)

local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = local_optim_algorithm,
xtol_rel = xtol_rel,
maxeval = maxeval,
maxtime = maxtime)

multistart_method = MultistartOptimization.TikTak(N_starts)
ret_mo = MultistartOptimization.multistart_minimization(multistart_method,
    local_method, prob)
