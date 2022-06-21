

using LinearAlgebra, FFTW
import BSON, Statistics, PyPlot, Random

import NMRSpectraSimulator

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

import MultistartOptimization, NLopt

#import Clustering

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# run prep_full.jl first.

r = 1
T = Float64
N_starts = 500
local_optim_algorithm = NLopt.LN_BOBYQA
xtol_rel = 1e-9
maxeval = 200 # 2 # 50
maxtime = Inf
β_optim_algorithm = :GN_DIRECT_L
w_lb_default = 1e-1
w_ub_default = 100.0
β_max_iters = 500 # 2 # 500
β_xtol_rel = 1e-9
β_ftol_rel = 1e-9
β_maxtime = Inf


N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
N_β = sum( NMRCalibrate.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )

shift_lb = -ones(N_d)
shift_ub = ones(N_d)

Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )

a_setp, b_setp, minxs,
    rets = NMRCalibrate.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#

y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
LS_inds = 1:length(U_cost)

# prepare.
#N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
@assert length(shift_ub) == length(shift_lb) == N_d

U_rad_cost = U_cost .* (2*π)

# setup inner optim over β.
q0, updatedfunc, getshiftfunc, N_vars_set,
run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc,
q_β = NMRCalibrate.setupcostnesteddwarpw(Es, Bs, As, fs, SW, LS_inds, U_rad_cost,
    y_cost, Δsys_cs, a_setp, b_setp, κs_β_DOFs, κs_β_orderings;
    β_optim_algorithm = β_optim_algorithm,
    w_lb_default = w_lb_default,
    w_ub_default = w_ub_default,
    β_max_iters = β_max_iters,
    β_xtol_rel = β_xtol_rel,
    β_ftol_rel = β_ftol_rel,
    β_maxtime = β_maxtime)

# set up outer optim over shifts.
N_β = sum( NMRCalibrate.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )
#N_β = sum( getNβ(Bs[n]) for n = 1:length(Bs) )
p_β = zeros(T, N_β) # persistant buffer.

obj_func = pp->NMRCalibrate.costnesteddw(U_rad_cost, y_cost, updatedfunc,
updatewfunc, pp, Es, Bs, κs_β_orderings, κs_β_DOFs, run_optim,
E_BLS, w_BLS, b_BLS, p_β)

#p = zeros(N_d)
p = [0.02;] # Serine 700 MHz BMRB.
#p = [0.13513301735885624;]
updatedfunc(p)

fill!(p_β, zero(T)) # always start from 0-phase?
minf, minx, ret, N_evals = run_optim(p_β)
p_β[:] = minx # take out?

NMRCalibrate.updateβ!(Bs, κs_β_orderings, κs_β_DOFs, p_β, 1)
updatewfunc(1.0)

cost = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2

#w_test = copy(ws)
w_test = copy(w_BLS) #[5.5;]

#### visualize.
N_viz = 50000

U = LinRange(u_min, u_max, N_viz)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w_test )
q_U = q2.(U_rad)

plotregion(P, U, q_U, P_y, y, P_cost, y_cost, 0.0, 1,
    "", "test", "";
    save_plot_flag = false,
    display_plot_flag = true,
    canvas_size = (1000,400))

# revisit cost again.

function evalcost(q2, U_cost_rad, y_cost)
    cost2 = 0.0
    for m = 1:length(U_cost_rad)
        cost2 += abs2(q2(U_cost_rad[m]) - y_cost[m])
    end

    return cost2
end

U_cost_rad = U_cost .* (2*π)
cost2 = evalcost(q2, U_cost_rad, y_cost)

#norm(b_BLS - reinterpret(T, y_cost))
@assert 12==3

q3, updateβfunc3, updatewfunc3, E_BLS3, w_BLS3, b_BLS3,
getβfunc3 = NMRCalibrate.setupcostβLSw(Es, Bs, As, LS_inds, U_rad_cost, y_cost,
    κs_β_DOFs, κs_β_orderings;
    w_lb_default = w_lb_default,
    w_ub_default = w_ub_default)
#
updatewfunc3(1.0)
E_BLS3*w_BLS3 - E_BLS*w_BLS
norm(E_BLS3- E_BLS)

(E_BLS3*w_BLS3)[1]
u_rad = U_cost[1]*2*pi
println("q0(u_rad) = ", q0(u_rad))
println("q3(u_rad) = ", q3(u_rad))
