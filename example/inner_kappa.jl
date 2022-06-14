

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

r = 2
T = Float64


y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
N_β = sum( NMRCalibrate.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )

shift_lb = -ones(N_d)
shift_ub = ones(N_d)

Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )

a_setp, b_setp, minxs,
    rets = NMRCalibrate.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#


d_test = zeros(N_d)

U_rad_cost = U_cost .* (2*π)

LS_inds = 1:length(U_cost)
optim_algorithm = :GN_DIRECT_L
β0 = zeros(N_β)

# # κ
# run_optim, obj_func3, E_BLS3, κ_BLS3, b_BLS3, updateβfunc3,
# q3 = NMRCalibrate.setupβLSsolver(optim_algorithm,
#     Es, Bs, As, LS_inds, U_rad_cost, y_cost;
#     κ_lb_default = 0.2,
#     κ_ub_default = 5.0,
#     max_iters = 500,
#     xtol_rel = 1e-9,
#     ftol_rel = 1e-9,
#     maxtime = Inf);
# @time run_optim(β0)
# display(κ_BLS3)

# # w
# run_optimw, obj_func, E_BLS, w_BLS, b_BLS, updateβfunc,
# q = NMRCalibrate.setupβLSsolverw(optim_algorithm,
#     Es, Bs, As, LS_inds, U_rad_cost, y_cost;
#     w_lb_default = 0.2,
#     w_ub_default = 5.0,
#     max_iters = 500,
#     xtol_rel = 1e-9,
#     ftol_rel = 1e-9,
#     maxtime = Inf);
# @time run_optimw(β0)
# display(w_BLS)


### debug.
w_lb_default = 0.2
w_ub_default = 3.0

N_β = sum( NMRCalibrate.getNβ(Bs[n]) for n = 1:length(Bs) )

β_lb = ones(T, N_β) .* (-π)
β_ub = ones(T, N_β) .* (π)

p_lb = β_lb
p_ub = β_ub

q, updateβfunc, updatewfunc, E_BLS, w_BLS, b_BLS,
getβfunc = NMRCalibrate.setupcostβLSw(Es, Bs, As, LS_inds, U_rad_cost, y_cost;
    w_lb_default = w_lb_default,
    w_ub_default = w_ub_default)

display(w)
updatewfunc(1.0)
display(w)

# I am here. improve timing.
@assert 66==5

Working on region 1
ret_mo = (value = 0.5192247403187571, location = [1.0, 1.0, 1.0], ret = :XTOL_REACHED)
Working on region 2
ret_mo = (value = 0.12538890087070964, location = [0.9945040434629705, 0.06125239652484532, 0.09882763402500577], ret = :MAXEVAL_REACHED)
Working on region 3
ret_mo = (value = 1.0022368527753995, location = [0.04556125088664971, 0.08617822515030829, 0.12157829228679165], ret = :MAXEVAL_REACHED)
3302.990934 seconds (53.47 G allocations: 1.179 TiB, 4.42% gc time, 0.24% compilation time)


# q is simulated spectra. f is cost function.
f = pp->NMRCalibrate.costβLSw(U_rad_cost, y_cost, updateβfunc, updatewfunc, pp, Es,
E_BLS, w_BLS, b_BLS, q)

df = xx->FiniteDiff.finite_difference_gradient(f, xx)

opt = NLopt.Opt(optim_algorithm, N_β)

run_optim = pp->NMRCalibrate.runNLopt!(opt,
    pp,
    f,
    df,
    p_lb,
    p_ub;
    max_iters = max_iters,
    xtol_rel = xtol_rel,
    ftol_rel = ftol_rel,
    maxtime = maxtime)

@assert 1==2

using BenchmarkTools

### BLS timings.
println("Timing: evaldesignmatrixκ!")
Q = zeros(Complex{Float64}, length(U_cost), N_κ + N_κ_singlets)
@btime NMRCalibrate.evaldesignmatrixκ!(Q, U_rad_cost, Es, w);

#
# println("q.(U_rad_cost)")
# @btime q.(U_rad_cost);
#
# println("obj_func.(p_star)")
# @btime obj_func(p_star);
#
# obj_func(p_star)
#
# q_U = q.(U_rad_cost)
#
# b_BLS = b = collect(reinterpret(Float64, y_cost))
#
# @btime B = reinterpret(Float64, Q);
# B = reinterpret(Float64, Q)
# @btime c = B*κ_BLS - b_BLS;
# c = B*κ_BLS - b_BLS
#
# norm(q.(U_rad_cost) - y_cost)^2 - obj_func(p_star)
#
# L = 50
# a = randn(Complex{Float64}, L)
# b = randn(Complex{Float64}, L)
#
# sum( abs2(a[i]-b[i]) for i in eachindex(a))
#
# sum( real(a[i]-b[i])^2 + imag(a[i]-b[i])^2 for i in eachindex(a) )

###########################


### reference.

# reference, zero shift, phase.
N_d = sum( NMRCalibrate.getNd(As[n]) for n = 1:length(As) )
N_β = sum( NMRCalibrate.getNβ(As[n]) for n = 1:length(As) )
N_λ = sum( NMRCalibrate.getNλ(As[n]) for n = 1:length(As) )
#shift_manual = zeros(T, N_d)
#β_manual = zeros(T, N_β)
#λ_manual = ones(T, N_λ)

#shift_manual = [-0.08; -0.1]
shift_manual = [-0.11; -0.11]
β_manual = copy(β_star)
λ_manual = copy(λ_star)

p_manual = [shift_manual; β_manual; λ_manual]

manual_cost = obj_func(p_manual)
#NMRCalibrate.parseκ!(Es, ones(T, length(κ_BLS)))
fill!(w, 1.0)
q_manual_U = q.(U_rad)
println("norm(q_manual_U) = ", norm(q_manual_U))


final_cost = obj_func(p_star)
q_star_U = q.(U_rad)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_manual_U), label = "q manual")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_star_U), "--", label = "q star")
#PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). data vs. manual")


### z.
Gs = collect( NMRSpectraSimulator.zCompoundFIDType(As[i]) for i = 1:length(As) )

g, updateβfunc5, updatezfunc5,
E_LS5, z_LS, b_LS5,
getβfunc5 = NMRCalibrate.setupcostβLS(Gs, As, LS_inds, U_rad_cost, y_cost)

println("Timing: updatezfunc() and parsez!()")
@time updatezfunc5(0.0)
@time NMRCalibrate.parsez!(Gs, z_LS)


g_star_U = g.(U_rad)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_manual_U), label = "q manual")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(g_star_U), "--", label = "g star")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). LS z vs data")

### fit β.

#@assert 1==2
#β_initial = ones(length(getβfunc(p_star)))
β_initial = copy(getβfunc(p_star))

# optim_algorithm = :LN_BOBYQA # good.
# β0 = copy(β_initial)

# optim_algorithm = :GN_ESCH
# β0 = copy(β_initial) # ok
#
# optim_algorithm = :GN_ESCH
# β0 = copy(β_initial)

# optim_algorithm = :GN_ISRES
# β0 = copy(β_initial) # not good.

optim_algorithm = :GN_DIRECT_L
β0 = copy(β_initial) # good.



# ## packaged up.
# run_optim, obj_func3, E_BLS3, κ_BLS3, b_BLS3, updateβfunc3,
# q3 = NMRCalibrate.setupβLSsolver(optim_algorithm,
#     Es, As, LS_inds, U_rad_cost, y_cost;
#     κ_lb_default = 0.2,
#     κ_ub_default = 5.0,
#     max_iters = 500,
#     xtol_rel = 1e-9,
#     ftol_rel = 1e-9,
#     maxtime = Inf);

optim_algorithm = NLopt.LN_BOBYQA
#optim_algorithm = NLopt.GN_ESCH
run_optim, obj_func3, E_BLS3, κ_BLS3, b_BLS3, updateβfunc3,
q3 = NMRCalibrate.setupβLSsolverMultistartoptim(optim_algorithm,
    Es, As, LS_inds, U_rad_cost, y_cost;
    κ_lb_default = 0.2,
    κ_ub_default = 5.0,
    max_iters = 50,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf,
    N_starts = 10,
    inner_xtol_rel = 1e-12,
    inner_maxeval = 1000,
    inner_maxtime = Inf);

println("obj_func3(β0) = ", obj_func3(β0))

updatedfunc(zeros(length(p_manual)))
updatedfunc(p_manual)
updateλfunc(p_manual)

# @btime run_optim(β0);
# @assert 5==3

println("Timing: run_optim")
p_β = copy(β0)

minf = NaN
minx = []
ret = :NaN
N_evals = -1
if typeof(optim_algorithm) == Symbol
    # using NLopt.
    @time minf, minx, ret, N_evals = run_optim(p_β)
else
    @time p = run_optim(p_β)

    minf = p.value
    minx = p.location
    ret = p.ret
end
println()

println("obj_func3(minx) = ", obj_func3(minx))
minx3 = copy(minx)
q3_star_U = q3.(U_rad)





PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
#PyPlot.plot(P, real.(q_manual_U), label = "q manual")
PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). refined vs data")





PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(g_star_U), "--", label = "z LS")
PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). refined vs LS z")

st_ind_d = 1
fin_ind_d = st_ind_d + N_d - 1

st_ind_β = fin_ind_d + 1
fin_ind_β = st_ind_β + N_β - 1

p3 = copy(p_manual)
p3[st_ind_β:fin_ind_β] = minx
refined_manual_cost = obj_func(p3)

println("last_version_cost   = ", last_version_cost)
println("refined_manual_cost = ", refined_manual_cost)





PyPlot.figure(fig_num)
fig_num += 1

#PyPlot.plot(P, real.(q_manual_U), label = "q manual")
PyPlot.plot(P_y, real.(y), label = "data")
PyPlot.plot(P, real.(q_star_U), "--", label = "last version fit")
PyPlot.plot(P, real.(q3_star_U), "--", label = "refined manual")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). refined vs last version (faulty) fit")
