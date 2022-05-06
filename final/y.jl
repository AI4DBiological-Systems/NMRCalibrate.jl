

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

# load_path = "/home/roy/MEGAsync/inputs/NMR/debug/test_As.bson"
# dict = BSON.load(load_path)
# As = collect( dict[:As][i] for i = 1:length(dict[:As]) )



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


ΩS_ppm = NMRSpectraSimulator.getPsnospininfo(As, hz2ppmfunc)
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
# need to speed up eval design matrix.
N_κ, N_κ_singlets = NMRCalibrate.countκ(Es)


q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
E_BLS, κ_BLS, b_BLS, getshiftfunc, getβfunc, getλfunc,
N_vars_set = NMRCalibrate.setupcostcLshiftLS(Es, As, fs, SW,
LS_inds, U_rad_cost, y_cost, Δ_shifts;
w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)

obj_func = pp->NMRCalibrate.costcLshift(U_rad_cost, y_cost,
updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es,
E_BLS, κ_BLS, b_BLS, q)

last_version_cost = obj_func(p_star)
q_star_U = q.(U_rad)

using BenchmarkTools

# ### BLS timings.
# println("Timing: evaldesignmatrixκ!")
# Q = zeros(Complex{Float64}, length(U_cost), N_κ + N_κ_singlets)
# @btime NMRCalibrate.evaldesignmatrixκ!(Q, U_rad_cost, Es, w);
# # Q = zeros(Float64, 2*length(U_cost), N_κ + N_κ_singlets)
# # @btime NMRCalibrate.evaldesignmatrixκ!(Q, U_rad_cost, Es, w);
#
# # println("Timing: evaldesignmatrixκ!")
# # C = zeros(Float64, 2*length(U_cost), N_κ + N_κ_singlets)
# # @btime NMRCalibrate.evaldesignmatrixκ!(C, U_rad_cost, Es, w);
# #
# # println("descrepancy: ", norm([real.(Q); imag.(Q)] - C))
# # println()
#
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







### unused.
# ### packaged up.
# run_optim, obj_func3, E_BLS, κ_BLS, b_BLS, updateβfunc,
# q3 = NMRCalibrate.setupβLSsolver(optim_algorithm,
#     Es, As, LS_inds, U_rad_cost, y_cost;
#     κ_lb_default = 0.2,
#     κ_ub_default = 5.0,
#     max_iters = 50,
#     xtol_rel = 1e-9,
#     ftol_rel = 1e-9,
#     maxtime = Inf);
#
# println("Timing: run_optim")
# p_β = copy(β0)
# @time minf, minx, ret, N_evals = run_optim(p_β)
# println()
#
#
# #@assert 1==44
#
# # force eval to update q2.
# println("f(minx) = ", f(minx))
# minx1 = copy(minx)
# q3_star_U = q3.(U_rad)
