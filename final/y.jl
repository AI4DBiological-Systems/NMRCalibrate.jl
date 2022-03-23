

using LinearAlgebra, FFTW
import BSON, Statistics, PyPlot, Random

import NMRSpectraSimulator

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

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


ΩS_ppm = NMRCalibrate.findfreqrange(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))
u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)


P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)

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


using BenchmarkTools

println("Timing: evaldesignmatrixκ!")
Q = zeros(Complex{Float64}, length(U_cost), N_κ + N_κ_singlets)
@btime NMRCalibrate.evaldesignmatrixκ!(Q, U_rad_cost, Es, w);
# Q = zeros(Float64, 2*length(U_cost), N_κ + N_κ_singlets)
# @btime NMRCalibrate.evaldesignmatrixκ!(Q, U_rad_cost, Es, w);

# println("Timing: evaldesignmatrixκ!")
# C = zeros(Float64, 2*length(U_cost), N_κ + N_κ_singlets)
# @btime NMRCalibrate.evaldesignmatrixκ!(C, U_rad_cost, Es, w);
#
# println("descrepancy: ", norm([real.(Q); imag.(Q)] - C))
# println()

println("q.(U_rad_cost)")
@btime q.(U_rad_cost);

println("obj_func.(p_star)")
@btime obj_func(p_star);

obj_func(p_star)

# TODO I am here.

q_U = q.(U_rad_cost)

b_BLS = b = collect(reinterpret(Float64, y_cost))

@btime B = reinterpret(Float64, Q);
@btime c = B*κ_BLS - b_BLS;

norm(q.(U_rad_cost) - y_cost)^2 - obj_func(p_star)

L = 50
a = randn(Complex{Float64}, L)
b = randn(Complex{Float64}, L)

sum( abs2(a[i]-b[i]) for i in eachindex(a))

sum( real(a[i]-b[i])^2 + imag(a[i]-b[i])^2 for i in eachindex(a) )

@assert 5==4

# complex to real, imag matrix, interlaced.
M1 = 5000
L1 = 50
W = randn(Complex{Float64}, M1, L1)
@btime Wr = reinterpret(Float64, W);
@btime Cr = [real.(W); imag.(W)];

# complex to real,imag array, interlaced.
V = randn(Complex{Float64}, M1)
@btime Vr = reinterpret(Float64, V);
@btime Br = [real.(V); imag.(V)];

Wr = reinterpret(Float64, W)
Vr = reinterpret(Float64, V)
@btime Wr'*Vr;

Cr = [real.(W); imag.(W)]
Br = [real.(V); imag.(V)]
@btime Cr'*Br;


# L = 500
# x = randn(L)
# y = randn(L)
# z = x + im .* y
#
# @btime sum(z-z);
# @btime sum(x-x) + sum(y-y);
