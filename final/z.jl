

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


### nested optim.
