
using LinearAlgebra, FFTW
import BSON, Statistics, Random
import PyPlot

import NMRSpectraSimulator

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

import PlotlyJS
import Plots
Plots.plotly()

#import Clustering

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


include("./helpers/final_helpers.jl")
include("./helpers/resonance_helpers.jl")

#projects_dir = save_folder_path
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-NRC-600"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/Nam2022_Serine"
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Serine-700"




projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-NRC-600"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/Nam2022_Serine"
#project_dit = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Serine-700"

T = Float64

### load block.
#load_path = joinpath(joinpath(projects_dir, project_name), "results_full.bson")
load_path = joinpath(projects_dir, "results_full.bson")
dict = BSON.load(load_path)
As = collect( dict[:As][i] for i = 1:length(dict[:As]) )
Es = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:length(As) )

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

p_star_set = dict[:p_star_set]
κ_lb_default = dict[:κ_lb_default]
κ_ub_default = dict[:κ_ub_default]
κ_star_set = dict[:κ_star_set]
d_star_set = dict[:d_star_set]
β_star_set = dict[:β_star_set]
λ_star_set = dict[:λ_star_set]
cost_inds_set = dict[:cost_inds_set]
w = dict[:w]

As = collect( Es[n].core for n = 1:length(Es))
### end load block.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)




### plot.


f = uu->NMRSpectraSimulator.evalmixture(uu, As)

# test params.
ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in As )
#ΩS_ppm_flat = NMRSpectraSimulator.combinevectors(ΩS_ppm)




P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)

## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As)

g = uu->evalitpproxycompoundresonance(uu, As[1])

f_U = f.(U)
q_U = q.(U)
g_U = g.(U)

g_qs_U, g_singlets_U = convertresonancetimeseries(g_U)

import Destruct
x, y = Destruct.destruct(g_U);
# I am here. think about how to use Destruct in resonance_helpers.jl


@assert 1==2

discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

## visualize.
PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, real.(q_U), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")
