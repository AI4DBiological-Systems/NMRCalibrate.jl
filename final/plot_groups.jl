### I Am here. move this to NMRSpectraSimulator.

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

import Destruct

include("./helpers/final_helpers.jl")
include("./helpers/resonance_helpers.jl")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])



#projects_dir = save_folder_path
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Serine-700"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-NRC-600"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-700"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Isoleucine-700"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Leucine-500-2mM"
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Glutamine-700"

plots_save_folder = joinpath(projects_dir, "plots")
isdir(plots_save_folder) || mkdir(plots_save_folder)
project_title = "resonance groupings"




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
# A.ss_params.d, A.ss_params.κs_λ, A.ss_params.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As[1:1])

#g = uu->evalitpproxycompoundresonance(uu, As[1])
A = As[1]
g1 = uu->evalitpproxysysresonance(A.qs, uu, A.ss_params.d, A.ss_params.κs_λ, A.ss_params.κs_β)
g2 = uu->evalsingletsresonance(uu, A.d_singlets, A.αs_singlets, A.Ωs_singlets,
A.β_singlets, A.λ0, A.κs_λ_singlets)

f_U = f.(U)
q_U = q.(U)
g1_U = g1.(U)
g2_U = g2.(U)

#g_qs_U, g_singlets_U = convertresonancetimeseries(g_U)


z1 = Destruct.destruct(g1_U);
q_singlets_evals = Destruct.destruct(g2_U);


# collect( length(z1[i]) for i = 1:length(z1))
# collect( length(z2[i]) for i = 1:length(z2))
qs_evals = assembleqsevals(z1, A.qs)


sanity_check = zeros(Complex{Float64}, length(q_U))
if !isempty(z1)
    sanity_check += sum(z1)
end
if !isempty(q_singlets_evals)
    sanity_check += sum(q_singlets_evals)
end
println("sanity: ", norm(sanity_check - q_U))
@assert norm(sanity_check - q_U) < 1e-14
println()



discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

## visualize.
PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, real.(q_U), label = "q")

#visualizesubgroups(P, qs_evals, q_singlets_evals, real)

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")

#@assert 1==2

plots_save_path = joinpath(plots_save_folder, "groups_real.html")
title_string = "$(project_title), real"
savefiggroups(plots_save_path, title_string,
    P, q_U, qs_evals, q_singlets_evals, real)

#
plots_save_path = joinpath(plots_save_folder, "groups_imag.html")
title_string = "$(project_title), imag"
savefiggroups(plots_save_path, title_string,
    P, q_U, qs_evals, q_singlets_evals, real)
#
plots_save_path = joinpath(plots_save_folder, "groups_modulus.html")
title_string = "$(project_title), modulus"
savefiggroups(plots_save_path, title_string,
    P, q_U, qs_evals, q_singlets_evals, real)
