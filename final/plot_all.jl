
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


#projects_dir = save_folder_path
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-NRC-600"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/Nam2022_Serine"

### user inputs.
# projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final"
#
# project_name = "D-(+)-Glucose-700"
# molecule_names = ["D-(+)-Glucose"; "DSS"]
#
# project_name = "D-(+)-Glucose-NRC-600"
# molecule_names = ["D-(+)-Glucose";]
#
# project_name = "L-Phenylalanine-700"
# molecule_names = ["L-Phenylalanine"; "DSS"]
#
# project_name = "L-Glutamine-700"
# molecule_names = ["L-Glutamine"; "DSS"]
#
# project_name = "D-(+)-Glucose-NRC-600"
# molecule_names = ["D-(+)-Glucose";]

### end user inputs.

### load block.
#load_path = joinpath(joinpath(projects_dir, project_name), "results_full.bson")
load_path = joinpath(projects_dir, "results_full.bson")
dict = BSON.load(load_path)
As = collect( dict[:As][i] for i = 1:length(dict[:As]) )
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

function graphall(dict, Δ_shifts, Es, y, U_y, P_y, fs, SW::T; fig_num::Int = 1) where T <: Real

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

    for r = 1:length(cost_inds_set)

        y_cost = y[cost_inds_set[r]]
        U_cost = U_y[cost_inds_set[r]]
        P_cost = P_y[cost_inds_set[r]]

        LS_inds = 1:length(U_cost)

        q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
        κ_BLS, getshiftfunc, getβfunc, getλfunc,
        N_vars_set = NMRCalibrate.setupcostcLshiftLS(Es, As, fs, SW,
        LS_inds, U_cost, y_cost, Δ_shifts;
        w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)

        obj_func = pp->NMRCalibrate.costcLshift(U_cost, y_cost,
        updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es, κ_BLS, q)



        #####

        ### end new.

        ### reference.

        # reference, zero shift, phase.
        p_test = zeros(sum(N_vars_set))
        initial_cost = obj_func(p_test)
        q_initial_U = q.(U)

        p_test = copy(p_star_set[r])
        final_cost = obj_func(p_test)
        q_final_U = q.(U)


        PyPlot.figure(fig_num)
        fig_num += 1

        PyPlot.plot(P_y, real.(y), label = "y")
        PyPlot.plot(P_cost, real.(y_cost), "x")

        PyPlot.legend()
        PyPlot.xlabel("ppm")
        PyPlot.ylabel("real")
        PyPlot.title("r = $(r). cost vs. data (y)")



        PyPlot.figure(fig_num)
        fig_num += 1

        PyPlot.plot(P, real.(q_initial_U), label = "q initial")
        PyPlot.plot(P_y, real.(y), label = "y")
        PyPlot.plot(P, real.(q_final_U), "--", label = "q final")
        #PyPlot.plot(P_cost, real.(y_cost), "x")

        PyPlot.legend()
        PyPlot.xlabel("ppm")
        PyPlot.ylabel("real")
        PyPlot.title("r = $(r). data (y) vs. fit")

    end

end

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)


Es = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:length(As) )
graphall(dict, Δ_shifts, Es, y, U_y, P_y, fs, SW; fig_num = 1)
