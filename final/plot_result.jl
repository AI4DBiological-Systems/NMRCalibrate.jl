
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

#projects_dir = save_folder_path
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-NRC-600"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/Nam2022_Serine"
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/L-Serine-700"

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
plots_save_folder = joinpath(projects_dir, "plots")
isdir(plots_save_folder) || mkdir(plots_save_folder)
project_title = "test" # TODO change later.

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

function graphall(dict, Δ_shifts, Es, y, U_y, P_y, fs, SW::T,
    project_title::String,
    plots_save_folder::String;
    fig_num::Int = 1) where T <: Real

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
        N_d = sum( length(As[n].ss_params.d) + length(As[n].d_singlets) for n = 1:length(As) )
        N_β = sum( NMRCalibrate.getNβ(As[n]) for n = 1:length(As) )
        N_λ = sum( NMRCalibrate.getNλ(As[n]) for n = 1:length(As) )
        shift_initial = zeros(T, N_d)
        β_initial = zeros(T, N_β)
        λ_initial = ones(T, N_λ)
        p_initial = [shift_initial; β_initial; λ_initial]

        initial_cost = obj_func(p_initial)
        NMRCalibrate.parseκ!(Es, ones(T, length(κ_BLS)))
        fill!(w, 1.0)
        q_initial_U = q.(U)
        #println("norm(q_initial_U) = ", norm(q_initial_U))

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

        plots_save_path = joinpath(plots_save_folder, "region_$(r)_real.html")
        title_string = "$(project_title): region $(r), real"
        savefigfitresult(plots_save_path, title_string,
         real.(q_final_U), P, P_cost, real.(y_cost);
        initial_fit = real.(q_initial_U))

        plots_save_path = joinpath(plots_save_folder, "region_$(r)_imaginary.html")
        title_string = "$(project_title): region $(r), imaginary"
        savefigfitresult(plots_save_path, title_string,
        imag.(q_final_U), P, P_cost, imag.(y_cost);
        initial_fit = imag.(q_initial_U))

        plots_save_path = joinpath(plots_save_folder, "region_$(r)_modulus.html")
        title_string = "$(project_title): region $(r), modulus"
        savefigfitresult(plots_save_path, title_string,
        abs.(q_final_U), P, P_cost, abs.(y_cost);
        initial_fit = abs.(q_initial_U))
    end

end

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)


Es = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:length(As) )
graphall(dict, Δ_shifts, Es, y, U_y, P_y, fs, SW,
project_title, plots_save_folder; fig_num = 1)
