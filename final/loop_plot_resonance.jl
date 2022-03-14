# not working.

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
include("./helpers/plot_resonance.jl")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])





# base_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/"
# project_set = filter(isdir, readdir(base_dir; join=true)) # all dirs in base_dir

project_set = Vector{String}(undef, 0)
tmp = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-NRC-600"
push!(project_set, tmp)

function loopplotresonance(project_set)

    for i = 1:length(project_set)
        projects_dir = project_set[i]

        println("Now on $(projects_dir)")
        plotresonancegroups(projects_dir)
        println()
    end

end


### batch.
loopplotresonance(project_set)
### end batch.
