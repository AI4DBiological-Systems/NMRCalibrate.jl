# packaged version of nested_optim_d_only.jl

import NMRDataSetup
import NMRSpectraSimulator
import NMRSpecifyRegions

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JSON

import OffsetArrays
import Interpolations

import Statistics
import Random

import NLopt
import MultistartOptimization
import FiniteDiff

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# for convinence.
T = Float64
combinevectors = NMRSpectraSimulator.combinevectors

### user inputs.
region_min_dist = 0.1 # in ppm.

fit_config_path = "/home/roy/Documents/repo/NMRData/input/fit_configs/calibrate_700MHz_type1_select_compounds.json"


project_name = "NRC-glucose-2018"

project_base_folder = "/home/roy/MEGAsync/outputs/NMR/calibrate/NRC"
project_folder = joinpath(project_base_folder, project_name)
#mkpath(project_folder)
### end inputs.


### load block.
load_path = joinpath(project_folder, "proxy.bson")
dict = BSON.load(load_path)


u_min = dict[:u_min]
u_max = dict[:u_max]
y = convert(Vector{Complex{Float64}}, dict[:y])
U_y = convert(Vector{Float64}, dict[:U_y])
SW = dict[:SW]
fs = dict[:fs]
ν_0ppm  = dict[:ν_0ppm]
λ0 = dict[:λ0]
molecule_names = dict[:molecule_names]
#w = dict[:w]
w = ones(length(dict[:As]))

As = collect( dict[:As][i] for i = 1:length(dict[:As]) )
Bs = collect( dict[:Bs][i] for i = 1:length(dict[:Bs]) )

### end load.

### alternatively, load from experiment.


### prepare positions.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
P_y = hz2ppmfunc.(U_y)

ΩS0 = NMRSpectraSimulator.getΩS(As)
ΩS0_ppm = NMRSpectraSimulator.getPs(ΩS0, hz2ppmfunc)

Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds, cost_inds_set,
    λ_lbs, λ_ubs, κs_β_orderings,
    κs_β_DOFs = NMRCalibrate.prepareoptim(fit_config_path, molecule_names,
    hz2ppmfunc, U_y, y, As; region_min_dist = region_min_dist)


# visualize cost regions.project_folder
# sort(P_cost_set[1]) # etc..
P_cost_set = collect( P_y[cost_inds_set[r]] for r = 1:length(cost_inds_set) )

a_setp, b_setp, minxs,
    rets = NMRCalibrate.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#

#include("align.jl")
#include("pkg_align.jl")
