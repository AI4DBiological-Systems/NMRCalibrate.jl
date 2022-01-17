

import NMRDataSetup
import NMRSpectraSimulator
import NMRCalibrate

# include("../src/NMRSpecifyRegions.jl")
# import .NMRSpecifyRegions
import NMRSpecifyRegions

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import Statistics

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/"

# project_name = "glucose-700"
# molecule_names = ["D-(+)-Glucose"; "DSS"]
# w = [20.0/4.6; 1.0] # BMRB-700 glucose: DSS is 0.0046 M = 4.6 mM.

project_name = "phenylalanine-700"
molecule_names = ["L-Phenylalanine"; "DSS"]
w = [20/0.5; 1.0] # BMRB-700 phenylalanine: DSS is 500 micro M.

# path to the GISSMO Julia storage folder.
base_path_JLD = "/home/roy/Documents/data/NMR/NMRData/src/input/molecules"

# proxy-related.
tol_coherence = 1e-2
α_relative_threshold = 0.05
λ0 = 3.4
Δcs_max = 0.2
κ_λ_lb = 0.5
κ_λ_ub = 2.5

### end inputs.


## load data.
load_folder_path = joinpath(projects_dir, project_name)
load_path = joinpath(load_folder_path, "$(project_name).bson")
dic = BSON.load(load_path)
s_t = dic[:s_t]
fs = dic[:fs]
SW = dic[:SW]
ν_0ppm = dic[:ν_0ppm]

# normalize.
s_t = s_t

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRCalibrate.getwraparoundDFTfreqs(N, fs, offset_Hz)

Z = maximum(abs.(DFT_s))
y = DFT_s[U_inds] ./ Z


#
S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "data spectrum, y")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data spectra, real part")



####### mixture proxy.


Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold)
As = mixture_params


ΩS_ppm = NMRCalibrate.findfreqrange(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(3.0)

u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

NMRSpectraSimulator.fitproxies!(As;
    κ_λ_lb = κ_λ_lb,
    κ_λ_ub = κ_λ_ub,
    u_min = u_min,
    u_max = u_max,
    Δr = 1.0,
    Δκ_λ = 0.05)


### plot.

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
#ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )






cs_config_path = "/home/roy/MEGAsync/inputs/NMR/configs/cs_config_reduced.txt"

cs_delta_group = NMRSpecifyRegions.extractinfofromconfig( cs_config_path, molecule_names)
Δsys_cs = NMRSpecifyRegions.condenseΔcsconfig(cs_delta_group)

ΩS0 = NMRSpecifyRegions.getΩS(As)
ΩS0_ppm = NMRSpecifyRegions.getPs(ΩS0, hz2ppmfunc)
exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_names, ΩS0_ppm, Δsys_cs;
min_dist = 0.1)



function testfunc(exp_info, P_cost0)

    band_inds = Vector{Int}(undef, 0)

    for i = 1:length(exp_info.regions)

        tmp = NMRSpecifyRegions.filterfreqpositions(P_cost0, [exp_info.regions[i].st;], [exp_info.regions[i].fin;])
        push!(band_inds, tmp...)
    end

    unique!(band_inds)
    return band_inds
end

U_cost0 = U_y
P_cost0 = hz2ppmfunc.(U_cost0)
y_cost0 = y

band_inds = testfunc(exp_info, P_cost0)

U_cost = U_cost0[band_inds]
P_cost = P_cost0[band_inds]

y_cost = y_cost0[band_inds]


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "data spectrum")
PyPlot.plot(P_cost, real.(y_cost), "^", label = "positions")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("positions against data spectrum, real part")
