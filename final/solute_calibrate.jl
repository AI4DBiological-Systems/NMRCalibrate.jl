

import NMRDataSetup
import NMRSpectraSimulator
import NMRSpecifyRegions

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JLD

#import Clustering
import Statistics
import Random

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.
# 0.1% DSS is 0.0046 M = 4.6 mM.
projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final"


## reboot. #######
project_name = "D-(+)-Glucose-700"
molecule_names = ["D-(+)-Glucose"; "DSS"]
w = [20.0/4.6; 1.0] # BMRB: DSS is 0.1 % => 4.6 mM

# project_name = "D-(+)-Glucose-NRC-600"
# molecule_names = ["D-(+)-Glucose";]
# w = [1.0; ]

# project_name = "L-Phenylalanine-700"
# molecule_names = ["L-Phenylalanine"; "DSS"]
# w = [20/0.5; 1.0] # BMRB-700 phenylalanine: DSS is 500 micro M.

# project_name = "L-Glutamine-700"
# molecule_names = ["L-Glutamine"; "DSS"]
# w = [20.0/0.5; 1.0] # BMRB: DSS is 500 uM => 0.5 mM

# project_name = "L-Glutamic acid-700"
# molecule_names = ["L-Glutamic acid"; "DSS"]
# w = [20.0/4.6; 1.0] # BMRB: DSS is 0.1% => 4.6 mM

# project_name = "L-Histidine-700"
# molecule_names = ["L-Histidine"; "DSS"]
# w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
#
# project_name = "L-Isoleucine-700"
# molecule_names = ["L-Isoleucine"; "DSS"]
# w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM

# project_name = "L-Serine-700"
# molecule_names = ["L-Serine"; "DSS"]
# w = [20.0/46; 1.0] # BMRB: DSS is 1 % => 46 mM
#
#
# project_name = "L-Alanine-700"
# molecule_names = ["L-Alanine"; "DSS"]
# w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
#
# project_name = "L-Threonine-700"
# molecule_names = ["L-Threonine"; "DSS"]
# w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
#
# project_name = "L-Tryptophan-700"
# molecule_names = ["L-Tryptophan"; "DSS"]
# w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM
#
# project_name = "L-Valine-700"
# molecule_names = ["L-Valine"; "DSS"]
# w = [20.0/0.5; 1.0] # BMRB: DSS is 500uM => 0.5 mM




w = w ./ norm(w) # since the fit data, y, is normalized.

# path to the GISSMO Julia storage folder.
base_path_JLD = "/home/roy/Documents/repo/NMRData//src/input/molecules"

# proxy-related.
tol_coherence = 1e-2
α_relative_threshold = 0.05
#Δc_partition_radius = 1e-1 # 17 element partition for glucose spin group 2.
Δc_partition_radius_candidates = [1e-1; 0.3; 0.5; 0.7; 0.8; 0.9]
#Δc_partition_radius_candidates = [0.9;]
λ0 = 3.4
Δcs_max = 0.2 # for proxy.
κ_λ_lb = 0.5
κ_λ_ub = 2.5

# if a Δc_partition_radius candidate gives max partition sizes less than of
#   equal to this value, then use that candidate.
early_exit_part_size = 7

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

PyPlot.plot(P_y, real.(S_U), label = "data spectrum")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data spectra")



####### mixture proxy.

function trydiffΔcradius(Δc_partition_radius_candidates::Vector{T},
    molecule_names, base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc,
    fs, SW, λ0, ν_0ppm, early_exit_part_size) where T <: Real

    @assert early_exit_part_size > 0
    @assert all(Δc_partition_radius_candidates .> zero(T))

    Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

    for Δc_partition_radius in Δc_partition_radius_candidates[1:end-1]

        mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
            base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0, ν_0ppm;
            tol_coherence = tol_coherence,
            α_relative_threshold = α_relative_threshold,
            Δc_partition_radius = Δc_partition_radius)
        As = mixture_params

        if all( all(NMRCalibrate.displaypartitionsizes(As[n]) .<= early_exit_part_size) for n = 1:length(As) )
            return mixture_params, Δc_partition_radius
        end
    end

    mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0, ν_0ppm;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius_candidates[end])

    return mixture_params, Δc_partition_radius_candidates[end]
end

Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

println("Timing: trydiffΔcradius()")
@time mixture_params, Δc_partition_radius = trydiffΔcradius(Δc_partition_radius_candidates,
    molecule_names, base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc,
    fs, SW, λ0, ν_0ppm, early_exit_part_size)
As = mixture_params

ΩS_ppm = NMRCalibrate.findfreqrange(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))

println("$(project_name): Partition sizes:")
display(NMRCalibrate.displaypartitionsizes(As[1]))
println("Δc_partition_radius = ", Δc_partition_radius)
println()

#@assert 1==2

# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(3.0)

u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

println("Timing: fitproxies!()")
@time NMRSpectraSimulator.fitproxies!(As;
#NMRSpectraSimulator.fitproxiessimple!(As;
    κ_λ_lb = κ_λ_lb,
    κ_λ_ub = κ_λ_ub,
    u_min = u_min,
    u_max = u_max,
    Δr = 1.0,
    Δκ_λ = 0.05)

### cost func.
combinevectors = NMRSpectraSimulator.combinevectors

#cs_config_path = "/home/roy/MEGAsync/inputs/NMR/configs/reduced_cs_config.txt"
cs_config_path = "/home/roy/Documents/repo/NMRData/src/input/reduced_cs_config.txt"


Δsys_cs, y_cost_all, U_cost_all, P_cost_all, exp_info, cost_inds,
cost_inds_set = NMRCalibrate.prepareoptim(cs_config_path, molecule_names, hz2ppmfunc,
U_y, y, As; region_min_dist = 0.1)

# visualize cost regions.
# sort(P_cost_set[1]) # etc..
P_cost_set = collect( P_y[cost_inds_set[r]] for r = 1:length(cost_inds_set) )

include("./final_helper.jl")


### optim all regions. # 900 secs.
y_cost = y_cost_all
U_cost = U_cost_all
P_cost = P_cost_all

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "data spectrum")
PyPlot.plot(P_cost, real.(y_cost), "^", label = "positions")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("positions against data spectrum, real part")

#@assert 1==2

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)

##### set up updates.

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
#ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )



## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
# purposely perturb κ.

Es = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:length(As) )


κ_lb_default = 0.2
κ_ub_default = 50.0

#@assert 1==2

## fit model.
println("Timing: calibrateregions()")
@time cost_inds_set, p_star_set, κ_BLS_set, d_star_set, β_star_set, λ_star_set,
proxies_set = calibrateregions(y, U_y, P_y, cost_inds_set,
Δ_shifts, As, fs, SW, w;
max_iters = 50000,
xtol_rel = 1e-3,
ftol_rel = 1e-6,
κ_lb_default = κ_lb_default,
κ_ub_default = κ_ub_default,
λ_each_lb = 0.9,
λ_each_ub = 1.1)

### save block.
save_path = joinpath(joinpath(projects_dir, project_name), "results_full.bson")
BSON.bson(save_path,
p_star_set = p_star_set,
κ_lb_default = κ_lb_default,
κ_ub_default = κ_ub_default,
κ_star_set = κ_BLS_set,
d_star_set = d_star_set,
β_star_set = β_star_set,
λ_star_set = λ_star_set,
cost_inds_set = cost_inds_set,
w = w,
# proxy setup-related below.
Δc_partition_radius = Δc_partition_radius,
tol_coherence = tol_coherence,
α_relative_threshold = α_relative_threshold,
λ0 = λ0,
Δcs_max = Δcs_max,
κ_λ_lb = κ_λ_lb,
κ_λ_ub = κ_λ_ub)
## end save block.


q_U_set = collect( proxies_set[r].(U) for r = 1:length(proxies_set) )
