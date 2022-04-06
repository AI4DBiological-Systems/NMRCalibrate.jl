

import NMRDataSetup
import NMRSpectraSimulator

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

using LinearAlgebra
using FFTW
import PyPlot
import BSON
#import JLD

#import Clustering
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
base_path_JLD = "/home/roy/Documents/repo/NMRData//src/input/molecules"

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
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(N, fs, offset_Hz)

Z = maximum(abs.(DFT_s))
y = DFT_s[U_inds] ./ Z

##############################

S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(S_U), label = "data spectrum")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data spectra")

#@assert 1==2


Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

mixture_params = NMRSpectraSimulator.setupmixtureproxies(molecule_names,
    base_path_JLD, Δcs_max_mixture, hz2ppmfunc, ppm2hzfunc, fs, SW, λ0,
    ν_0ppm, dummy_SSFID;
    tol_coherence = tol_coherence,
    α_relative_threshold = α_relative_threshold)
As = mixture_params



u_min = ppm2hzfunc(-0.5)
u_max = ppm2hzfunc(3.0)

NMRSpectraSimulator.fitproxies!(As;
    κ_λ_lb = κ_λ_lb,
    κ_λ_ub = κ_λ_ub,
    u_min = u_min,
    u_max = u_max,
    Δr = 1.0,
    Δκ_λ = 0.05)


### plot.


f = uu->NMRSpectraSimulator.evalmixture(uu, mixture_params)

# test params.
ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )
#ΩS_ppm_flat = NMRSpectraSimulator.combinevectors(ΩS_ppm)




P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)

## parameters that affect qs.
# A.ss_params.d, A.ss_params.κs_λ, A.ss_params.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, mixture_params)

f_U = f.(U)
q_U = q.(U)

discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

## visualize.
PyPlot.figure(fig_num)
fig_num += 1


PyPlot.plot(P_y, real.(y), label = "data spectrum")
PyPlot.plot(P, real.(q_U), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data vs q")



# manual for now. DSS, last molecule.
Bs = As[end:end]
Δ_shifts = ones(2) .* 0.05
w_lb = ones(1) .* 1e-5
w_ub = ones(1) .* 1000.0
# end manual.

N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )

st_ind = 1
updatedfunc = pp->NMRCalibrate.updatemixtured!(Bs, pp, st_ind, fs, SW, Δ_shifts)

st_ind_β = N_d + 1
updateβfunc = pp->NMRCalibrate.updateβ!(Bs, pp, st_ind_β)
N_β = sum( NMRCalibrate.getNβ(Bs[n]) for n = 1:length(Bs) )
NMRCalibrate.getNβ(Bs[1])

N_vars = N_d + N_β
A_BLS, w_BLS = NMRCalibrate.setupupdatew(length(U), length(Bs))

# reference, zero shift, phase.
p_test = zeros(N_vars) # reset.
fill!(w_BLS, 1.0) # reset.

updatedfunc(p_test)
updateβfunc(p_test)

q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Bs; w = w_BLS)
q_U_ref = q.(U)


# shift.
p_test[1] = -1.0 # test.
p_test[end] = -π
updatedfunc(p_test)
updateβfunc(p_test)

##tmp = y[LS_inds]
#tmp = q.(U) .* 5
tmp = q.(U)
U_LS = U
b_BLS = [real.(tmp); imag.(tmp)]
NMRCalibrate.updatew!(A_BLS, b_BLS, w_BLS, U_LS, Bs, w_lb, w_ub)


q_U = q.(U)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_U_ref), label = "reference")
PyPlot.plot(P, real.(q_U), label = "shifted")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("q")
