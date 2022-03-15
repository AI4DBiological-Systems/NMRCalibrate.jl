

import NMRDataSetup
import NMRSpectraSimulator
# include("/home/roy/Documents/repo/NMRSpectraSimulator/src/NMRSpectraSimulator.jl")
# import .NMRSpectraSimulator

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
base_path_JLD = "/home/roy/Documents/repo/NMRData/src/input/molecules"

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

PyPlot.plot(P_y, real.(S_U), label = "data spectrum")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data spectra")



####### mixture proxy.


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

P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
#ΩS_ppm = collect( hz2ppmfunc.( NMRSpectraSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in mixture_params )


# lorentzian (oracle/ground truth)
f = uu->NMRSpectraSimulator.evalmixture(uu, mixture_params)

U_LS = U
w = ones(length(As))


## parameters that affect qs.
# A.ss_params.d, A.ss_params.κs_λ, A.ss_params.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
# purposely perturb κ.

As2 = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:length(As) )

Ag = As2[end]
#Ag.κ = collect( rand(length(Ag.κ[i])) for i = 1:length(Ag.κ) )
#Ag.κ_singlets = rand(length(Ag.κ_singlets))
Ag.κ[1][1] = 0.3
Ag.κ[1][2] = 0.7
Ag.κ[1][3] = 0.1
Ag.κ_singlets[1] = 0.4

# assemble proxy and compare against oracle at locations in U.
q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As2; w = w)

f_U = f.(U)
q_U = q.(U)

discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, real.(q_U), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")

# view the average Δc vector for each partition for DSS (last compound).
average_Δc_vectors_DSS = NMRCalibrate.viewaverageΔc(As[end])

############## test the shifting, β, and λ of DSS.

# manual for now. DSS, last molecule.
Bs = As[end:end]
# Bs[1].d
# As2[end].core.d

N_shifts = sum( length(Bs[n].d) + length(Bs[n].d_singlets) for n = 1:length(Bs) )
Δ_shifts = ones(N_shifts) .* 0.05
# w_lb = ones(1) .* 1e-5
# w_ub = ones(1) .* 1000.0
# end manual.

N_d = sum( getNd(Bs[n]) for n = 1:length(Bs) )

st_ind = 1
updatedfunc = pp->NMRCalibrate.updatemixtured!(Bs, pp, st_ind, fs, SW, Δ_shifts)

st_ind_β = N_d + 1
updateβfunc = pp->NMRCalibrate.updateβ!(Bs, pp, st_ind_β)
N_β = sum( NMRCalibrate.getNβ(Bs[n]) for n = 1:length(Bs) )

# λupdate.
st_ind_λ = st_ind_β + N_β
updateλfunc = pp->NMRCalibrate.updateλ!(Bs, pp, st_ind_λ)
N_λ = sum( NMRCalibrate.getNλ(Bs[n]) for n = 1:length(Bs) )

N_vars = N_d + N_β + N_λ


# reference, zero shift, phase.
p_test = zeros(N_vars)
p_test[1:st_ind_β-1] .= 0.0
p_test[st_ind_β:st_ind_λ-1] .= 0.0
p_test[st_ind_λ:end] .= 1.0
updatedfunc(p_test)
updateβfunc(p_test)
updateλfunc(p_test)

#q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Bs; w = w)
#q = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Bs)
q_U_ref = q.(U)

# shift.
p_test[1] = -1.0 # this is the spin group for DSS' small peaks.
#p_test[end-1:end] .= -π
p_test[st_ind_λ-1] = -π # this entry affects the phase for DSS' main peak.
p_test[end] = 2.0
updatedfunc(p_test)
updateβfunc(p_test)
updateλfunc(p_test)

q_U = q.(U)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_U_ref), label = "reference q")
PyPlot.plot(P, real.(q_U), label = "shifted")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("q")


# reset shifts and phases.
p_test[1:st_ind_β-1] .= 0.0
p_test[st_ind_β:st_ind_λ-1] .= 0.0
p_test[st_ind_λ:end] .= 1.0
updatedfunc(p_test)
updateβfunc(p_test)
updateλfunc(p_test)

# reset κ.
Ag = As2[end]
#Ag.κ = collect( rand(length(Ag.κ[i])) for i = 1:length(Ag.κ) )
#Ag.κ_singlets = rand(length(Ag.κ_singlets))
Ag.κ[1][1] = 1.0
Ag.κ[1][2] = 1.0
Ag.κ[1][3] = 1.0
Ag.κ_singlets[1] = 1.0
q_U_ref = q.(U)


####### perturb κ
N_κ, N_κ_singlets = NMRCalibrate.countκ(As2)
N_κ_vars = N_κ + N_κ_singlets
E_BLS, κ_BLS = NMRCalibrate.setupupdatew(length(U_LS), N_κ_vars)

κ_lb = ones(N_κ_vars) .* 0.2
κ_ub = ones(N_κ_vars) .* 50.0
fill!(κ_BLS, 12.0) # reset.

# update As2 with κ_BLS.
NMRCalibrate.parseκ!(As2, κ_BLS)

check_sys = collect( As2[i].κ for i = 1:length(As2) )
check_singlets = collect( As2[i].κ_singlets for i = 1:length(As2) )
# should be all 12's.

Ag = As2[end]
#Ag.κ = collect( rand(length(Ag.κ[i])) for i = 1:length(Ag.κ) )
#Ag.κ_singlets = rand(length(Ag.κ_singlets))
Ag.κ[1][1] = 2.0
Ag.κ[1][2] = 0.7
Ag.κ[1][3] = 1.5
Ag.κ_singlets[1] = 1.23
q_oracle = q.(U)

# least square solve on κ.
b_BLS = [real.(q_oracle); imag.(q_oracle)]
NMRCalibrate.updateκ!(E_BLS, b_BLS, κ_BLS, U_LS, As2, κ_lb, κ_ub)


# visualize.
q_U = q.(U)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_oracle), label = "oracle q")
PyPlot.plot(P, real.(q_U), "--", label = "LS kappa")
PyPlot.plot(P, real.(q_U_ref), "x", label = "reference q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("estimating kappa")
