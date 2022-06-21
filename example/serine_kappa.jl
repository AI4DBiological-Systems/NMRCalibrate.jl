

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

import NLopt
import Statistics
import Random

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])
#

# First, run prep_full.jl without optimization.

r = 1
T = Float64

y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
# LS_inds = 1:length(U_cost)

Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )

Es[1].κs_α[1][:] = [3.6547181133373705, 7.4585069055388304, 5.106122985644498]
#Es[1].core.λ0 = 2.3752439850648
#Es[1].core.ss_params.κs_λ[1] = 1.0
Es[1].core.ss_params.κs_β[1] = [1.7457371427914214, 2.9647985946238675, 2.5466218556915283] # calibration answer.

#Es[1].core.ss_params.d[1] = 5.4534685586024585 # calibration answer


#tmp = (1.7457371427914214 + 2.9647985946238675)/2
Es[1].core.ss_params.κs_β[1][2] = Es[1].core.ss_params.κs_β[1][1]
Es[1].κs_α[1][2] = Es[1].κs_α[1][1]

# # originals.
# As[1].Δc_bar[1][1] = [-0.04550785321752465, -0.9527933956628394, -0.0016987511196359661]
# As[1].Δc_bar[1][2] = [-0.9542311847361354, -0.04525616042734937, -0.0005126548365153027]
# As[1].Δc_bar[1][3] = [-0.0002609620463399076, -0.0019504439098113058, -0.9977885940438487]

## modify/test.
#As[1].Δc_bar[1][1] = As[1].Δc_bar[1][2]

# http://gissmo.nmrfam.wisc.edu/entry/bmse000867/simulation_1
# http://gissmo.nmrfam.wisc.edu/entry/bmse000048/simulation_1

βs = collect( dot(Es[1].core.ss_params.κs_β[1], As[1].Δc_bar[1][k]) for k = 1:length(As[1].Δc_bar[1]))
vs = collect( Es[1].core.ss_params.κs_β[1] .* As[1].Δc_bar[1][k] for k = 1:length(As[1].Δc_bar[1]))



# for calibration, set w to all ones.
q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = ones(length(Es)))


# eval.
N_viz = 50000
U = LinRange(u_min, u_max, N_viz)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

q_U = q2.(U_rad)

# viz.

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_U), label = "fit")
PyPlot.plot(P_cost, real.(y_cost), label = "data")
PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.gca().invert_xaxis()
PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data vs. bandpassed")
