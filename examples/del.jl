

PyPlot.close("all")

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_initial_U), label = "q initial")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_final_U), "--", label = "q final")
#PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data (y) vs. fit")

#κ2 = κ_BLS .* 2
final_cost = obj_func(p_star)
κ2 = copy(κ_BLS)
# κ2[1] *= 1
# κ2[2] *= 1
# κ2[2] *= 1
NMRCalibrate.parseκ!(Es, κ_BLS)

q2_evals = q.(U_y)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P_y, real.(q2_evals), label = "q2")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data (y) vs. q2")


f2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, Es; w = w)
f2_evals = f2.(U_y)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P_y, real.(f2_evals), label = "f2")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("data (y) vs. f2")

@assert 1==2

# view the average Δc vector for each partition for DSS (last compound).
#average_Δc_vectors_DSS = NMRCalibrate.viewaverageΔc(As[end])

y21 = y .* exp(-im*Es[1].core.κs_β[1][1])
y22 = y .* exp(-im*Es[1].core.κs_β[2][1])

# visualize.
#q_U = q.(U)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(y21), label = "y2")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("y2")

# just optim the phase offline, since no overlap!

ΩS0 = NMRSpecifyRegions.getΩS(As)
ΩS0_ppm = NMRSpecifyRegions.getPs(ΩS0, hz2ppmfunc)




q_U21 = q_final_U .* exp(-im*Es[1].core.κs_β[1][1])
q_U22 = q_final_U .* exp(-im*Es[1].core.κs_β[2][1])

q_U2 = q_final_U .* exp(-im*β_star[1])
#q_U2 = q_final_U .* exp(-π*2/4)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_U2), label = "q2")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("q2")
