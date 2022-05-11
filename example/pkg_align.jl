# packaged version of align.jl. Run align_prep.jl or prep_full.jl first.

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

import NLopt

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


save_flag = false


N_d = sum( NMRCalibrate.getNd(Bs[n]) for n = 1:length(Bs) )
N_β = sum( NMRCalibrate.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )

shift_lb = -ones(N_d)
shift_ub = ones(N_d)

Es = collect( NMRSpectraSimulator.καFIDModelType(Bs[i]) for i = 1:length(Bs) )

a_setp, b_setp, minxs,
    rets = NMRCalibrate.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#

w = ones(Float64, length(As))
#LS_inds = 1:length(U_cost)
loop_range = 1:length(cost_inds_set)

println("Timing:")
@time obj_funcs, minfs, minxs, rets = NMRCalibrate.aligncompound(y,
    U_y,
    P_y,
    As,
    Bs,
    Es,
    fs,
    SW,
    Δsys_cs,
    a_setp, b_setp, κs_β_DOFs, κs_β_orderings,
    shift_lb,
    shift_ub,
    cost_inds_set;
    loop_range = loop_range,
    w = w,
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-3,
    maxeval = 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    κ_lb_default = 1e-2,
    κ_ub_default = 1e2,
    β_max_iters = 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf);

if save_flag
    save_path = joinpath(project_folder, "alignment_results.bson")
    BSON.bson(save_path, region_min_dist = region_min_dist,
    minfs = minfs,
    minxs = minxs,
    rets = rets)
end
dummy = 1


#### visualize.
U = LinRange(u_min, u_max, 50000)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

for r in loop_range

    q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)

    obj_funcs[r](minxs[r])
    q_U = q2.(U_rad)

    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(P_y, real.(y), label = "y")
    PyPlot.plot(P, real.(q_U), label = "q")
    #PyPlot.plot(P, real.(q3_star_U), "--", label = "run_optim")

    PyPlot.legend()
    PyPlot.xlabel("ppm")
    PyPlot.ylabel("real")
    PyPlot.title("Alignment result: region $(r)")

end

for r in loop_range
    println("region $(r):")
    println("objective: $(minfs[r]), return status: $(rets)")
    println("shift variable:")
    display(minxs)
    println()
end
