# packaged version of align.jl. Run align_prep.jl or prep_full.jl first.

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

import NLopt
import PlotlyJS
using Plots; plotly()

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

save_BSON_flag = true

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

dummy = 1

if save_BSON_flag
    save_path = joinpath(project_folder, "alignment_results.bson")
    BSON.bson(save_path, region_min_dist = region_min_dist,
    minfs = minfs,
    minxs = minxs,
    rets = rets)
end


function plotalignmentresults(As, Es, w, save_folder,
    P_y, y,
    region_min_dist,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor;
    canvas_size = (1000, 400),
    display_flag = false,
    save_plot_flag = true,
    N_viz = 50000)

    U = LinRange(u_min, u_max, N_viz)
    P = hz2ppmfunc.(U)
    U_rad = U .* (2*π)

    for r in loop_range

        q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)

        obj_funcs[r](minxs[r])
        q_U = q2.(U_rad)

        # reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.
        inds, _ = NMRSpectraSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
        P_display = P[inds]
        U_display = U[inds]
        q_U_display = q_U[inds]

        # plot.
        isdir(save_folder) || mkpath(save_folder)

        plots_save_path = joinpath(save_folder, "results_real.html")
        title_string = "$(project_name) alignment results, region $(r), real part"

        plot_obj = Plots.plot( P_display,
            real.(q_U_display),
            title = title_string,
            label = "model",
            seriestype = :line,
            ticks = :native,
            xlims = (P_display[1],P_display[end]),
            hover = P_display,
            linewidth = 4,
            xlabel = "ppm",
            ylabel = "real part of spectrum",
            xflip = true,
            size = canvas_size)

        Plots.plot!(plot_obj, P_y, real.(y), label = "data",
            seriestype = :line,
            linestyle = :dot,
            xflip = true,
            linewidth = 4)


        Plots.savefig(plot_obj, plots_save_path)

        if save_plot_flag
            display(plot_obj)
        end

    end

    return nothing
end


#### visualize.
display_reduction_factor = 100
display_threshold_factor = 0.05/10
save_plot_flag = true
display_flag = true

plotalignmentresults(As, Es, w, project_folder,
    P_y, y,
    region_min_dist,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor;
    canvas_size = (1000, 400),
    display_flag = true,
    save_plot_flag = true,
    N_viz = 50000)

for r in loop_range
    println("region $(r):")
    println("objective: $(minfs[r]), return status: $(rets)")
    println("shift variable:")
    display(minxs)
    println()
end
