

# for GISSMO L-Phenylalanine.

# if save_BSON_flag
#     save_path = joinpath(project_folder, "alignment_results.bson")
#     BSON.bson(save_path, region_min_dist = region_min_dist,
#     minfs = minfs,
#     minxs = minxs,
#     rets = rets)
# end


function plotregion(P, U, q_fixed_U, q_free_U, P_y, y, P_cost, y_cost, display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = true,
    display_plot_flag = true,
    canvas_size = (1000, 400))

    # # reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.
    # inds, _ = NMRSpectraSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
    # P_display = P[inds]
    # U_display = U[inds]
    # q_U_display = q_U[inds]

    P_display = P
    U_display = U
    q_free_U_display = q_free_U
    q_fixed_U_display = q_fixed_U

    # plot.
    isdir(save_folder) || mkpath(save_folder)

    plots_save_path = joinpath(save_folder, file_name)
    #title_string = "$(project_name) alignment results, region $(r), real part"

    plot_obj = Plots.plot( P_display,
        real.(q_free_U_display),
        title = title_string,
        label = "free",
        seriestype = :line,
        ticks = :native,
        xlims = (P_display[1],P_display[end]),
        hover = P_display,
        linewidth = 4,
        xlabel = "ppm",
        ylabel = "real part of spectrum",
        xflip = true,
        size = canvas_size)

    Plots.plot!(plot_obj, P_display, real.(q_fixed_U_display), label = "fixed",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)

    Plots.plot!(plot_obj, P_cost, real.(y_cost), label = "fit data",
        markershape = :circle,
        seriestype = :scatter,
        xflip = true)

    if save_plot_flag
        Plots.savefig(plot_obj, plots_save_path)
    end

    if display_plot_flag
        display(plot_obj)
    end


    return nothing
end


N_viz = 50000
canvas_size = (1000, 400)

U = LinRange(u_min, u_max, N_viz)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)



#### visualize.
# minxs[1][1] = 0.000
# minxs[2][1] = 0.000
display_reduction_factor = 100
display_threshold_factor = 0.05/10
if "L-Isoleucine" in molecule_names
    display_reduction_factor = 1
    display_threshold_factor = 0.001/10
end
save_plot_flag = true
display_flag = true

#println("objective: $(minf), return status: $(rets)")
println("shift variable:")
display(minx)
println()

q2 = uu->NMRSpectraSimulator.evalitpproxymixture(uu, As, Es; w = w)

obj_func(minx)
q_free_U = q2.(U_rad)

# file_name = "fitted_kappa_alpha.html"
# title_string = "L-Phenylalanine, fitted κ_α, real part"

fill_val = Es[1].κs_α[1][1]
fill!(Es[1].κs_α[1], fill_val)
fill!(Es[1].κs_α[2], fill_val)
q_fixed_U = q2.(U_rad)

file_name = "kappa_alpha.html"
title_string = "$(molecule_names[1]), fixed vs. free κ_α, real part"

save_folder = project_folder


plotregion(P, U, q_fixed_U, q_free_U, P_y, y, P_cost, y_cost,
    display_threshold_factor, display_reduction_factor,
    save_folder, title_string, file_name;
    save_plot_flag = save_plot_flag,
    display_plot_flag = display_flag,
    canvas_size = canvas_size)
