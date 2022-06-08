# packaged version of align.jl. Run align_prep.jl or prep_full.jl first.

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

import NLopt
import PlotlyJS
using Plots; plotly()


include("helper.jl")

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

#LS_inds = 1:length(U_cost)
loop_range = 1:length(cost_inds_set)

println("Timing:")
@time obj_funcs, minfs, minxs, rets, ws = NMRCalibrate.alignquantificationcompound(y,
    U_y,
    P_y,
    As,
    Bs,
    Es,
    fs,
    SW,
    Δsys_cs,
    a_setp, b_setp, #κs_β_DOFs, κs_β_orderings,
    shift_lb,
    shift_ub,
    cost_inds_set;
    loop_range = loop_range,
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-9,
    maxeval = 50, # 2, # 50,
    maxtime = Inf,
    β_optim_algorithm = :GN_DIRECT_L,
    w_lb_default = 1e-1,
    w_ub_default = 100.0,
    β_max_iters = 500, # 2, # 500,
    β_xtol_rel = 1e-9,
    β_ftol_rel = 1e-9,
    β_maxtime = Inf);

dummy = 1

#minxs[3][1] = 0.0

#@assert 99==44

if save_BSON_flag
    save_path = joinpath(project_folder, "quantify_results.bson")
    BSON.bson(save_path, region_min_dist = region_min_dist,
    minfs = minfs,
    minxs = minxs,
    rets = rets)
end

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

# save data plot.
file_name = "data.html"
plots_save_path = joinpath(project_folder, file_name)



plotquantificationresults(As, Es, ws, project_folder,
    P_y, y,
    region_min_dist,
    obj_funcs, minfs, minxs, rets,
    display_reduction_factor, display_threshold_factor,
    cost_inds_set, loop_range;
    canvas_size = (1000, 400),
    display_flag = display_flag,
    save_plot_flag = save_plot_flag,
    N_viz = 50000)

for r in loop_range
    println("region $(r):")
    println("objective: $(minfs[r]), return status: $(rets)")
    println("shift variable:")
    display(minxs)
    println()
end
