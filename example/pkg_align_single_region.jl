
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

println("Timing:")
@time obj_func, minf, minx, ret, P_cost, y_cost = NMRCalibrate.aligncompoundsingleregion(y,
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
    w = w,
    N_starts = 100,
    local_optim_algorithm = NLopt.LN_BOBYQA,
    #xtol_rel = 1e-3,
    xtol_rel = 1e-9,
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
