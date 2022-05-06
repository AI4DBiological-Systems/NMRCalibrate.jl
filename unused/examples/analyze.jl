
PyPlot.close("all")

Random.seed!(25)

### load block.
load_path = joinpath(joinpath(projects_dir, project_name), "results.bson")
dict = BSON.load(load_path)
p_star = dict[:p_star]
κ_lb_default = dict[:κ_lb_default]
κ_ub_default = dict[:κ_ub_default]

As = collect( Es[n].core for n = 1:length(As))
q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
κ_BLS, getshiftfunc, getβfunc, getλfunc, N_vars_set = NMRCalibrate.setupcostcLshiftLS(Es, As, fs, SW, LS_inds, U_cost, y_cost, Δ_shifts; w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)

obj_func = pp->NMRCalibrate.costcLshift(U_cost, y_cost,
updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es, κ_BLS, q)

### end new.

### reference.

# reference, zero shift, phase.
p_test = zeros(sum(N_vars_set))
initial_cost = obj_func(p_test)
q_initial_U = q.(U)

p_test = copy(p_star)
final_cost = obj_func(p_test)
q_final_U = q.(U)


PyPlot.figure(fig_num)
fig_num += 1


PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("cost vs. data (y)")


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
