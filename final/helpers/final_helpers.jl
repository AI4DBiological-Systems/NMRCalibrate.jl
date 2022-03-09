
function calibrateregions( y::Vector{Complex{T}},
    U_y, P_y, cost_inds_set, Δ_shifts, As, fs, SW::T, w;
    max_iters = 5000,
    xtol_rel = 1e-7,
    ftol_rel = 1e-12,
    κ_lb_default = 0.2,
    κ_ub_default = 50.0,
    λ_each_lb = 0.9,
    λ_each_ub = 1.1 ) where T <: Real

    #
    N_compounds = length(As)
    N_regions = length(cost_inds_set)

    Es = collect( NMRSpectraSimulator.κCompoundFIDType(As[i]) for i = 1:N_compounds )

    p_star_set = Vector{Vector{T}}(undef, N_regions)
    κ_BLS_set = Vector{Vector{T}}(undef, N_regions)

    d_star_set = Vector{Vector{T}}(undef, N_regions)
    β_star_set = Vector{Vector{T}}(undef, N_regions)
    λ_star_set = Vector{Vector{T}}(undef, N_regions)
    proxies_set = Vector{Function}(undef, N_regions)

    for r = 1:N_regions

        y_cost = y[cost_inds_set[r]]
        U_cost = U_y[cost_inds_set[r]]
        P_cost = P_y[cost_inds_set[r]]

        ###
        LS_inds = 1:length(U_cost)

        println("Timing: runalignment(), r = $(r)")
        @time p_star, q, κ_BLS, getshiftfunc, getβfunc, getλfunc,
        obj_func, N_vars_set = NMRCalibrate.runalignment(Δ_shifts,
        U_cost, y_cost, LS_inds, Es, As, fs, SW;
        max_iters = max_iters,
        xtol_rel = xtol_rel,
        ftol_rel = ftol_rel,
        w = w,
        κ_lb_default = κ_lb_default,
        κ_ub_default = κ_ub_default,
        λ_each_lb = λ_each_lb,
        λ_each_ub = λ_each_ub)

        p_star_set[r] = p_star
        κ_BLS_set[r] = κ_BLS
        proxies_set[r] = q

        d_star_set[r] = getshiftfunc(p_star)
        β_star_set[r] = getβfunc(p_star)
        λ_star_set[r] = getλfunc(p_star)
    end

    return cost_inds_set, p_star_set, κ_BLS_set, d_star_set, β_star_set, λ_star_set, proxies_set
end


function savefigfitresult(save_path::String,
    title_string::String,
    q_U::Vector{T},
    P,
    P_cost,
    y_cost::Vector{T};
    initial_fit::Vector{T} = zeros(T, 0),
    P_y = P_cost,
    y = y_cost) where T

    plot_obj = Plots.plot( P_y,
        y,
        title = title_string,
        label = "Data",
        seriestype = :line,
        ticks = :native,
        xlims = (P[1],P[end]),
        hover = P,
        linewidth = 4,
        xflip = true,
        size = (1600, 900))

    Plots.plot!(plot_obj, P, q_U, label = "Fitted model",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)

    if !isempty(initial_fit)
        Plots.plot!(plot_obj, P, initial_fit, label = "Initial model",
            seriestype = :line,
            linestyle = :dash,
            xflip = true,
            linewidth = 4)
    end

    if !isempty(y_cost)
        Plots.plot!(plot_obj, P_cost, y_cost, label = "Fit positions",
            markershape = :circle,
            seriestype = :scatter)
    end

    Plots.savefig(plot_obj, save_path)
end
