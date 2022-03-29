"""
window ∈ (0,1)
optim_algorithm can be :GN_ESCH, :GN_ISRES, :LN_BOBYQA, :GN_DIRECT_L
"""
function setupitpab(window::T, N_itp_samples::Int, input_range_percentage::T;
    p0 = [0.5; 0.0],
    p_lb = [0.1; -5.0],
    p_ub = [0.6; 5.0],
    max_iters = 5000,
    xtol_rel = 1e-5,
    ftol_rel = 1e-5,
    maxtime = Inf,
    optim_algorithm = :LN_BOBYQA) where T <: Real

    # get piece-wise linear monotone maps.
    infos, zs, p_range = MonotoneMaps.prepareboxboundwarping(zero(T), one(T), window; N_itp_samples = N_itp_samples, input_range_percentage = input_range_percentage)

    # get compact sigmoid parameters fitted to each of the piece-wise linear maps.
    costfuncs, minxs, rets = MonotoneMaps.getcompactsigmoidparameters(infos;
    p0 = p0, p_lb = p_lb, p_ub = p_ub, optim_algorithm = optim_algorithm, max_iters = 5000)
    #qs = collect( tt->MonotoneMaps.evalcompositelogisticprobit(tt, minxs[i][1], minxs[i][2]) for i = 1:length(minxs) )

    Δp = p_range[2]-p_range[1]
    itp_range = p_range[1]:Δp:p_range[end]

    # itp.
    a_samples = collect( minxs[i][1] for i = 1:length(minxs) )
    b_samples = collect( minxs[i][2] for i = 1:length(minxs) )

    a_itp = Interpolations.interpolate(a_samples, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    a_sitp = Interpolations.scale(a_itp, itp_range)
    a_setp = Interpolations.extrapolate(a_sitp, Interpolations.Flat()) # zero outside interp range.

    b_itp = Interpolations.interpolate(b_samples, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    b_sitp = Interpolations.scale(b_itp, itp_range)
    b_setp = Interpolations.extrapolate(b_sitp, Interpolations.Flat()) # zero outside interp range.

    return a_setp, b_setp,
        minxs, rets # debug
end
