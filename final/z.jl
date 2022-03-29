import MultistartOptimization, NLopt

f(x) = sum(x -> abs2(x - 1), x)                     # objecive ∑(xᵢ-1)²
P = MultistartOptimization.MinimizationProblem(f, fill(-2, 4), fill(2, 4)) # search in [-2, 2]⁴


#local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = NLopt.GN_ESCH,
local_method = MultistartOptimization.NLoptLocalMethod(; algorithm = NLopt.LN_BOBYQA,
    xtol_rel = 1e-12,
    maxeval = 10,
    maxtime = Inf) # loaded on demand with `using NLopt`


multistart_method = MultistartOptimization.TikTak(100)
p = MultistartOptimization.multistart_minimization(multistart_method,
    local_method, P)
