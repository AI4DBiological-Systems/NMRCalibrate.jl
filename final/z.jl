

using LinearAlgebra, FFTW
import BSON, Statistics, PyPlot, Random

import NMRSpectraSimulator

include("../src/NMRCalibrate.jl")
import .NMRCalibrate

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

#import Clustering

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# load_path = "/home/roy/MEGAsync/inputs/NMR/debug/test_As.bson"
# dict = BSON.load(load_path)
# As = collect( dict[:As][i] for i = 1:length(dict[:As]) )



projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/test_glucose1/"
#projects_dir = "/home/roy/MEGAsync/outputs/NMR/calibrate/final/D-(+)-Glucose-700-r3/"
r = 3

T = Float64

### load block.
#load_path = joinpath(joinpath(projects_dir, project_name), "results_full.bson")
load_path = joinpath(projects_dir, "results_$(r).bson")
dict = BSON.load(load_path)

Es = collect( dict[:Es][i] for i = 1:length(dict[:Es]) )
As = collect( Es[n].core for n = 1:length(Es))

Δsys_cs = convert(Vector{Vector{Float64}}, dict[:Δsys_cs])
y = convert(Vector{Complex{Float64}}, dict[:y])
U_y = convert(Vector{Float64}, dict[:U_y])
SW = dict[:SW]
fs = dict[:fs]
ν_0ppm = dict[:ν_0ppm]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
P_y = hz2ppmfunc.(U_y)


ΩS_ppm = NMRCalibrate.findfreqrange(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSpectraSimulator.combinevectors(ΩS_ppm))
u_offset = 0.5
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)


P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)

p_star = dict[:p_star]
κ_lb_default = dict[:κ_lb_default]
κ_ub_default = dict[:κ_ub_default]
κ_star = dict[:κ_star]
d_star = dict[:d_star]
β_star = dict[:β_star]
λ_star = dict[:λ_star]
cost_inds_set = dict[:cost_inds_set]
w = dict[:w]


### end load block.

Δ_shifts = NMRSpectraSimulator.combinevectors(Δsys_cs)

y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]

LS_inds = 1:length(U_cost)

q, updatedfunc, updateβfunc, updateλfunc, updateκfunc,
κ_BLS, getshiftfunc, getβfunc, getλfunc,
N_vars_set = NMRCalibrate.setupcostcLshiftLS(Es, As, fs, SW,
LS_inds, U_cost, y_cost, Δ_shifts;
w = w, κ_lb_default = κ_lb_default, κ_ub_default = κ_ub_default)

obj_func = pp->NMRCalibrate.costcLshift(U_cost, y_cost,
updatedfunc, updateβfunc, updateλfunc, updateκfunc, pp, Es, κ_BLS, q)




#####

### end new.

### reference.

# reference, zero shift, phase.
N_d = sum( NMRCalibrate.getNd(As[n]) for n = 1:length(As) )
N_β = sum( NMRCalibrate.getNβ(As[n]) for n = 1:length(As) )
N_λ = sum( NMRCalibrate.getNλ(As[n]) for n = 1:length(As) )
#shift_manual = zeros(T, N_d)
#β_manual = zeros(T, N_β)
#λ_manual = ones(T, N_λ)

shift_manual = [-0.08; -0.1]
β_manual = copy(β_star)
λ_manual = copy(λ_star)

p_manual = [shift_manual; β_manual; λ_manual]

manual_cost = obj_func(p_manual)
#NMRCalibrate.parseκ!(Es, ones(T, length(κ_BLS)))
fill!(w, 1.0)
q_manual_U = q.(U)
println("norm(q_manual_U) = ", norm(q_manual_U))

# TODO I need to flip the phase a bit more. write optim code for just phase and lambda, given shift.


final_cost = obj_func(p_star)
q_star_U = q.(U)


PyPlot.figure(fig_num)
fig_num += 1

#PyPlot.plot(P, real.(q_manual_U), label = "q manual")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q_star_U), "--", label = "q star")
#PyPlot.plot(P_cost, real.(y_cost), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). data (y) vs. fit vs. manual")


############# LS fit on manual.

N_d = sum( NMRCalibrate.getNd(As[n]) for n = 1:length(As) )
N_β = sum( NMRCalibrate.getNβ(As[n]) for n = 1:length(As) )
N_λ = sum( NMRCalibrate.getNλ(As[n]) for n = 1:length(As) )

st_ind_d = 1
fin_ind_d = st_ind_d + N_d - 1

st_ind_β = fin_ind_d + 1
fin_ind_β = st_ind_β + N_β - 1

st_ind_λ = fin_ind_β + 1
fin_ind_λ = st_ind_λ + N_λ -1

NMRCalibrate.updatemixtured!(As, p_manual, st_ind_d, fs, SW, Δ_shifts)


#### the z idea.
# TODO do the κ idea.

Gs = collect( NMRSpectraSimulator.zCompoundFIDType(As[i]) for i = 1:length(As) )


β_lb = ones(T, N_β) .* (-π)
β_ub = ones(T, N_β) .* (π)
β_initial = zeros(T, N_β)



###
g, updateβfunc, updatezfunc,
z_BLS, getβfunc = NMRCalibrate.setupcostβLS(Gs, As, LS_inds, U_cost, y_cost)

println("Timing: updatezfunc() and parsez!()")
@time updatezfunc(0.0)
@time NMRCalibrate.parsez!(Gs, z_BLS)


g_star_U = g.(U)

# PyPlot.figure(fig_num)
# fig_num += 1
#
# PyPlot.plot(P, real.(q_manual_U), label = "q manual")
# PyPlot.plot(P_y, real.(y), label = "y")
# PyPlot.plot(P, real.(g_star_U), "--", label = "g star")
#
# PyPlot.legend()
# PyPlot.xlabel("ppm")
# PyPlot.ylabel("real")
# PyPlot.title("r = $(r). data (y) vs. fit vs. g")

###


function getmeanΔc(Δc_m_compound::Vector{Vector{Vector{T}}},
    part_inds_compound::Vector{Vector{Vector{Int}}}) where T

    @assert length(part_inds_compound) == length(Δc_m_compound)

    out = Vector{Vector{Vector{T}}}(undef, length(part_inds_compound))

    for i = 1:length(part_inds_compound)

        N_parition_elements = length(part_inds_compound[i])
        out[i] = Vector{Vector{T}}(undef, N_parition_elements)

        for k = 1:N_parition_elements
            inds = part_inds_compound[i][k]

            out[i][k] = Statistics.mean(Δc_m_compound[i][inds])
        end
    end

    return out
end

A = As[1]
Δc_m_compound_avg = getmeanΔc(A.Δc_m_compound, A.part_inds_compound)
#CNMRSpectraSimulator.combinevectors(Δc_m_compound_avg)

function estimateβ(inds_set::Vector{Vector{Int}},
    c_all::Vector{Vector{T}},
    phase_vec::Vector{T};
    st_ind::Int = 1) where T
    # loop through As, assemble. Δc matrix.
    #N_rows = sum( sum( length(As[n].part_inds_compound[i]) for i = 1:length(As[n].part_inds_compound) ) for n = 1:length(As) )

    @assert !isempty(c_all)
    N_partition_elements = length(inds_set)
    N_spins = length(c_all[1])

    C = Matrix{T}(undef, N_partition_elements, N_spins)
    y = Vector{T}(undef, N_partition_elements)

    for k = 1:N_partition_elements
        inds = inds_set[k]

        C[k,:] = Statistics.mean( c_all[inds] )
        y[k] = phase_vec[st_ind+k-1]
    end

    p_β = C\y

    fin_ind = st_ind + N_spins - 1

    return p_β, fin_ind, C, y
end

# z_flatten is the sum of all partition elements, of all spin groups in mixture.
# p_β is the sum of all β parameters, which is sum of all spins in each spin group in mixture.
function getinitialβ!(p_β::Vector{T},
    z_flatten::Vector{Complex{T}}) where T <: Real

    fill!(p_β, Inf)
    st_ind = 1

    for n = 1:length(As)
        A = As[n]

        for i = 1:length(A.part_inds_compound)
            p_β_i, fin_ind, C, y = estimateβ(A.part_inds_compound[i],
                A.Δc_m_compound[i], angle.(z_flatten); st_ind = st_ind)

            println("norm(C*p_β_i-y) = ", norm(C*p_β_i-y)) # debug.
            #@assert length(p_β_i) == length(st_ind:fin_ind)
            p_β[st_ind:fin_ind] = p_β_i
            st_ind = fin_ind + 1
        end
    end

    return nothing
end

β_z = similar(β_initial)
fill!(β_z, Inf)
getinitialβ!(β_z, z_BLS)
display([β_z; ])
clamp!(β_z, -π+1e-5, π -1e-5)


#optim_algorithm = :LN_BOBYQA
#optim_algorithm = :GN_ESCH

# initial guess

optim_algorithm = :LN_BOBYQA # good.
β0 = copy(β_initial)

optim_algorithm = :GN_ESCH
β0 = copy(β_initial) # good for 500.

optim_algorithm = :GN_ESCH
β0 = copy(β_z) # bad for 500, good for 5000.

optim_algorithm = :GN_ISRES
β0 = copy(β_z) # bad for 500. good for 5000

optim_algorithm = :GN_DIRECT_L
β0 = copy(β_z) # good for 500.


println("optim_algorithm = ", optim_algorithm)
println("β0 = ", β0)


### packaged up.
run_optim, f, κ_BLS, updateβfunc,
q3 = NMRCalibrate.setupβLSsolver(optim_algorithm,
    Es, As, LS_inds, U_cost, y_cost;
    κ_lb_default = 0.2,
    κ_ub_default = 5.0,
    max_iters = 50,
    xtol_rel = 1e-9,
    ftol_rel = 1e-9,
    maxtime = Inf);

println("Timing: run_optim")
p_β = copy(β0)
@time minf, minx, ret, N_evals = run_optim(p_β)




# force eval to update q2.
println("f(minx) = ", f(minx))
minx1 = copy(minx)
q3_star_U = q3.(U)

## unpackaged.
println("Timing: fitβ")
#β0 = β_initial
@time minx, q2, κ_BLS, getβfunc2,
obj_func2 = NMRCalibrate.fitβLSκ(U_cost, y_cost, LS_inds, Es, As,
    β0, β_lb, β_ub;
    optim_algorithm = optim_algorithm,
    max_iters = 500)

q2_star_U = q2.(U)
println("norm(q2_star_U-q3_star_U) = ", norm(q2_star_U-q3_star_U))
println()

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_manual_U), label = "q manual")
PyPlot.plot(P_y, real.(y), label = "y")
PyPlot.plot(P, real.(q2_star_U), "--", label = "q2 star")
PyPlot.plot(P, real.(g_star_U), "--", label = "g star")
PyPlot.plot(P, real.(q3_star_U), "--", label = "q3 star")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("r = $(r). data (y) vs. fit vs. q")



##### make q faster:
using BenchmarkTools

m = 1

@btime q3(U[m]) # 3 to 4 microsec.
@btime q2(U[m]) # 3 to 4 microsec.
@btime q(U[m]) # 3 to 4 microsec.

r0 = 2*π*U[m] - Es[1].core.ss_params.d[1]
@btime A.qs[1][1](r0)

# TODO have q(r,ξ,bb) use ss_params, to get q(rr) in NMRSpectraSimulator.
#   see if speed up.
