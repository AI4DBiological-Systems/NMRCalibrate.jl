function evalsingletsresonance(u::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T}) where T <: Real

    u_rad = 2*π*u

    out = zeros(Complex{T}, length(αs_singlets))
    for i = 1:length(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out[i] = αs_singlets[i]*exp(im*βs_singlets[i])/(λ+im*(τ-Ω))
    end

    #return out # type stable.
    return tuple(out...) # type unstable.
end

# ::Vector{Vector{Complex{T}}} # if using type stable.
# type unstable.
function evalitpproxysysresonance(qs::Vector{Vector{Function}},
    u::T, d::Vector{T}, κs_λ::Vector{T},
    κs_β::Vector{Vector{T}}) where T <: Real

    @assert length(d) == length(qs)

    #out = collect( zeros(Complex{T}, length(qs[i])) for i = 1:length(qs) ) # type stable.

    N = sum( length(qs[i]) for i = 1:length(qs) )
    out = zeros(Complex{T}, N)
    n = 0

    u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])
            #out[i][k] = qs[i][k](r, κs_λ[i], κs_β[i]) # type stable

            n += 1
            out[n] = qs[i][k](r, κs_λ[i])
        end
    end

    #return out # type stable.
    return tuple(out...) # type unstable.
end

function assembleqsevals(z1::NTuple{N, Vector{Complex{T}}}, qs::Vector{Vector{Function}}) where {N,T}
    #
    out_qs = Vector{Vector{Vector{Complex{T}}}}(undef, length(qs))

    n = 0
    for i = 1:length(qs)
        out_qs[i] = Vector{Vector{Complex{T}}}(undef, length(qs[i]))

        for k = 1:length(qs[i])

            n += 1
            out_qs[i][k] = z1[n]
        end
    end

    return out_qs
end

#
function visualizesubgroups(P, qs_evals, q_singlets_evals, f::Function)

    for i = 1:length(qs_evals)
        for k = 1:length(qs_evals[i])

            PyPlot.plot(P, f.(qs_evals[i][k]), label = "sub-group ($(i),$(k))")
        end
    end

    for k = 1:length(q_singlets_evals)

        PyPlot.plot(P, f.(q_singlets_evals[k]), label = "singlet $(k)")
    end
end

"""
save resonance groupings of a compound. f is real, imag, or abs.
"""
function savefiggroups(save_path::String,
    title_string::String,
    P,
    q_U::Vector{Complex{T}},
    qs_U,
    q_singlets_U,
    f::Function) where T

    plot_obj = Plots.plot( P,
        f.(q_U),
        title = title_string,
        label = "sum of sub-models",
        seriestype = :line,
        ticks = :native,
        xlims = (P[1],P[end]),
        hover = P,
        linewidth = 4,
        xflip = true,
        size = (1600, 900))

    #
    for i = 1:length(qs_evals)
        for k = 1:length(qs_evals[i])

            Plots.plot!(plot_obj, P, f.(qs_evals[i][k]), label = "sub-group ($(i),$(k))",
                seriestype = :line,
                linestyle = :dot,
                xflip = true,
                linewidth = 4)
        end
    end

    for k = 1:length(q_singlets_evals)

        Plots.plot!(plot_obj, P, f.(q_singlets_evals[k]), label = "singlet $(k)",
            seriestype = :line,
            linestyle = :dot,
            xflip = true,
            linewidth = 4)
    end

    Plots.savefig(plot_obj, save_path)
end
