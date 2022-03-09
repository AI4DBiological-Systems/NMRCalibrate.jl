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

    return out
end

function evalitpproxysysresonance(qs::Vector{Vector{Function}},
    u::T, d::Vector{T}, κs_λ::Vector{T},
    κs_β::Vector{Vector{T}})::Vector{Vector{Complex{T}}} where T <: Real

    @assert length(d) == length(qs)

    out = collect( zeros(Complex{T}, length(qs[i])) for i = 1:length(qs) )

    u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])
            out[i][k] = qs[i][k](r, κs_λ[i], κs_β[i])
        end
    end

    return out
end

# with proxy.
function evalitpproxycompoundresonance(u,
    A::NMRSpectraSimulator.CompoundFIDType{T}) where T <: Real

    u_rad = 2*π*u

    out_sys = evalitpproxysysresonance(A.qs, u, A.d, A.κs_λ, A.κs_β)

    out_singlets = evalsingletsresonance(u, A.d_singlets, A.αs_singlets, A.Ωs_singlets,
    A.β_singlets, A.λ0, A.κs_λ_singlets)

    return out_sys, out_singlets
end


function convertresonancetimeseries(g_U)
    @assert !isempty(g_U)
    M = length(g_U)

    out_qs = Vector{Vector{Vector{Float64}}}(undef, M)
    out_singlets = collect( Vector{Float64}(undef, M) for l = 1:N_singlets )

    for m = 1:M
        qs_eval, singlets_eval = g_U[m]

        for i = 1:length(qs_eval)
            for k = 1:length(qs_eval[i])
                out_qs[i][k][m] = qs_eval[i][k]
            end
        end

        for l = 1:length(singlets_eval)
            out_singlets[l][m] = singlets_eval[l]
        end
    end

    return out_qs, out_singlets
end
