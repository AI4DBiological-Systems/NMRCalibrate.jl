
mutable struct FiniteSignalType{T}
    s0::Vector{T} # default signal.
    a::Float64
    b::Float64

    s::Vector{T} # buffer.
    #a::Float64
    #b::Float64
    a_ind
    b_ind
end

function FiniteSignalType(s, a, b)
    return FiniteSignalType(copy(s), a, b, copy(s), 1, length(s))
end

function resetsignal!(S::FiniteSignalType{T}) where T

    S.a_ind = 1
    S.b_ind = length(S.s0)
    S.s[:] = S.s0

    return nothing
end

# positive is delay, or shift towards more positive.
function applyintegershift!(S::FiniteSignalType{T}, d::Int) where T
    #
    S.a_ind = S.a_ind-d
    S.b_ind = S.b_ind-d

    return nothing
end

function evalfinitesignal(s::Vector{T}, i::Int)::T where T
    if 1 <= i <= length(s)
        return s[i]
    end

    return zero(T)
end

function evalsignal(S::FiniteSignalType{T}, n::Int)::T where T

    #ind_t_Float = convertcompactdomain(t, S.a, S.b, S.a_ind, S.b_ind)
    #ind_t_Float = round(Int, ind_t_Float)
    #@assert 1 <= n <= length(S.s)

    return evalfinitesignal(S.s, S.a_ind + n-1)
end

#  https://lcondat.github.io/publis/condat_icip08_rotations.pdf
function getthiriancoeffs!(b::Vector{T}, τ::T) where T <: Real

    #resize!(b, N)
    N = length(b)
    for k = 1:N
        b[k] = (-1)^k * binomial(N,k) * prod( (τ-n)/(τ-n-k) for n = 0:N)
    end

    return nothing
end

# https://lcondat.github.io/publis/condat_icip08_rotations.pdf
# τ ∈ [0, 0.5]
# mutates s and b_buffer.
function applydelaythirian!(b_buffer::Vector{T}, s::Vector{T}, τ::T) where T <: Real

    getthiriancoeffs!(b_buffer, τ)

    N = length(b_buffer)
    M = length(s)

    for i = M-1:-1:0
        for k = 1:N

            #s[i+1] += b_buffer[k]*(s[i-k+1]-s[i+k+1])

            tmp = evalfinitesignal(s, i-k+1) - evalfinitesignal(s, i+k+1)

            s[i+1] += b_buffer[k]*tmp
        end
    end

    return nothing
end

function applyadvancethirian!(b_buffer::Vector{T}, s::Vector{T}, τ::T) where T <: Real

    getthiriancoeffs!(b_buffer, τ)

    N = length(b_buffer)
    M = length(s)

    #reverse!(s)
    for i = M-1:-1:0
        for k = 1:N

            # tmp = evalfinitesignal(s, i-k+1) - evalfinitesignal(s, i+k+1)
            # s[i+1] += b_buffer[k]*tmp


            # reverse: if we have the ORDERED lists:
            # U = {1, 2, ..., M}, V = {M, M-1, ... , 1}, then
            # j = M-i+1, where i is a given element from U, and j is the
            #   corresponding element in V.

            j = M-i #j = M-(i+1)+1
            j1 = j+k #j1 = M-(i-k+1)+1
            j2 = j-k #j2 = M-(i+k+1)+1
            tmp = evalfinitesignal(s, j1) - evalfinitesignal(s, j2)
            s[j] += b_buffer[k]*tmp
        end
    end
    #reverse!(s)

    return nothing
end

function applythirian!(b_buffer)

    if τ > 0.5
        return nothing # shift, then apply thirian.
    end

    return nothing # applythirian!()
end


function applyshift!(S::FiniteSignalType{T}, shift::T, b_thirian) where T

    d = round(Int, shift)
    τ = shift-d
    println(τ)

    applyintegershift!(S, d)
    if τ > 0
        applydelaythirian!(b_thirian, S.s, τ)
    else
        applyadvancethirian!(b_thirian, S.s, abs(τ))
    end

    return nothing
end


# other interesting articles about fractional delay: https://www.dsprelated.com/showarticle/1327.php
