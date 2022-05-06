using LinearAlgebra
using FFTW
import PyPlot
import BSON
import JLD

#import Clustering
import Statistics
import Random
import TotalLeastSquares

PyPlot.close("all")
fig_num = 1


#f = (xx,uu)->(sinc(uu*xx)+cos(xx))

a = 2.3
f = (xx,uu)->(a/(a+(xx-uu)^2))

L = 20
λs = rand(L) .* 20
c = rand(L)

function getA(U,f,λs)
    M = length(U)
    L = length(λs)

    A = Matrix{Float64}(undef, M, L)
    for m = 1:M
        for l = 1:L
            A[m,l] = f(U[m], λs[l])
        end
    end

    return A
end

function evaly(u, f, λs, c::Vector{T})::T where T<: Real
    @assert length(c) == length(λs)
    return sum(c[i]*f(u,λs[i]) for i = 1:length(c))
end

### set up least squares.
U = LinRange(-10,10, 500)
Y = collect( evaly(U[m], f, λs, c) for m = 1:length(U))

A = getA(U,f,λs)
norm(A*c - Y)

### set up total least squares.
@time c1,c2,c3 = TotalLeastSquares.flts(A, Y, return_set = true)
