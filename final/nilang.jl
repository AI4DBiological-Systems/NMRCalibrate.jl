using NiLang, BenchmarkTools

using LinearAlgebra, FFTW
import BSON, Statistics, PyPlot, Random


# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

import MultistartOptimization, NLopt
import MonotoneMaps
import Interpolations

PyPlot.close("all")
fig_num = 1

Random.seed!(25)
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])



#@i f(a, b) = a += b
#f(3,6)
#@btime (~f)(-13.23,9.34);

@i function r_axpy!(a::T, x::AbstractVector{T}, y!::AbstractVector{T}) where T
    @safe @assert length(x) == length(y!)
    for i=1:length(x)
        y![i] += a * x[i]
    end
end

@i function r_loss(out!, a, x, y!, z)
    r_axpy!(a, x, y!)
    for i=1:length(z)
    	out! += z[i] * y![i]
    end
end


out, a, x, y, z = 0.0, 2.0, randn(3), randn(3), randn(3)

# out, a, x, y, z = r_loss(out, a, x, y, z)
# out, a, x, y, z = (~r_loss)(out, a, x, y, z)
@instr r_loss(out, a, x, y, z)
@instr (~r_loss)(out, a, x, y, z)



######### box-muller.
@i function boxmuller(x::T, y::T) where T
    @routine @invcheckoff begin
        @zeros T θ logx _2logx
        θ += 2π * y
        logx += log(x)
        _2logx += -2 * logx
    end

    # store results
    z1 ← zero(T)
    z2 ← zero(T)
    z1 += _2logx ^ 0.5
    ROT(z1, z2, θ) # https://docs.juliahub.com/NiLang/U2llx/0.8.1/instructions/
    ~@routine

    SWAP(x, z1)
    SWAP(y, z2)

    # arithmetic uncomputing: recomputing the original values of `x` and `y` to deallocate z1 and z2
    @routine @invcheckoff begin
        @zeros T at sq _halfsq
        at += atan(y, x)
        if (y < 0, ~)
            at += T(2π)
        end
        sq += x ^ 2
        sq += y ^ 2
        _halfsq -= sq / 2
    end
    z1 -= exp(_halfsq)
    z2 -= at / (2π)
    @invcheckoff z1 → zero(T)
    @invcheckoff z2 → zero(T)
    ~@routine
end

# I am here. do a algebraic sigmoid version, and see if nilang can reverse that.



@i function boxmuller(x::T, y::T) where T
    @routine @invcheckoff begin
        @zeros T θ logx _2logx
        θ += 2π * y
        logx += log(x)
        _2logx += -2 * logx
    end

    # store results
    z1 ← zero(T)
    z2 ← zero(T)
    z1 += _2logx ^ 0.5
    ROT(z1, z2, θ) # https://docs.juliahub.com/NiLang/U2llx/0.8.1/instructions/
    ~@routine
end
