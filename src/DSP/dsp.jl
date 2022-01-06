
function getwraparoundDFTfreqs(N::Int, fs::T, ν_begin::T) where T

    U0 = getDFTfreqrange(N, fs)
    out, inds = wrapfreqrange(U0, ν_begin, fs)

    return U0, out, inds
end


# Assumes U0 is sorted in ascending order.
function wrapfreqrange(U0, ν_begin::T, fs::T) where T <: Real

    N = length(U0)
    ind = findfirst(xx->(xx>ν_begin), U0)

    # TODO handle these exceptions with more grace.
    @assert typeof(ind) == Int
    @assert ind <= N

    out = Vector{T}(undef, N)
    #out[1:ind] = U0[1:ind] .+ fs
    #out[ind+1:end] = U0[ind+1:end]

    M = N-ind
    out[1:M] = U0[ind+1:end]
    out[M+1:end] = U0[1:ind] .+ fs

    inds = collect(1:N)
    inds[1:M] = collect(ind+1:N)
    inds[M+1:end] = collect(1:ind)

    return out, inds
end

function gettimerange(N::Int, fs::T) where T
    Ts::T = 1/fs

    return zero(T):Ts:(N-1)*Ts
 end

 function getDFTfreqrange(N::Int, fs::T)::LinRange{T} where T
     a = zero(T)
     b = fs-fs/N

     return LinRange(a, b, N)
 end

 # case 1D.
function computeDTFTch3eq29(h, u::T, Λ)::Complex{T} where T <: Real

    # debug_array = zeros(Complex{T}, length(Λ))
    # debug_array2 = zeros(Complex{T}, length(Λ))
    # fs = 1/(Λ[2]-Λ[1])

    running_sum = zero(T)
    for i = 1:length(Λ)
        x = Λ[i]

        running_sum += h[i]*exp(-im*2*π*u*x)


        # # debug.
        # debug_array[i] = h[i]*exp(-im*2*π*u*x)
        #
        # #debug_array2[i] = h[i]*exp(-im*2*π*u*x +im*35*2*pi)
        # debug_array2[i] = h[i]*exp(-im*2*π*(u+3*fs)*x)
    end

    #println("debug: norm(debug_array -debug_array2 ) = ", norm(debug_array -debug_array2 ) )
    return running_sum
end
