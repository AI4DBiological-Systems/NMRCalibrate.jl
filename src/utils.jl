function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end