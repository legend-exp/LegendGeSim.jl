_fast_ustrip(x::Number) = x
_fast_ustrip(A::AbstractArray{T}) where {T<:Real} = A
_fast_ustrip(A::AbstractRange{<:Quantity{T}}) where {T<:Number} = ustrip.(A)
_fast_ustrip(A::AbstractArray{<:Quantity{T}}) where {T<:Number} = reinterpret(T, A)
