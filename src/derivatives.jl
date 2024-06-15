
struct FirstDerivative{T} <: AbstractVector{T}
	x::AbstractVector{<:Real}
	y::AbstractVector{T}
	bc::NTuple{2}{T}
	direction::Symbol
	function FirstDerivative{T}(x, y, bc, direction) where {T}
		size(x) == size(y) || throw(DimensionMismatch(
			"cannot match grid of length $(length(x)) with vector of length $(length(y))"))
		direction ∈ (:upward, :downward) || throw(ArgumentError("direction must be :upward or :downward"))
		return new(x, y, bc, direction)
	end
end

function FirstDerivative(x::AbstractVector, y::AbstractVector; bc = (0, 0), direction = :upward)
	FirstDerivative{eltype(y)}(x, y, bc, direction)
end

Base.size(d::FirstDerivative) = (length(d.x), 1)

Base.IndexStyle(d::FirstDerivative) = IndexLinear()

function Base.getindex(d::FirstDerivative{T}, i::Int) where {T}
	x, y, bc, direction = d.x, d.y, d.bc, d.direction
	if direction == :upward
		if i == length(x)
			return convert(T, bc[end])
		else
			Δxp = x[min(i, length(x)-1)+1] - x[min(i, length(x)-1)]
			return convert(T, (y[i+1] - y[i]) / Δxp)
		end
	else
		if i == 1
			return convert(T, bc[1])
		else
			Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
			return convert(T, (y[i] - y[i-1]) / Δxm)
		end
	end
end


struct SecondDerivative{T} <: AbstractVector{T}
	x::AbstractVector{<:Real}
	y::AbstractVector{T}
	bc::NTuple{2}{T}
	function SecondDerivative{T}(x, y, bc) where {T}
		length(x) == length(y) || throw(DimensionMismatch(
			"cannot match grid of length $(length(x)) with vector of length $(length(y))"))
		return new(x, y, bc)
	end
end

function SecondDerivative(x::AbstractVector, y::AbstractVector; bc = (0, 0))
	SecondDerivative{eltype(y)}(x, y, bc)
end

Base.size(d::SecondDerivative) = (length(d.x), 1)

Base.IndexStyle(d::SecondDerivative) = IndexLinear()

function Base.getindex(d::SecondDerivative{T}, i::Int) where {T}
	x, y, bc = d.x, d.y, d.bc
    Δxp = x[min(i, length(x)-1)+1] - x[min(i, length(x)-1)]
    Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
    Δx = (Δxm + Δxp) / 2
    if i == 1 
    	return convert(T, y[2] / (Δxp * Δx) + (y[1] - bc[1] * Δxm) / (Δxm * Δx) - 2 * y[1] / (Δxp * Δxm))
    elseif i ==  length(x)
    	return convert(T, (y[end] + bc[end] * Δxp) / (Δxp * Δx) + y[end - 1] / (Δxm * Δx) - 2 * y[end] / (Δxp * Δxm))
    else
    	return convert(T, y[i + 1] / (Δxp * Δx) + y[i - 1] / (Δxm * Δx) - 2 * y[i] / (Δxp * Δxm))
    end
end



