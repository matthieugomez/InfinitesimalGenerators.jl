# bc corresponds to value of first derivative outside the grid. Default to (0, 0), which corresponds to reflecting.
struct FirstDerivative{T, X <: AbstractVector{<:Real}, Y <: AbstractVector{<: Real}} <: AbstractVector{T}
	x::X
	y::Y
	bc::NTuple{2}{T}
	direction::Symbol
	function FirstDerivative(x, y, bc, direction)
		size(x) == size(y) || throw(DimensionMismatch(
			"cannot match grid of length $(length(x)) with vector of length $(length(y))"))
		direction ∈ (:forward, :backward) || throw(ArgumentError("direction must be :forward or :backward"))
		return new{float(eltype(y)), typeof(x), typeof(y)}(x, y, bc, direction)
	end
end

FirstDerivative(x, y; bc = (0, 0), direction = :forward) = FirstDerivative(x, y, bc, direction)

Base.size(d::FirstDerivative) = (length(d.x), 1)

Base.IndexStyle(d::FirstDerivative) = IndexLinear()

function Base.getindex(d::FirstDerivative{T}, i::Int) where {T}
	x, y, bc, direction = d.x, d.y, d.bc, d.direction
	if direction == :forward
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


struct SecondDerivative{T, X <: AbstractVector{<:Real}, Y <: AbstractVector{<: Real}} <: AbstractVector{T}
	x::X
	y::Y
	bc::NTuple{2}{T}
	function SecondDerivative(x, y, bc)
		length(x) == length(y) || throw(DimensionMismatch(
			"cannot match grid of length $(length(x)) with vector of length $(length(y))"))
		return new{float(eltype(y)), typeof(x), typeof(y)}(x, y, bc)
	end
end

SecondDerivative(x, y; bc = (0, 0)) = SecondDerivative(x, y, bc)

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



