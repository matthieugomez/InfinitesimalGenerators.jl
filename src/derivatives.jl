"""
    FirstDerivative(x, f; direction = :forward, bc = (0, 0))
    FirstDerivative(grid, F, dim; direction = :forward, bc = (0, 0))

Lazily compute first-order derivatives (using finite-difference scheme) on a grid.

In one dimension, `x` is a grid vector and `f` is a vector with `length(f) == length(x)`.
In multiple dimensions, `grid` is a `NamedTuple` such as `(; x = xs, y = ys)`,
`F` is an array with size `(length(xs), length(ys))`, and `dim` is the dimension
name to differentiate.

`bc` is the value of the first derivative at the lower and upper boundaries.
"""
struct FirstDerivative{T, X <: AbstractVector{<:Real}, Y <: AbstractVector{<: Real}} <: AbstractVector{T}
	x::X
	y::Y
	bc::NTuple{2, T}
	direction::Symbol
	function FirstDerivative(x, y, bc, direction)
		size(x) == size(y) || throw(DimensionMismatch(
			"cannot match grid of length $(length(x)) with vector of length $(length(y))"))
		direction ∈ (:forward, :backward, :upward, :downward) || throw(ArgumentError("direction must be :forward/:upward or :backward/:downward"))
		return new{float(eltype(y)), typeof(x), typeof(y)}(x, y, bc, direction)
	end
end

FirstDerivative(x, y; bc = (0, 0), direction = :forward) = FirstDerivative(x, y, bc, direction)

Base.size(d::FirstDerivative) = (length(d.x),)

Base.IndexStyle(d::FirstDerivative) = IndexLinear()

function Base.getindex(d::FirstDerivative{T}, i::Int) where {T}
	x, y, bc, direction = d.x, d.y, d.bc, d.direction
	if direction ∈ (:forward, :upward)
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


"""
    SecondDerivative(x, f; bc = (0, 0))
    SecondDerivative(grid, F, dim; bc = (0, 0))
    SecondDerivative(grid, F, dim1, dim2; direction = :up, bc = (0, 0))

Lazily compute second-order derivatives (using finite-difference scheme) on a grid

In one dimension, `x` is a grid vector and `f` is a vector with `length(f) == length(x)`.
In multiple dimensions, `grid` is a `NamedTuple` such as `(; x = xs, y = ys)`,
and `F` is an array with size `(length(xs), length(ys))`.

Use `SecondDerivative(grid, F, :x, :x)` for own second derivatives and
`SecondDerivative(grid, F, :x, :y; direction = :up)` for directional cross
derivatives. The `:up` direction is the main-diagonal stencil; `:down` is the
anti-diagonal stencil.
"""
struct SecondDerivative{T, X <: AbstractVector{<:Real}, Y <: AbstractVector{<: Real}} <: AbstractVector{T}
	x::X
	y::Y
	bc::NTuple{2, T}
	function SecondDerivative(x, y, bc)
		length(x) == length(y) || throw(DimensionMismatch(
			"cannot match grid of length $(length(x)) with vector of length $(length(y))"))
		return new{float(eltype(y)), typeof(x), typeof(y)}(x, y, bc)
	end
end

SecondDerivative(x, y; bc = (0, 0)) = SecondDerivative(x, y, bc)

Base.size(d::SecondDerivative) = (length(d.x),)

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


function FirstDerivative(grid::NamedTuple, y::AbstractArray, dim::Symbol; bc = (0, 0), direction = :forward)
    direction ∈ (:forward, :backward, :upward, :downward) || throw(ArgumentError("direction must be :forward/:upward or :backward/:downward"))
    x, d, shape = _fd_dimension(grid, y, dim)
    dy = Array{float(eltype(y))}(undef, shape)

    for I in CartesianIndices(shape)
        i = I[d]
        if direction ∈ (:forward, :upward)
            if i == length(x)
                dy[I] = bc[end]
            else
                Ip = _fd_replace_index(I, d, i + 1)
                dxp = x[i + 1] - x[i]
                dy[I] = (y[Ip] - y[I]) / dxp
            end
        else
            if i == 1
                dy[I] = bc[1]
            else
                Im = _fd_replace_index(I, d, i - 1)
                dxm = x[i] - x[i - 1]
                dy[I] = (y[I] - y[Im]) / dxm
            end
        end
    end
    return dy
end

SecondDerivative(grid::NamedTuple, y::AbstractArray, dim::Symbol; bc = (0, 0)) =
    SecondDerivative(grid, y, dim, dim; bc = bc)

function SecondDerivative(grid::NamedTuple, y::AbstractArray, dim1::Symbol, dim2::Symbol; bc = (0, 0), direction = :up)
    if dim1 != dim2
        direction ∈ (:up, :down, :upward, :downward) || throw(ArgumentError("direction must be :up/:upward or :down/:downward"))
        x1, d1, shape = _fd_dimension(grid, y, dim1)
        x2, d2, _ = _fd_dimension(grid, y, dim2)
        d12y = Array{float(eltype(y))}(undef, shape)

        for I in CartesianIndices(shape)
            i1, i2 = I[d1], I[d2]
            dx1m, dx1p, _ = _fd_grid_steps(x1, i1)
            dx2m, dx2p, _ = _fd_grid_steps(x2, i2)

            I00 = I
            Ip0 = _fd_corner(I, d1, i1 + 1, d2, i2, shape)
            Im0 = _fd_corner(I, d1, i1 - 1, d2, i2, shape)
            I0p = _fd_corner(I, d1, i1, d2, i2 + 1, shape)
            I0m = _fd_corner(I, d1, i1, d2, i2 - 1, shape)
            Ipp = _fd_corner(I, d1, i1 + 1, d2, i2 + 1, shape)
            Ipm = _fd_corner(I, d1, i1 + 1, d2, i2 - 1, shape)
            Imp = _fd_corner(I, d1, i1 - 1, d2, i2 + 1, shape)
            Imm = _fd_corner(I, d1, i1 - 1, d2, i2 - 1, shape)

            if direction ∈ (:up, :upward)
                d12y[I] = (y[Ipp] - y[Ip0] - y[I0p] + y[I00]) / (2 * dx1p * dx2p) +
                    (y[I00] - y[Im0] - y[I0m] + y[Imm]) / (2 * dx1m * dx2m)
            else
                d12y[I] = -(y[Ipm] - y[Ip0] - y[I0m] + y[I00]) / (2 * dx1p * dx2m) -
                    (y[I00] - y[Im0] - y[I0p] + y[Imp]) / (2 * dx1m * dx2p)
            end
        end
        return d12y
    end

    dim = dim1
    x, d, shape = _fd_dimension(grid, y, dim)
    d2y = Array{float(eltype(y))}(undef, shape)

    for I in CartesianIndices(shape)
        i = I[d]
        dxm, dxp, dx = _fd_grid_steps(x, i)
        if i == 1
            Ip = _fd_replace_index(I, d, 2)
            d2y[I] = y[Ip] / (dxp * dx) + (y[I] - bc[1] * dxm) / (dxm * dx) - 2 * y[I] / (dxp * dxm)
        elseif i == length(x)
            Im = _fd_replace_index(I, d, length(x) - 1)
            d2y[I] = (y[I] + bc[end] * dxp) / (dxp * dx) + y[Im] / (dxm * dx) - 2 * y[I] / (dxp * dxm)
        else
            Ip = _fd_replace_index(I, d, i + 1)
            Im = _fd_replace_index(I, d, i - 1)
            d2y[I] = y[Ip] / (dxp * dx) + y[Im] / (dxm * dx) - 2 * y[I] / (dxp * dxm)
        end
    end
    return d2y
end

function _fd_dimension(grid::NamedTuple, y::AbstractArray, dim::Symbol)
    haskey(grid, dim) || throw(ArgumentError("grid has no dimension `$dim`; valid dimensions are $(collect(keys(grid)))"))
    for (name, x) in pairs(grid)
        x isa AbstractVector || throw(ArgumentError("grid entry `$name` must be a vector"))
        length(x) >= 2 || throw(ArgumentError("grid entry `$name` must contain at least two points"))
        all(x[i] < x[i + 1] for i in 1:(length(x) - 1)) ||
            throw(ArgumentError("grid entry `$name` must be strictly increasing"))
    end
    shape = ntuple(i -> length(grid[i]), length(grid))
    size(y) == shape || throw(DimensionMismatch("array has size $(size(y)) but the grid has size $shape"))
    d = findfirst(==(dim), keys(grid))
    x = grid[dim]
    return x, d, shape
end

_fd_replace_index(I::CartesianIndex, d::Int, value::Int) =
    CartesianIndex(ntuple(k -> k == d ? value : I[k], length(I)))

function _fd_grid_steps(x::AbstractVector, i::Int)
    dxp = x[min(i, length(x) - 1) + 1] - x[min(i, length(x) - 1)]
    dxm = x[max(i - 1, 1) + 1] - x[max(i - 1, 1)]
    dx = (dxm + dxp) / 2
    return dxm, dxp, dx
end

function _fd_corner(I::CartesianIndex, d1::Int, i1::Int, d2::Int, i2::Int, shape)
    CartesianIndex(ntuple(k -> k == d1 ? clamp(i1, 1, shape[k]) : (k == d2 ? clamp(i2, 1, shape[k]) : I[k]), length(I)))
end
