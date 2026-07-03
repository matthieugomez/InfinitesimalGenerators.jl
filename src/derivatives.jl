function _fd_direction(direction::Symbol)
    direction ∈ (:forward, :up, :upward) && return :up
    direction ∈ (:backward, :down, :downward) && return :down
    throw(ArgumentError("direction must be :forward/:up or :backward/:down"))
end

"""
    FirstDerivative(x, f; direction = :forward, bc = (0, 0))
    FirstDerivative(grid, F, dim; direction = :forward, bc = (0, 0))

Lazily compute first-order derivatives (using finite-difference scheme) on a grid.
The result is an `AbstractArray` whose entries are computed on demand; use `collect`
to materialize it.

In one dimension, `x` is a grid vector and `f` is a vector with `length(f) == length(x)`.
In multiple dimensions, `grid` is a tuple of axis vectors such as `(xs, ys)`,
`F` is an array with size `(length(xs), length(ys))`, and `dim` is the dimension
number to differentiate. `NamedTuple` grids with symbolic dimensions are also
accepted as a convenience.

`direction` selects the one-sided difference: `:forward` (synonyms `:up`, `:upward`)
or `:backward` (synonyms `:down`, `:downward`).

`bc` is the value of the first derivative at the lower and upper boundaries.
"""
struct FirstDerivative{T, N, G <: Tuple, Y <: AbstractArray} <: AbstractArray{T, N}
    grid::G
    y::Y
    dim::Int
    bc::NTuple{2, T}
    direction::Symbol   # normalized to :up or :down
    function FirstDerivative(grid::Tuple, y::AbstractArray, dim::Integer; bc = (0, 0), direction = :forward)
        _fd_dimension(grid, y, dim)
        T = float(eltype(y))
        return new{T, ndims(y), typeof(grid), typeof(y)}(grid, y, Int(dim), NTuple{2, T}(bc), _fd_direction(direction))
    end
end

FirstDerivative(x::AbstractVector, y::AbstractVector; bc = (0, 0), direction = :forward) =
    FirstDerivative((x,), y, 1; bc = bc, direction = direction)

FirstDerivative(grid::NamedTuple, y::AbstractArray, dim::Integer; bc = (0, 0), direction = :forward) =
    FirstDerivative(Tuple(grid), y, dim; bc = bc, direction = direction)

FirstDerivative(grid::NamedTuple, y::AbstractArray, dim::Symbol; bc = (0, 0), direction = :forward) =
    FirstDerivative(Tuple(grid), y, _fd_dimension_index(grid, dim); bc = bc, direction = direction)

Base.size(d::FirstDerivative) = size(d.y)

Base.IndexStyle(::Type{<:FirstDerivative}) = IndexCartesian()

function Base.getindex(d::FirstDerivative{T, N}, I::Vararg{Int, N}) where {T, N}
    x, y = d.grid[d.dim], d.y
    i = I[d.dim]
    if d.direction === :up
        i == length(x) && return d.bc[2]
        Ip = ntuple(k -> k == d.dim ? I[k] + 1 : I[k], Val(N))
        return convert(T, (y[Ip...] - y[I...]) / (x[i + 1] - x[i]))
    else
        i == 1 && return d.bc[1]
        Im = ntuple(k -> k == d.dim ? I[k] - 1 : I[k], Val(N))
        return convert(T, (y[I...] - y[Im...]) / (x[i] - x[i - 1]))
    end
end

"""
    SecondDerivative(x, f; bc = (0, 0))
    SecondDerivative(grid, F, dim; bc = (0, 0))
    SecondDerivative(grid, F, dim1, dim2; direction = :up)

Lazily compute second-order derivatives (using finite-difference scheme) on a grid.
The result is an `AbstractArray` whose entries are computed on demand; use `collect`
to materialize it.

In one dimension, `x` is a grid vector and `f` is a vector with `length(f) == length(x)`.
In multiple dimensions, `grid` is a tuple of axis vectors such as `(xs, ys)`,
and `F` is an array with size `(length(xs), length(ys))`. `NamedTuple` grids
with symbolic dimensions are also accepted as a convenience.

Use `SecondDerivative(grid, F, 1, 1)` for own second derivatives and
`SecondDerivative(grid, F, 1, 2; direction = :up)` for directional cross
derivatives. The `:up` direction (synonyms `:forward`, `:upward`) is the
main-diagonal stencil; `:down` (synonyms `:backward`, `:downward`) is the
anti-diagonal stencil.

For own second derivatives, `bc` is the value of the *first* derivative at the
lower and upper boundaries. Cross derivatives do not accept boundary conditions
(the stencil is clamped at the edges of the grid) and throw if a nonzero `bc`
is passed.
"""
struct SecondDerivative{T, N, G <: Tuple, Y <: AbstractArray} <: AbstractArray{T, N}
    grid::G
    y::Y
    dim1::Int
    dim2::Int
    bc::NTuple{2, T}
    direction::Symbol   # normalized to :up or :down; only used when dim1 != dim2
    function SecondDerivative(grid::Tuple, y::AbstractArray, dim1::Integer, dim2::Integer; bc = (0, 0), direction = :up)
        _fd_dimension(grid, y, dim1)
        _fd_dimension(grid, y, dim2)
        if dim1 != dim2 && bc != (0, 0)
            throw(ArgumentError("cross derivatives do not accept boundary conditions; `bc` only applies to own second derivatives"))
        end
        T = float(eltype(y))
        return new{T, ndims(y), typeof(grid), typeof(y)}(grid, y, Int(dim1), Int(dim2), NTuple{2, T}(bc), _fd_direction(direction))
    end
end

SecondDerivative(x::AbstractVector, y::AbstractVector; bc = (0, 0)) =
    SecondDerivative((x,), y, 1, 1; bc = bc)

SecondDerivative(grid::Tuple, y::AbstractArray, dim::Integer; bc = (0, 0)) =
    SecondDerivative(grid, y, dim, dim; bc = bc)

SecondDerivative(grid::NamedTuple, y::AbstractArray, dim::Integer; bc = (0, 0)) =
    SecondDerivative(Tuple(grid), y, dim; bc = bc)

SecondDerivative(grid::NamedTuple, y::AbstractArray, dim::Symbol; bc = (0, 0)) =
    SecondDerivative(Tuple(grid), y, _fd_dimension_index(grid, dim); bc = bc)

SecondDerivative(grid::NamedTuple, y::AbstractArray, dim1::Integer, dim2::Integer; bc = (0, 0), direction = :up) =
    SecondDerivative(Tuple(grid), y, dim1, dim2; bc = bc, direction = direction)

SecondDerivative(grid::NamedTuple, y::AbstractArray, dim1::Symbol, dim2::Symbol; bc = (0, 0), direction = :up) =
    SecondDerivative(Tuple(grid), y, _fd_dimension_index(grid, dim1), _fd_dimension_index(grid, dim2); bc = bc, direction = direction)

Base.size(d::SecondDerivative) = size(d.y)

Base.IndexStyle(::Type{<:SecondDerivative}) = IndexCartesian()

function Base.getindex(d::SecondDerivative{T, N}, I::Vararg{Int, N}) where {T, N}
    y = d.y
    if d.dim1 == d.dim2
        x = d.grid[d.dim1]
        i = I[d.dim1]
        dxm, dxp, dx = _fd_grid_steps(x, i)
        if i == 1
            Ip = ntuple(k -> k == d.dim1 ? 2 : I[k], Val(N))
            return convert(T, y[Ip...] / (dxp * dx) + (y[I...] - d.bc[1] * dxm) / (dxm * dx) - 2 * y[I...] / (dxp * dxm))
        elseif i == length(x)
            Im = ntuple(k -> k == d.dim1 ? length(x) - 1 : I[k], Val(N))
            return convert(T, (y[I...] + d.bc[2] * dxp) / (dxp * dx) + y[Im...] / (dxm * dx) - 2 * y[I...] / (dxp * dxm))
        else
            Ip = ntuple(k -> k == d.dim1 ? I[k] + 1 : I[k], Val(N))
            Im = ntuple(k -> k == d.dim1 ? I[k] - 1 : I[k], Val(N))
            return convert(T, y[Ip...] / (dxp * dx) + y[Im...] / (dxm * dx) - 2 * y[I...] / (dxp * dxm))
        end
    else
        CI = CartesianIndex(I)
        shape = size(y)
        i1, i2 = I[d.dim1], I[d.dim2]
        x1, x2 = d.grid[d.dim1], d.grid[d.dim2]
        dx1m, dx1p, _ = _fd_grid_steps(x1, i1)
        dx2m, dx2p, _ = _fd_grid_steps(x2, i2)

        I00 = CI
        Ip0 = _fd_corner(CI, d.dim1, i1 + 1, d.dim2, i2, shape)
        Im0 = _fd_corner(CI, d.dim1, i1 - 1, d.dim2, i2, shape)
        I0p = _fd_corner(CI, d.dim1, i1, d.dim2, i2 + 1, shape)
        I0m = _fd_corner(CI, d.dim1, i1, d.dim2, i2 - 1, shape)
        Ipp = _fd_corner(CI, d.dim1, i1 + 1, d.dim2, i2 + 1, shape)
        Ipm = _fd_corner(CI, d.dim1, i1 + 1, d.dim2, i2 - 1, shape)
        Imp = _fd_corner(CI, d.dim1, i1 - 1, d.dim2, i2 + 1, shape)
        Imm = _fd_corner(CI, d.dim1, i1 - 1, d.dim2, i2 - 1, shape)

        if d.direction === :up
            return convert(T, (y[Ipp] - y[Ip0] - y[I0p] + y[I00]) / (2 * dx1p * dx2p) +
                (y[I00] - y[Im0] - y[I0m] + y[Imm]) / (2 * dx1m * dx2m))
        else
            return convert(T, -(y[Ipm] - y[Ip0] - y[I0m] + y[I00]) / (2 * dx1p * dx2m) -
                (y[I00] - y[Im0] - y[I0p] + y[Imp]) / (2 * dx1m * dx2p))
        end
    end
end

function _fd_dimension_index(grid::NamedTuple, dim::Symbol)
    haskey(grid, dim) || throw(ArgumentError("grid has no dimension `$dim`; valid dimensions are $(collect(keys(grid)))"))
    return findfirst(==(dim), keys(grid))
end

function _fd_dimension(grid::Tuple, y::AbstractArray, dim::Integer)
    1 <= dim <= length(grid) || throw(ArgumentError("dimension must be between 1 and $(length(grid)); got $dim"))
    for (d, x) in pairs(grid)
        x isa AbstractVector || throw(ArgumentError("grid axis $d must be a vector"))
        length(x) >= 2 || throw(ArgumentError("grid axis $d must contain at least two points"))
        all(x[i] < x[i + 1] for i in 1:(length(x) - 1)) ||
            throw(ArgumentError("grid axis $d must be strictly increasing"))
    end
    shape = ntuple(i -> length(grid[i]), length(grid))
    size(y) == shape || throw(DimensionMismatch("array has size $(size(y)) but the grid has size $shape"))
    d = Int(dim)
    x = grid[d]
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
