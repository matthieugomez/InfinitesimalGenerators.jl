"""
    MultivariateDiffusionProcess(grid; drift, variance, covariance = nothing)

Returns a diffusion process on a tensor-product grid.

`grid` is a `NamedTuple` mapping state names to strictly increasing grid vectors.
`drift` has the same names as `grid`; each entry is either a scalar or an array
with the same size as the tensor grid. `variance` has the same names and gives
the diagonal entries of the instantaneous covariance matrix. Cross terms are
optional and are passed as pair names such as `(; xy = axy)`.
"""
struct MultivariateDiffusionProcess{N, G, D, V, C} <: ContinuousTimeMarkovProcess{N}
    grid::G
    drift::D
    variance::V
    covariance::C

    function MultivariateDiffusionProcess(grid::NamedTuple, drift, variance, covariance)
        # Validate grid.
        for (name, x) in pairs(grid)
            x isa AbstractVector || throw(ArgumentError("grid entry `$name` must be a vector"))
            length(x) >= 2 || throw(ArgumentError("grid entry `$name` must contain at least two points"))
            all(x[i] < x[i + 1] for i in 1:(length(x) - 1)) ||
                throw(ArgumentError("grid entry `$name` must be strictly increasing"))
        end

        shape = ntuple(i -> length(grid[i]), length(grid))
        names = keys(grid)

        drift_arrays = _mvd_coefficient_arrays(drift, names, shape, :drift)
        variance_arrays = _mvd_coefficient_arrays(variance, names, shape, :variance)
        for (name, a) in pairs(variance_arrays)
            any(<(0), a) && throw(ArgumentError("`variance.$name` must be nonnegative"))
        end

        # Collect off-diagonal covariance entries.
        covpairs = _mvd_covariance_pairs(names)
        covnames = Tuple(map(last, covpairs))
        if covariance === nothing
            covariance_arrays = NamedTuple{covnames}(ntuple(_ -> zeros(Float64, shape), length(covnames)))
        else
            covariance isa NamedTuple ||
                throw(ArgumentError("`covariance` must be a NamedTuple with pair keys such as `(; xy = covxy)`"))
            invalid = setdiff(keys(covariance), covnames)
            isempty(invalid) ||
                throw(ArgumentError("unknown covariance entries $(collect(invalid)); valid entries are $(collect(covnames))"))
            covariance_arrays = NamedTuple{covnames}(map(covnames) do key
                haskey(covariance, key) ? _mvd_coefficient_array(covariance[key], shape, Symbol(:covariance_, key)) : zeros(Float64, shape)
            end)
        end

        if length(names) > 1
            for I in CartesianIndices(shape)
                Σ = zeros(Float64, length(names), length(names))
                for d in eachindex(names)
                    Σ[d, d] = variance_arrays[names[d]][I]
                end
                for (d1, d2, key) in covpairs
                    Σ[d1, d2] = covariance_arrays[key][I]
                    Σ[d2, d1] = covariance_arrays[key][I]
                end
                scale = max(1.0, maximum(abs, Σ))
                eigmin(Symmetric(Σ)) >= -sqrt(eps(Float64)) * scale ||
                    throw(ArgumentError("instantaneous covariance matrix must be positive semidefinite at grid index $I"))
            end
        end

        return new{length(grid), typeof(grid), typeof(drift_arrays), typeof(variance_arrays), typeof(covariance_arrays)}(
            grid, drift_arrays, variance_arrays, covariance_arrays)
    end
end

function MultivariateDiffusionProcess(grid; drift, variance, covariance = nothing)
    throw(ArgumentError("`grid` must be a NamedTuple, e.g. `(; x = xs, y = ys)`"))
end

function MultivariateDiffusionProcess(grid::NamedTuple; drift, variance, covariance = nothing)
    MultivariateDiffusionProcess(grid, drift, variance, covariance)
end

state_space(X::MultivariateDiffusionProcess) = Tuple(X.grid)

function _mvd_coefficient_arrays(x, names, shape, label::Symbol)
    throw(ArgumentError("`$label` must be a NamedTuple with the same names as the grid"))
end

function _mvd_coefficient_arrays(x::NamedTuple, names, shape, label::Symbol)
    Set(keys(x)) == Set(names) ||
        throw(ArgumentError("`$label` must have the same names as the grid: got $(collect(keys(x))), expected $(collect(names))"))
    NamedTuple{names}(map(name -> _mvd_coefficient_array(x[name], shape, Symbol(label, :_, name)), names))
end

_mvd_coefficient_array(x::Number, shape, label::Symbol) = fill(float(x), shape)

function _mvd_coefficient_array(x, shape, label::Symbol)
    size(x) == shape || throw(ArgumentError("coefficient `$label` has size $(size(x)) but the state grid has size $shape"))
    collect(float.(x))
end

function _mvd_covariance_pairs(names)
    pairs = Tuple((i, j, Symbol(names[i], names[j])) for i in 1:length(names) for j in (i + 1):length(names))
    pairnames = map(last, pairs)
    if length(unique(pairnames)) != length(pairnames)
        throw(ArgumentError("ambiguous covariance pair names from state names $(collect(names)); rename states"))
    end
    return pairs
end

function _mvd_add_weight!(rows, cols, vals, row::Int, col::Int, weight)
    iszero(weight) && return nothing
    push!(rows, row)
    push!(cols, col)
    push!(vals, weight)
    return nothing
end

function _mvd_add_first_derivative!(rows, cols, vals, L, row, I, d, x, coefficient)
    i = I[d]
    dxm, dxp, _ = _fd_grid_steps(x, i)
    if coefficient >= 0
        i < length(x) || return nothing
        Ip = _fd_replace_index(I, d, i + 1)
        _mvd_add_weight!(rows, cols, vals, row, L[Ip], coefficient / dxp)
        _mvd_add_weight!(rows, cols, vals, row, row, -coefficient / dxp)
    else
        i > 1 || return nothing
        Im = _fd_replace_index(I, d, i - 1)
        _mvd_add_weight!(rows, cols, vals, row, row, coefficient / dxm)
        _mvd_add_weight!(rows, cols, vals, row, L[Im], -coefficient / dxm)
    end
    return nothing
end

function _mvd_add_second_derivative!(rows, cols, vals, L, row, I, d, x, coefficient)
    iszero(coefficient) && return nothing
    i = I[d]
    dxm, dxp, dx = _fd_grid_steps(x, i)
    if 1 < i < length(x)
        Ip = _fd_replace_index(I, d, i + 1)
        Im = _fd_replace_index(I, d, i - 1)
        _mvd_add_weight!(rows, cols, vals, row, L[Ip], coefficient / (dxp * dx))
        _mvd_add_weight!(rows, cols, vals, row, L[Im], coefficient / (dxm * dx))
        _mvd_add_weight!(rows, cols, vals, row, row, -2 * coefficient / (dxp * dxm))
    elseif i == 1
        Ip = _fd_replace_index(I, d, 2)
        _mvd_add_weight!(rows, cols, vals, row, L[Ip], coefficient / (dxp * dx))
        _mvd_add_weight!(rows, cols, vals, row, row, coefficient / (dxm * dx) - 2 * coefficient / (dxp * dxm))
    else
        Im = _fd_replace_index(I, d, length(x) - 1)
        _mvd_add_weight!(rows, cols, vals, row, row, coefficient / (dxp * dx) - 2 * coefficient / (dxp * dxm))
        _mvd_add_weight!(rows, cols, vals, row, L[Im], coefficient / (dxm * dx))
    end
    return nothing
end

function _mvd_add_cross_derivative!(rows, cols, vals, L, row, I, d1, d2, x1, x2, shape, coefficient)
    iszero(coefficient) && return nothing
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

    if coefficient >= 0
        w1 = coefficient / (2 * dx1p * dx2p)
        _mvd_add_weight!(rows, cols, vals, row, L[Ipp], w1)
        _mvd_add_weight!(rows, cols, vals, row, L[Ip0], -w1)
        _mvd_add_weight!(rows, cols, vals, row, L[I0p], -w1)
        _mvd_add_weight!(rows, cols, vals, row, L[I00], w1)

        w2 = coefficient / (2 * dx1m * dx2m)
        _mvd_add_weight!(rows, cols, vals, row, L[I00], w2)
        _mvd_add_weight!(rows, cols, vals, row, L[Im0], -w2)
        _mvd_add_weight!(rows, cols, vals, row, L[I0m], -w2)
        _mvd_add_weight!(rows, cols, vals, row, L[Imm], w2)
    else
        w1 = -coefficient / (2 * dx1p * dx2m)
        _mvd_add_weight!(rows, cols, vals, row, L[Ipm], w1)
        _mvd_add_weight!(rows, cols, vals, row, L[Ip0], -w1)
        _mvd_add_weight!(rows, cols, vals, row, L[I0m], -w1)
        _mvd_add_weight!(rows, cols, vals, row, L[I00], w1)

        w2 = -coefficient / (2 * dx1m * dx2p)
        _mvd_add_weight!(rows, cols, vals, row, L[I00], w2)
        _mvd_add_weight!(rows, cols, vals, row, L[Im0], -w2)
        _mvd_add_weight!(rows, cols, vals, row, L[I0p], -w2)
        _mvd_add_weight!(rows, cols, vals, row, L[Imp], w2)
    end
    return nothing
end

"""
    generator(X::MultivariateDiffusionProcess; check = :throw, check_tol = 1e-12)

Assemble the sparse infinitesimal generator matrix. The matrix acts on `vec(f)`,
where `f` is an array on the tensor grid described by `state_space(X)`.

By default, throws an error if the finite-difference stencil creates a negative
off-diagonal entry. Use `check = :warn` to return the matrix with a warning, or
`check = false` to skip the check.
"""
function generator(X::MultivariateDiffusionProcess; check = :throw, check_tol = 1e-12)
    grid = X.grid
    names = keys(grid)
    N = length(names)
    shape = size(X)
    cartesian_indices = CartesianIndices(shape)
    linear_indices = LinearIndices(shape)
    n = length(cartesian_indices)

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for I in cartesian_indices
        row = linear_indices[I]
        for d in 1:N
            name = names[d]
            _mvd_add_first_derivative!(rows, cols, vals, linear_indices, row, I, d, grid[name], X.drift[name][I])
            _mvd_add_second_derivative!(rows, cols, vals, linear_indices, row, I, d, grid[name], 0.5 * X.variance[name][I])
        end
        for (d1, d2, key) in _mvd_covariance_pairs(names)
            _mvd_add_cross_derivative!(rows, cols, vals, linear_indices, row, I, d1, d2, grid[names[d1]], grid[names[d2]], shape, X.covariance[key][I])
        end
    end

    A = sparse(rows, cols, vals, n, n)
    rowsum = vec(sum(A, dims = 2))
    A = A - sparse(1:n, 1:n, rowsum, n, n)

    check_mode = check === true ? :throw : check
    if check_mode !== false
        check_mode ∈ (:throw, :warn) ||
            throw(ArgumentError("`check` must be `:throw`, `:warn`, `true`, or `false`"))
        check_rows, check_cols, check_vals = findnz(A)
        for k in eachindex(check_vals)
            if check_rows[k] != check_cols[k] && check_vals[k] < -check_tol
                rowI = cartesian_indices[check_rows[k]]
                colI = cartesian_indices[check_cols[k]]
                message = string(
                    "MultivariateDiffusionProcess generator has a negative off-diagonal entry ",
                    "(row index $rowI, column index $colI, value $(check_vals[k])) and is not a valid ",
                    "Markov generator. Rows still sum to zero, but a negative off-diagonal entry is a ",
                    "negative transition rate. This usually happens with correlated states when the ",
                    "instantaneous covariance matrix is positive semidefinite but the grid-scaled ",
                    "diagonal-dominance condition for the directional cross-derivative stencil fails. ",
                    "For a two-state pair (x, y), the local condition is roughly ",
                    "abs(covariance.xy) <= variance.x * Δy / Δx and ",
                    "abs(covariance.xy) <= variance.y * Δx / Δy; on equally spaced grids this reduces ",
                    "to abs(covariance.xy) <= min(variance.x, variance.y). To fix it, refine or rescale ",
                    "the grids so Δx / Δy is closer to sqrt(variance.x / variance.y), use a state ",
                    "transformation that balances local volatilities, or reduce the covariance. Use ",
                    "`generator(X; check = :warn)` or `check = false` only if you intentionally want ",
                    "the raw finite-difference operator without a Markov-process interpretation.")
                if check_mode === :warn
                    @warn message row = check_rows[k] col = check_cols[k] value = check_vals[k]
                else
                    throw(ArgumentError(message))
                end
                break
            end
        end
    end
    return A
end
