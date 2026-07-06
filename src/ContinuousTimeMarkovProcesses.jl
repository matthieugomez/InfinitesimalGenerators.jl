"""
    ContinuousTimeMarkovProcess{N}

Abstract supertype for continuous-time Markov processes with `N` tensor-product
state-space axes.

Concrete subtypes must define `state_space(X)` and `generator(X)`.
`state_space(X)` must return a tuple of state-space axes, including for
one-dimensional processes. `size(X)` is derived from those axes, and `ndims(X)`
returns the type parameter `N`.
"""
abstract type ContinuousTimeMarkovProcess{N} end
# This type should define generator(), which returns a generator matrix 𝔸 such that
# 𝔸f = lim_{t→0} (E[f(x_t)|x_0=x] - f(x))/t

"""
    generator(X::ContinuousTimeMarkovProcess)

Return the generator (transition-rate) matrix `𝔸` of the discretized process — the matrix
representation of the infinitesimal generator

    𝔸: f ↦ lim_{t→0} (E[f(x_t) | x_0 = x] - f(x)) / t

acting on the state space flattened in column-major order. Rows sum to zero and off-diagonal
entries are non-negative transition rates.
"""
function generator end

"""
    state_space(X::ContinuousTimeMarkovProcess)

Return the state-space axes of `X` as a tuple.

The return value is always a tuple, including for one-dimensional processes:
`state_space(X)[1]` is the first axis, and `only(state_space(X))` is the grid of
a one-dimensional process.
"""
function state_space end

"""
    size(X::ContinuousTimeMarkovProcess)

Return the tensor shape of arrays defined on the state space of `X`.

This is the length of each state-space axis. Operators such as
`stationary_distribution(X)` and `feynman_kac(X, ...)` use this shape for
process-level inputs and outputs.
"""
Base.size(X::ContinuousTimeMarkovProcess) = map(length, state_space(X))

"""
    length(X::ContinuousTimeMarkovProcess)

Return the total number of grid points or finite states in `X`, equal to
`prod(size(X))`. This is the length of flattened state vectors and the row/column
dimension of `generator(X)`.
"""
Base.length(X::ContinuousTimeMarkovProcess) = prod(size(X))

"""
    ndims(X::ContinuousTimeMarkovProcess)

Return the number of tensor-product state-space axes of `X`.
"""
Base.ndims(::ContinuousTimeMarkovProcess{N}) where {N} = N

function _flatten_state_vector_argument(X::ContinuousTimeMarkovProcess, x, label::Symbol)
    length(x) == length(X) ||
        throw(DimensionMismatch("`$label` has length $(length(x)) but process has length $(length(X))"))
    return vec(x)
end

function _flatten_state_time_argument(X::ContinuousTimeMarkovProcess, x, ts, label::Symbol)
    state_shape = size(X)
    if x isa AbstractVector
        length(x) == length(X) ||
            throw(DimensionMismatch("`$label` has length $(length(x)) but process has length $(length(X))"))
        return x
    elseif size(x) == state_shape
        return vec(x)
    elseif x isa AbstractMatrix && size(x, 1) == length(X) && size(x, 2) in (1, length(ts))
        return x
    elseif length(size(x)) == length(state_shape) + 1 &&
            size(x)[1:end - 1] == state_shape &&
            size(x, ndims(x)) in (1, length(ts))
        return reshape(x, length(X), size(x, ndims(x)))
    else
        throw(DimensionMismatch(
            "`$label` must have shape $(state_shape), length $(length(X)), " *
            "or time-varying shape ($(state_shape)..., $(length(ts)))"))
    end
end

function _reshape_state_output(X::ContinuousTimeMarkovProcess, x::AbstractVector)
    length(x) == length(X) ||
        throw(DimensionMismatch("state vector has length $(length(x)) but process has length $(length(X))"))
    length(size(X)) == 1 && return x
    return reshape(x, size(X))
end

function _reshape_state_time_output(X::ContinuousTimeMarkovProcess, x::AbstractMatrix)
    size(x, 1) == length(X) ||
        throw(DimensionMismatch("state-time array has $(size(x, 1)) rows but process has length $(length(X))"))
    length(size(X)) == 1 && return x
    return reshape(x, size(X)..., size(x, 2))
end

"""
    stationary_distribution(X::ContinuousTimeMarkovProcess; δ = 0.0, rebirth = Ones(length(X)), kwargs...)

Compute the stationary distribution of the Markov process `X`, returned as an array of
probability *masses* shaped like `size(X)`: the array sums to one with no grid
weights, and expectations are unweighted dot products `sum(g .* f)`. To get a
*density* for a diffusion, divide by the cell widths of the grid — essential when the
grid is non-uniform.

The keywords `δ` and `rebirth` have the same meaning as in the matrix form; `rebirth`
is an array shaped like `size(X)`. Remaining keyword arguments are forwarded to
`generator(X)` — e.g. `check = :warn` for a [`MultivariateDiffusionProcess`](@ref).
"""
function stationary_distribution(X::ContinuousTimeMarkovProcess; δ = 0.0, rebirth = nothing, ψ = nothing, kwargs...)
    rebirth = _resolve_rebirth_argument(rebirth, ψ, Ones(length(X)))
    size(rebirth) == size(X) || size(rebirth) == (length(X),) ||
        throw(DimensionMismatch("`rebirth` has size $(size(rebirth)) but the state space has size $(size(X))"))
    _reshape_state_output(X, stationary_distribution(generator(X; kwargs...); δ = δ, rebirth = vec(rebirth)))
end

"""
    feynman_kac(X::ContinuousTimeMarkovProcess, ts; f = nothing, ψ = nothing, v = nothing, direction = :backward, kwargs...)

Solve the Feynman–Kac PDE associated with a Markov process `X` on the time grid
`ts`, using implicit Euler time steps.

For the process form, `f`, `ψ`, and `v` are arrays shaped like the state space (`size(X)`);
`f` and `v` may also carry a trailing time dimension of length `length(ts)`. The result has
one slice per date in `ts`. Remaining keyword arguments are forwarded to `generator(X)` —
e.g. `check = :warn` for a [`MultivariateDiffusionProcess`](@ref).
"""
function feynman_kac(X::ContinuousTimeMarkovProcess, ts; f = nothing, ψ = nothing, v = nothing, direction = :backward, kwargs...)
    𝔸 = generator(X; kwargs...)
    f_flat = f === nothing ? zeros(eltype(𝔸), length(X)) :
        _flatten_state_time_argument(X, f, ts, :f)
    ψ_flat = ψ === nothing ? zeros(eltype(𝔸), length(X)) :
        _flatten_state_vector_argument(X, ψ, :ψ)
    v_flat = v === nothing ? zeros(eltype(𝔸), length(X)) :
        _flatten_state_time_argument(X, v, ts, :v)
    u = feynman_kac(𝔸, ts; f = f_flat, ψ = ψ_flat, v = v_flat, direction = direction)
    return _reshape_state_time_output(X, u)
end
