"""
    AdditiveFunctional

Abstract supertype for additive functionals `m` of a Markov process — cumulative quantities
whose increments depend on the state, such as cumulative log growth. The defining verb is
[`tilted_generator`](@ref)`(m, ξ)`; the operators [`cgf`](@ref), [`tail_index`](@ref), and
[`principal_eigenvalue`](@ref) are computed from it.
"""
abstract type AdditiveFunctional end

"""
    tilted_generator(m::AdditiveFunctional, ξ)

Return the tilted generator matrix `𝔸_ξ` of the additive functional `m`, i.e. the matrix
representation of the operator

    𝔸_ξ: f ↦ lim_{t→0} (E[e^{ξ mₜ} f(xₜ) | x₀ = x] - f(x)) / t.

`tilted_generator(m, 0)` is the generator of the underlying process. Everything the package
computes about additive functionals derives from this matrix: [`cgf`](@ref) and
[`tail_index`](@ref) use its principal eigenvalue, and it can be passed directly to
[`feynman_kac`](@ref) to compute finite-horizon moments `E[e^{ξ mₜ} ψ(xₜ)]`. To define a
custom additive functional, subtype `AdditiveFunctional` and define this function.

A tilted generator is not a generator in the sense of [`check_generator`](@ref): its
off-diagonal entries must be nonnegative (a *Metzler* matrix — equivalently, a generator
plus a diagonal), but its rows need not sum to zero. Metzlerity is what makes the principal
eigenvalue real and its eigenvector positive; a custom implementation must preserve it
([`principal_eigenvalue`](@ref) checks it on entry and warns if it fails).
"""
function tilted_generator end

# fallback for subtypes that define the legacy one-argument closure form
tilted_generator(m::AdditiveFunctional, ξ::Number) = tilted_generator(m)(ξ)

"""
    cgf(m::AdditiveFunctional, ξ; r0 = nothing, η0 = nothing)

Return the long-run scaled cumulant generating function of `m` evaluated at `ξ`,

    Λ(ξ) = lim_{t→∞} log(E[e^{ξ mₜ}]) / t,

computed as the principal eigenvalue of the tilted generator [`tilted_generator`](@ref)`(m, ξ)`
by inverse iteration (Hansen and Scheinkman 2009). `r0` is an initial guess for the
eigenvector (defaults to a vector of ones) and `η0` an initial guess for the eigenvalue.
Use [`cgf_eigenvector`](@ref) to also obtain the associated eigenvector.
"""
function cgf(m::AdditiveFunctional, ξ::Number; r0 = nothing, η0 = nothing)
    𝔸 = tilted_generator(m, ξ)
    principal_eigenvalue(𝔸; r0 = r0 === nothing ? Ones(size(𝔸, 1)) : r0, η0 = η0)[1]
end

"""
    cgf_eigenvector(m::AdditiveFunctional, ξ, side = :right; r0 = nothing, η0 = nothing)

Return the pair `(Λ(ξ), vector)`: the long-run scaled cumulant generating function of `m`
at `ξ` (as in [`cgf`](@ref)) together with the associated principal eigenvector of the
tilted generator.

With `side = :right`, the right eigenvector — the Hansen–Scheinkman eigenfunction. With
`side = :left`, the left eigenvector, normalized to sum to one — the stationary
distribution of the twisted process, i.e. the distribution over states from which the
long-run growth of `E[e^{ξ mₜ}]` is achieved.
"""
function cgf_eigenvector(m::AdditiveFunctional, ξ::Number, side::Symbol = :right; r0 = nothing, η0 = nothing)
    𝔸 = tilted_generator(m, ξ)
    r0 = r0 === nothing ? Ones(size(𝔸, 1)) : r0
    if side == :right
        return principal_eigenvalue(𝔸; r0 = r0, η0 = η0)
    elseif side == :left
        η, l = principal_eigenvalue(𝔸'; r0 = r0, η0 = η0)
        return η, l ./ sum(l)
    else
        throw(ArgumentError("side must be :right or :left"))
    end
end

# deprecated closure form: cgf(m)(ξ) returning (η, eigenvector)
function cgf(m::AdditiveFunctional; eigenvector = :right, r0 = nothing, η0 = nothing)
    Base.depwarn("`cgf(m)(ξ)` is deprecated; use `cgf(m, ξ)` for the value and " *
        "`cgf_eigenvector(m, ξ, :right)` or `cgf_eigenvector(m, ξ, :left)` for eigenvectors", :cgf)
    ξ -> begin
        eigenvector in (:right, :left) ||
            throw(ArgumentError("the keyword argument eigenvector can only take the value :right or :left"))
        cgf_eigenvector(m, ξ, eigenvector; r0 = r0, η0 = η0)
    end
end

"""
    tail_index(m::AdditiveFunctional; δ = 0, bracket = (1e-5, 1e3), xatol = 1e-4)

Compute the tail index of the stationary distribution of `e^m` when units die (are reset) at
rate `δ`, i.e. the ζ such that `cgf(m, ζ) = δ`.

The root is searched for in `bracket`; if `cgf(m, ξ) - δ` has the same sign at both ends,
an `ArgumentError` reports the two values so the bracket can be moved. Remaining keyword
arguments are passed to `Roots.fzero`.
"""
function tail_index(m::AdditiveFunctional; δ = 0, verbose = false, r0 = nothing, xatol = 1e-4, bracket = (1e-5, 1e3), kwargs...)
    r0 !== nothing && Base.depwarn("the `r0` keyword argument is deprecated and has no effect", :tail_index)
    Λ = ξ -> begin
        η = cgf(m, ξ)
        verbose && @show (:LR, ξ, η)
        return η - δ
    end
    ξlo, ξhi = bracket
    flo, fhi = Λ(ξlo), Λ(ξhi)
    flo * fhi <= 0 ||
        throw(ArgumentError("`cgf(m, ξ) - δ` has the same sign at both ends of `bracket = $bracket` " *
            "($flo at ξ = $ξlo and $fhi at ξ = $ξhi), so no tail index was bracketed. " *
            "Pass a `bracket` whose endpoints straddle the root of `cgf(m, ξ) = δ`."))
    fzero(Λ, ξlo, ξhi; xatol = xatol, kwargs...)
end

"""
    tail_index(μ::Number, σ::Number; δ = 0)

Closed form for constant coefficients: the tail index of the stationary distribution of a
size `w` growing as `dw/w = μ dt + σ dZ` (note that `μ` is the arithmetic growth rate of `w`
itself, equal to the drift of `log w` plus `σ^2/2`) with death rate `δ`.
"""
function tail_index(μ::Number, σ::Number; δ::Number = 0)
    if σ > 0
        (1 - 2 * μ / σ^2 + sqrt((1- 2 * μ / σ^2)^2 + 8 * δ / σ^2)) / 2
    else
        δ / μ
    end
end

#=======================================================================================

The additive functional
    dm_t = μm(x_t)dt + σm(x_t)dZ^m_t
of a Markov process x, parameterized by canonical coefficients:
    drift      = μm
    variance   = σm²
    covariance = cov(dm, dx)   (requires a diffusion state; nothing = idiosyncratic noise)

========================================================================================#

struct AdditiveFunctionalProcess{TX <: ContinuousTimeMarkovProcess, Tμ <: AbstractVector{<:Number}, Tσ2 <: AbstractVector{<:Number}, TC} <: AdditiveFunctional
    X::TX
    drift::Tμ
    variance::Tσ2
    covariance::TC
end

_af_coefficient(x::Number, X, label::Symbol) = fill(float(x), length(X))

function _af_coefficient(x::AbstractArray, X, label::Symbol)
    (size(x) == size(X) || size(x) == (length(X),)) ||
        throw(DimensionMismatch("`$label` has size $(size(x)) but the state space has size $(size(X))"))
    collect(float.(vec(x)))
end

_af_covariance(::Nothing, X::ContinuousTimeMarkovProcess, variance) = nothing

function _af_covariance(c::Union{Number, AbstractArray}, X::DiffusionProcess, variance)
    cv = _af_coefficient(c, X, :covariance)
    for i in eachindex(cv)
        bound = variance[i] * X.σx[i]^2
        cv[i]^2 <= bound + sqrt(eps(Float64)) * max(1.0, bound) ||
            throw(ArgumentError("the joint covariance matrix of (dm, dx) must be positive " *
                "semidefinite: covariance(x)² ≤ variance(x) σ(x)² fails at grid index $i"))
    end
    cv
end

function _af_covariance(c::NamedTuple, X::MultivariateDiffusionProcess, variance)
    names = keys(X.grid)
    issubset(keys(c), names) ||
        throw(ArgumentError("`covariance` names $(collect(keys(c))) must be a subset of the grid names $(collect(names))"))
    shape = size(X)
    arrays = NamedTuple{names}(map(names) do name
        haskey(c, name) ? _mvd_coefficient_array(c[name], shape, Symbol(:covariance_, name)) : zeros(Float64, shape)
    end)
    # positive semidefiniteness of the joint covariance matrix of (dx, dm)
    N = length(names)
    covpairs = _mvd_covariance_pairs(names)
    varm = reshape(variance, shape)
    TΣ = promote_type(eltype(varm), map(eltype, Tuple(arrays))...,
        map(eltype, Tuple(X.variance))..., map(eltype, Tuple(X.covariance))...)
    Σ = zeros(TΣ, N + 1, N + 1)
    for I in CartesianIndices(shape)
        for d in 1:N
            Σ[d, d] = X.variance[names[d]][I]
        end
        for (d1, d2, key) in covpairs
            Σ[d1, d2] = Σ[d2, d1] = X.covariance[key][I]
        end
        for d in 1:N
            Σ[d, N + 1] = Σ[N + 1, d] = arrays[names[d]][I]
        end
        Σ[N + 1, N + 1] = varm[I]
        _check_psd(Σ) ||
            throw(ArgumentError("the joint covariance matrix of (dx, dm) must be positive semidefinite at grid index $I"))
    end
    arrays
end

_af_covariance(c, X::ContinuousTimeMarkovProcess, variance) =
    throw(ArgumentError("`covariance` requires a diffusion state: a vector or scalar for a " *
        "DiffusionProcess, or a NamedTuple such as `(; x = cx)` for a MultivariateDiffusionProcess"))

"""
    AdditiveFunctional(X::ContinuousTimeMarkovProcess; drift = 0, variance = 0, covariance = nothing)
    AdditiveFunctional(X, μm, σm; ρ = 0.0)

An additive functional `m` of the Markov process `X`, defined by

    dmₜ = μm(xₜ) dt + σm(xₜ) dZᵐₜ.

`X` can be any process in the package — diffusions, chains, products, switching processes.
`drift` (`= μm`) and `variance` (`= σm²`) are scalars or arrays shaped like `size(X)`.
`covariance` specifies `cov(dmₜ, dxₜ)` when the noise of `m` is correlated with the
innovations of the state; it requires a diffusion state — a scalar or vector for a
`DiffusionProcess`, a `NamedTuple` with grid names (e.g. `(; x = cx)`) for a
`MultivariateDiffusionProcess`. Omitted, the noise of `m` is idiosyncratic.

The positional form takes the volatility `σm` rather than the variance, and `ρ` — the
correlation between `dZᵐ` and the innovations of a univariate diffusion state — rather than
the covariance (`covariance = ρ σm σx`).

Typical use: `m = log w` for a size `w` (wealth, firm size) growing at a state-dependent
rate; then [`cgf`](@ref) gives its long-run CGF and [`tail_index`](@ref) the Pareto exponent
of its stationary distribution.
"""
function AdditiveFunctional(X::ContinuousTimeMarkovProcess; drift = 0, variance = 0, covariance = nothing)
    driftv = _af_coefficient(drift, X, :drift)
    variancev = _af_coefficient(variance, X, :variance)
    all(>=(0), variancev) || throw(ArgumentError("`variance` must be nonnegative"))
    AdditiveFunctionalProcess(X, driftv, variancev, _af_covariance(covariance, X, variancev))
end

_af_square(σ::Number) = float(σ)^2
_af_square(σ::AbstractArray) = float.(σ) .^ 2

function AdditiveFunctional(X::ContinuousTimeMarkovProcess, μm, σm)
    AdditiveFunctional(X; drift = μm, variance = _af_square(σm))
end

function AdditiveFunctional(X::DiffusionProcess, μm, σm; ρ = 0.0)
    if iszero(ρ)
        AdditiveFunctional(X; drift = μm, variance = _af_square(σm))
    else
        AdditiveFunctional(X; drift = μm, variance = _af_square(σm), covariance = ρ .* σm .* X.σx)
    end
end

function tilted_generator(m::AdditiveFunctionalProcess, ξ::Number)
    Diagonal(ξ .* m.drift .+ 0.5 .* ξ^2 .* m.variance) + _af_state_generator(m.X, m.covariance, ξ)
end

function _af_state_generator(X::ContinuousTimeMarkovProcess, ::Nothing, ξ::Number)
    𝔸 = generator(X)
    𝔸 isa Tridiagonal ? 𝔸 : sparse(𝔸)
end

function _af_state_generator(X::DiffusionProcess, c::AbstractVector, ξ::Number)
    generator(X.x, X.μx .+ ξ .* c, X.σx)
end

function _af_state_generator(X::MultivariateDiffusionProcess, c::NamedTuple, ξ::Number)
    names = keys(X.grid)
    drift = NamedTuple{names}(map(name -> X.drift[name] .+ ξ .* c[name], names))
    generator(MultivariateDiffusionProcess(X.grid; drift = drift, variance = X.variance, covariance = X.covariance))
end

generator(m::AdditiveFunctionalProcess) = tilted_generator(m, 1)

#=======================================================================================

Legacy type: AdditiveFunctional constructs an equivalent functional; kept for
backward compatibility.
dx_t = μ(x)dt + σ(x) dZ_t
dm_t = μm(x)dt + σm(x)dZ^m_t
with
corr(dZ^m_t, dZ_t) = ρ

========================================================================================#

"""
    AdditiveFunctionalDiffusion(X::DiffusionProcess, μm, σm; ρ = 0.0)

An additive functional `m` of the diffusion `X`, defined by

    dmₜ = μm(xₜ) dt + σm(xₜ) dZᵐₜ,        corr(dZᵐₜ, dZₜ) = ρ,

where `μm` and `σm` are vectors evaluated on the grid of `X`.

Legacy type: [`AdditiveFunctional`](@ref)`(X, μm, σm; ρ = ρ)` constructs an equivalent
functional, and is the preferred entry point; you rarely need to construct this type
directly.
"""
struct AdditiveFunctionalDiffusion{TX <: DiffusionProcess, Tμ <: AbstractVector{<:Number}, Tσ <: AbstractVector{<:Number}, TR <: Number} <: AdditiveFunctional
    X::TX
    μm::Tμ
    σm::Tσ
    ρ::TR
end

function AdditiveFunctionalDiffusion(X::TX, μm::Tμ, σm::Tσ; ρ::TR = 0.0) where {TX <: DiffusionProcess, Tμ <: AbstractVector{<:Number}, Tσ <: AbstractVector{<:Number}, TR <: Number}
    length(X.x) == length(μm) == length(σm) || throw(ArgumentError("Vector for grid, drift, and volatility should have the same size"))
    AdditiveFunctionalDiffusion{TX, Tμ, Tσ, TR}(X, μm, σm, ρ)
end

function generator(M::AdditiveFunctionalDiffusion)
    Diagonal(M.μm .+ 0.5 .* M.σm.^2) + generator(M.X.x, M.X.μx .+ M.ρ .* M.σm .* M.X.σx, M.X.σx)
end

function tilted_generator(M::AdditiveFunctionalDiffusion, ξ::Number)
    Diagonal(ξ .* M.μm .+ 0.5 .* ξ^2 .* M.σm.^2) + generator(M.X.x, M.X.μx .+ ξ .* M.ρ .* M.σm .* M.X.σx, M.X.σx)
end
