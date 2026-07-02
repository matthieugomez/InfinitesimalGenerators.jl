abstract type AdditiveFunctional end

# Should define generator which is a generator matrix T such that
# Tf = lim_{t→0} E[e^{ξ * m_t} f(x_t)|x_0=x]/t


"""
    cgf(m::AdditiveFunctional; eigenvector = :right)

Return the long-run scaled cumulant generating function of `m`, i.e. the function

    ξ ⭌ lim_{t→∞} log(E[e^{ξ * m_t}])/t

computed as the principal eigenvalue of the ξ-tilted generator (Hansen and Scheinkman 2009).
The returned function gives a tuple `(η, x)` with the eigenvalue `η` and the associated right
eigenvector (or, with `eigenvector = :left`, the left eigenvector normalized to sum to one).
"""
function cgf(m::AdditiveFunctional; eigenvector = :right, r0 = Ones(length(m.X)), η0 = nothing)
    ξ -> begin
        if eigenvector == :right
            principal_eigenvalue(tilted_generator(m)(ξ); r0 = r0, η0 = η0)
        elseif eigenvector == :left
            η, l = principal_eigenvalue(tilted_generator(m)(ξ)'; r0 = r0, η0 = η0)
            return η, l ./ sum(l)
        else
            throw(ArgumentError("the keyword argument eigenvector can only take the value :right or :left"))
        end
    end
end

"""
    tail_index(m::AdditiveFunctional; δ = 0)

Compute the tail index of the stationary distribution of `e^m` when units die (are reset) at
rate `δ`, i.e. the ζ such that `cgf(m)(ζ) = δ`.
"""
function tail_index(m::AdditiveFunctional; δ = 0, verbose = false, r0 = nothing, xatol = 1e-4, kwargs...)
    r0 !== nothing && Base.depwarn("the `r0` keyword argument is deprecated and has no effect", :tail_index)
    fzero((1e-5, 1e3); xatol = xatol, kwargs...) do ξ
        η, _ = cgf(m)(ξ)
        verbose && @show (:LR, ξ, η)
        return η - δ
    end
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

Diffusion Case
dx_t = μ(x)dt + σ(x) dZ_t
dm_t = μm(x)dt + σm(x)dZ^m_t
with
corr(dZ^m_t, dZ_t) = ρ

========================================================================================#

"""
    AdditiveFunctionalDiffusion(X::DiffusionProcess, μm, σm; ρ = 0.0)

An additive functional `m` of the diffusion `X`, defined by

    dmₜ = μm(xₜ) dt + σm(xₜ) dZᵐₜ,        corr(dZᵐₜ, dZₜ) = ρ,

where `μm` and `σm` are vectors evaluated on the grid of `X`. Typical use: `m = log w` for a
size `w` (wealth, firm size) growing at a state-dependent rate; then `cgf` gives its long-run
CGF and `tail_index` the Pareto exponent of its stationary distribution.
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

function tilted_generator(M::AdditiveFunctionalDiffusion)
    ξ -> generator(AdditiveFunctionalDiffusion(M.X, ξ .* M.μm, ξ .* M.σm, M.ρ))
end
