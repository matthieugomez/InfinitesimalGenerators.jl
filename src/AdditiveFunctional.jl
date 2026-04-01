abstract type AdditiveFunctional end

# Should define generator which is a transition matrix T such that
# Tf = lim_{t→0} E[e^{ξ * m_t} f(x_t)|x_0=x]/t


"""
Compute the long run cgf(m), i.e. the function
ξ ⭌ lim_{t→∞} log(E[e^{ξ * m_t}])/t
"""
function cgf(m::AdditiveFunctional; eigenvector = :right, r0 = Ones(length(m.X.x)), η0 = nothing)
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
compute the tail index of the stationary distribution of e^{m}, i.e.
ζ such that cgf(m)(ζ) = δ
"""
function tail_index(m::AdditiveFunctional; δ = 0, verbose = false, r0 = nothing, xatol = 1e-4, kwargs...)
    r0 !== nothing && Base.depwarn("the `r0` keyword argument is deprecated and has no effect", :tail_index)
    fzero((1e-5, 1e3); xatol = xatol, kwargs...) do ξ
        η, _ = cgf(m)(ξ)
        verbose && @show (:LR, ξ, η)
        return η - δ
    end
end

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

mutable struct AdditiveFunctionalDiffusion{TX <: DiffusionProcess, Tμ <: AbstractVector{<:Number}, Tσ <: AbstractVector{<:Number}, TR <: Number} <: AdditiveFunctional
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
