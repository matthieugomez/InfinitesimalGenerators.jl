abstract type AdditiveFunctional end
# Should define generator which is a transition matrix T such that
# Tf = lim_{t→0} E[e^{ξ * m_t} f(x_t)|x_0=x]/t


"""
Compute the long run cgf(m), i.e. the function
ξ ⭌ lim_{t→∞} log(E[e^{ξ * m_t}])/t
"""
function cgf(m::AdditiveFunctional; eigenvector = :right, r0 = Ones(length(m.X.x)))
    if eigenvector == :right
        ξ -> principal_eigenvalue(generator(m)(ξ); r0 = r0)
    elseif eigenvector == :left
        ξ -> begin
            η, l = principal_eigenvalue(generator(m)(ξ)'; r0 = r0)
            return η, l ./ sum(l)
        end
    else
        throw(ArgumentError("the keyword argument eigenvector can only take the value :right or :left"))
    end
end

"""
compute the tail index of the stationary distribution of e^{m}, i.e.
ζ such that cgf(m)(ζ) = δ
"""
function tail_index(m::AdditiveFunctional; δ = 0, verbose = false, r0 = Ones(length(m.X.x)), xatol = 1e-4, kwargs...)
    fzero((1e-5, 1e3); xatol = xatol, kwargs...) do ξ
        η, r0 = cgf(m; r0 = r0)(ξ)
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


#========================================================================================

Diffusion Case
dx_t = μ(x)dt + σ(x) dZ_t
dm_t = μm(x)dt + σm(x)dZ^m_t
with 
corr(dZ^m_t, dZ_t) = ρ

========================================================================================#

mutable struct AdditiveFunctionalDiffusion <: AdditiveFunctional
    X::DiffusionProcess
    μm::AbstractVector{<:Number}
    σm::AbstractVector{<:Number}
    ρ::Number
end

function AdditiveFunctionalDiffusion(X::DiffusionProcess, μm::AbstractVector{<:Number}, σm::AbstractVector{<:Number}; ρ::Number = 0.0)
    length(X.x) == length(μm) == length(σm) || throw(ArgumentError("Vector for grid, drift, and volatility should have the same size"))
    AdditiveFunctionalDiffusion(X, μm, σm, ρ)
end

function generator(M::AdditiveFunctionalDiffusion)
    ξ -> generator(M.X.x, ξ .* M.μm .+ 0.5 * ξ^2 .* M.σm.^2, M.X.μx .+ ξ .* M.ρ .* M.σm .* M.X.σx, M.X.σx)
end
