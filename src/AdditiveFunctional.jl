abstract type AdditiveFunctional end
# Should define generator which is a transition matrix T such that
# Tf = lim_{t→0} E[e^{ξ * m_t} f(x_t)|x_0=x]/t


"""
Compute the long run cgf(m), i.e. the function
ξ ⭌ lim_{t→∞} log(E[e^{ξ * m_t}])/t
"""
function cgf(m::AdditiveFunctional; eigenvector = :right, r0 = Ones(length(m.X.x)))
    ξ -> principal_eigenvalue(generator(m)(ξ); eigenvector = eigenvector, r0 = r0)
end

"""
compute the tail index of the stationary distribution of e^{m}, i.e.
ζ such that cgf(m)(ζ) = δ
"""
function tail_index(m::AdditiveFunctional; δ = 0.0, verbose = false, r0 = ones(length(m.X.x)), xatol = 1e-4, kwargs...)
    r0 = deepcopy(r0)
    fzero((1e-5, 1e3); xatol = xatol, kwargs...) do ξ
        _, η, f = cgf(m; r0 = r0)(ξ)
        copyto!(r0, f)
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
    𝕋::Tridiagonal
end

function AdditiveFunctionalDiffusion(X::DiffusionProcess, μm::AbstractVector{<:Number}, σm::AbstractVector{<:Number}; ρ::Number = 0.0)
    AdditiveFunctionalDiffusion(X, μm, σm, ρ, deepcopy(X.𝕋))
end

function generator(M::AdditiveFunctionalDiffusion)
    ξ -> Diagonal(ξ .* M.μm .+ 0.5 * ξ^2 .* M.σm.^2) + generator!(M.𝕋, M.X.x, M.X.μx .+ ξ .* M.ρ .* M.σm .* M.X.σx, M.X.σx)
end
