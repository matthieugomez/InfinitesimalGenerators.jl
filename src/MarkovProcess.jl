
abstract type MarkovProcess end
# This type sould define generator(), which is a transition matrix 𝕋 such that
# 𝕋f = lim_{t→0} E[f(x_t)|x_0=x]/t

"""
 computes the stationary distribution corresponding to the MarkovProcess X
"""
function stationary_distribution(X::MarkovProcess; δ = 0.0, ψ = Zeros(length(X.x)))
    δ >= 0 ||      throw(ArgumentError("δ needs to be positive"))
    if δ > 0
        return clean_eigenvector_left((δ * I - generator(X)') \ (δ * ψ))
    else
        g, η, _ = principal_eigenvalue(generator(X); eigenvector = :left)
        abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
        return g
    end
end


#========================================================================================

Application for Diffusion Process x_t defined by:
dx = μ(x) dt + σ(x) dZ_t

========================================================================================#

mutable struct DiffusionProcess <: MarkovProcess
    x::AbstractVector{<:Real}
    μx::AbstractVector{<:Real}
    σx::AbstractVector{<:Real}
    𝕋::Tridiagonal
end

function DiffusionProcess(x::AbstractVector{<:Real}, μx::AbstractVector{<:Real}, σx::AbstractVector{<:Real})
    length(x) == length(μx) || error("Vector for grid, drift, and volatility should have the same size")
    length(μx) == length(σx) || error("Vector for grid, drift, and volatility should have the same size")
    n = length(x)
    𝕋 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    generator!(𝕋, x, μx, σx)
    DiffusionProcess(x, μx, σx, 𝕋)
end

generator(X::DiffusionProcess) = X.𝕋
state_space(X::DiffusionProcess) = X.x

# Compute the generator 
function generator!(𝕋, x, μx::AbstractVector, σx::AbstractVector)
    # The key is that sum of each row = 0.0 and off diagonals are positive
    n = length(x)
    fill!(𝕋, 0)
    @inbounds for i in 1:n
        Δxp =x[min(i, n-1)+1] - x[min(i, n-1)]
        Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
        Δx = (Δxm + Δxp) / 2
        # upwinding to ensure off diagonals are posititive
        if μx[i] >= 0
            𝕋[i, min(i + 1, n)] += μx[i] / Δxp
            𝕋[i, i] -= μx[i] / Δxp
        else
            𝕋[i, i] += μx[i] / Δxm
            𝕋[i, max(i - 1, 1)] -= μx[i] / Δxm
        end
        𝕋[i, max(i - 1, 1)] += 0.5 * σx[i]^2 / (Δxm * Δx)
        𝕋[i, i] -= 0.5 * σx[i]^2 * 2 / (Δxm * Δxp)
        𝕋[i, min(i + 1, n)] += 0.5 * σx[i]^2 / (Δxp * Δx)
    end
    # ensure machine precision
    c = sum(adjoint(𝕋), dims = 1)
    for i in 1:n
        𝕋[i, i] -= c[i]
    end
    return 𝕋
end




# Special cases.
# it's important to take low p to have the right tail index of Additive functional (see tests)
function OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, 
    xmin = quantile(Normal(xbar, σ / sqrt(2 * κ)), p), xmax = quantile(Normal(xbar, σ / sqrt(2 * κ)), 1 - p))
    if xmin > 0
        x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    else
        x = range(xmin, stop = xmax, length = length)
    end
    μx = κ .* (xbar .- x)
    σx = σ .* Ones(Base.length(x))
    DiffusionProcess(x, μx, σx)
end

function CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, α = 2 * κ * xbar / σ^2, β = σ^2 / (2 * κ), xmin = quantile(Gamma(α, β), p), xmax = quantile(Gamma(α, β), 1 - p), pow = 2)
    # check 0 is not attainable
    @assert (2 * κ * xbar) / σ^2 > 1
    x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    μx = κ .* (xbar .- x)
    σx = σ .* sqrt.(x)
    DiffusionProcess(x, μx, σx)
end


function ∂(X::DiffusionProcess)
    generator!(deepcopy(X.𝕋), X.x, X.μx, Zeros(length(X.x))) ./ X.μx
end
