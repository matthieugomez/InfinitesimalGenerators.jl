#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + v2 * ∂∂(f)
Note that
𝔸'g = v_0 * g - ∂(v1 * g) + ∂∂(v2 * g)

========================================================================================#
function operator!(𝔸, Δ, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    # The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(𝔸, 0)
    @inbounds for i in 1:n
        if v1[i] >= 0
            𝔸[min(i + 1, n), i] += v1[i] * invΔxp[i]
            𝔸[i, i] -= v1[i] * invΔxp[i]
        else
            𝔸[i, i] += v1[i] * invΔxm[i]
            𝔸[max(i - 1, 1), i] -= v1[i] * invΔxm[i]
        end
        𝔸[max(i - 1, 1), i] += v2[i] * invΔxm[i] * invΔx[i]
        𝔸[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        𝔸[min(i + 1, n), i] += v2[i] * invΔxp[i] * invΔx[i]
    end
    c = sum(𝔸, dims = 1)
    for i in 1:n
        𝔸[i, i] += v0[i] - c[i]
    end
    adjoint(𝔸)
end


#========================================================================================

Markov Diffusion 
dx = μ(x) dt + σ(x) dZ_t

========================================================================================#

mutable struct MarkovDiffusion <: MarkovProcess
    x::AbstractVector{<:Real}
    μx::AbstractVector{<:Real}
    σx::AbstractVector{<:Real}
    𝔸::Tridiagonal
    Δ::Tuple{<:AbstractVector, <:AbstractVector, <:AbstractVector, <:AbstractVector}
end

function MarkovDiffusion(x::AbstractVector{<:Real}, μx::AbstractVector{<:Real}, σx::AbstractVector{<:Real})
    length(x) == length(μx) || error("Vector for grid, drift, and volatility should have the same size")
    length(μx) == length(σx) || error("Vector for grid, drift, and volatility should have the same size")
    n = length(x)
    𝔸 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    Δ = make_Δ(x)
    MarkovDiffusion(x, μx, σx, 𝔸, Δ)
end

function make_Δ(x)
    n = length(x)
    Δxm = zero(x)
    Δxm[1] = x[2] - x[1]
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zero(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δxp[end] = x[n] - x[n-1]
    Δx = (Δxm .+ Δxp) / 2
    return x, 1 ./ Δx, 1 ./ Δxm, 1 ./ Δxp
end

generator!(X::MarkovDiffusion) = operator!(X.𝔸, X.Δ, Zeros(length(X.x)), X.μx, 0.5 * X.σx.^2)

Base.length(X::MarkovDiffusion) = length(X.x)


function ∂(X::MarkovDiffusion, f::AbstractVector)
    operator!(X.𝔸, X.Δ, Zeros(length(X.x)), Ones(length(X.x)), Zeros(length(X.x))) * f
end

# it's important to take 1e-6 to have the right tail index of multiplicative functional (see tests)
function OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, 
    xmin = quantile(Normal(xbar, σ / sqrt(2 * κ)), p), xmax = quantile(Normal(xbar, σ / sqrt(2 * κ)), 1 - p))
    if xmin > 0
        x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    else
        x = range(xmin, stop = xmax, length = length)
    end
    μx = κ .* (xbar .- x)
    σx = σ .* Ones(Base.length(x))
    MarkovDiffusion(x, μx, σx)
end

function CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, α = 2 * κ * xbar / σ^2, β = σ^2 / (2 * κ), xmin = quantile(Gamma(α, β), p), xmax = quantile(Gamma(α, β), 1 - p), pow = 2)
    x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    μx = κ .* (xbar .- x)
    σx = σ .* sqrt.(x)
    MarkovDiffusion(x, μx, σx)
end
#========================================================================================

Multiplicative functional M for diffusion
dM/M = μM(x) dt + σM(x) dZt
dx = μ(x) dt + σ(x) dZt

========================================================================================#

mutable struct MultiplicativeFunctionalDiffusion <: MultiplicativeFunctional
    X::MarkovDiffusion
    μM::AbstractVector{<:Number}
    σM::AbstractVector{<:Number}
    ρ::Number
    δ::Number
end

function MultiplicativeFunctionalDiffusion(X::MarkovDiffusion, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; ρ::Number = 0.0, δ::Number = 0.0)
    length(X.x) == length(μM) || error("Vector for grid and μM should have the same size")
    length(X.x) == length(σM) || error("Vector for grid and σM should have the same size")
    MultiplicativeFunctionalDiffusion(X, μM, σM, ρ, δ)
end

function generator!(M::MultiplicativeFunctionalDiffusion)
    ξ -> operator!(M.X.𝔸, M.X.Δ, ξ .* M.μM .+ 0.5 * ξ * (ξ - 1) .* M.σM.^2 .- M.δ,  M.X.μx .+ ξ .* M.σM .* M.ρ .* M.X.σx, 0.5 * M.X.σx.^2)
end
Base.length(M::MultiplicativeFunctionalDiffusion) = length(M.X)