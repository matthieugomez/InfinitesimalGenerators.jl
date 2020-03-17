#========================================================================================

Compute the operator 
Tf = v_0 * f + v1 * ∂(f) + v2 * ∂∂(f)


Mote that
T'g = v_0 * g - ∂(v1 * g) + ∂∂(v2 * g)

========================================================================================#
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

function operator_second_order(x::AbstractVector, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    n = length(x)
    T = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    Δ = make_Δ(x)
    operator_second_order!(T, Δ, v0, v1, v2)
end

function operator_second_order!(T, Δ, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    # The key is that sum of each row = 0.0 and off diagonals are positive (singular M-matrix)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(T, 0)
    @inbounds for i in 1:n
        if v1[i] >= 0
            T[i, min(i + 1, n)] += v1[i] * invΔxp[i]
            T[i, i] -= v1[i] * invΔxp[i]
        else
            T[i, i] += v1[i] * invΔxm[i]
            T[i, max(i - 1, 1)] -= v1[i] * invΔxm[i]
        end
        T[i, max(i - 1, 1)] += v2[i] * invΔxm[i] * invΔx[i]
        T[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        T[i, min(i + 1, n)] += v2[i] * invΔxp[i] * invΔx[i]
    end
    c = sum(adjoint(T), dims = 1)
    for i in 1:n
        T[i, i] += v0[i] - c[i]
    end
    return T
end

#========================================================================================

Diffusion Process
dx = μ(x) dt + σ(x) dZ_t

========================================================================================#

mutable struct DiffusionProcess <: MarkovProcess
    x::AbstractVector{<:Real}
    μx::AbstractVector{<:Real}
    σx::AbstractVector{<:Real}
    T::Tridiagonal
    Δ::Tuple{<:AbstractVector, <:AbstractVector, <:AbstractVector, <:AbstractVector}
end

function DiffusionProcess(x::AbstractVector{<:Real}, μx::AbstractVector{<:Real}, σx::AbstractVector{<:Real})
    length(x) == length(μx) || error("Vector for grid, drift, and volatility should have the same size")
    length(μx) == length(σx) || error("Vector for grid, drift, and volatility should have the same size")
    n = length(x)
    T = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    Δ = make_Δ(x)
    operator_second_order!(T, Δ, Zeros(length(x)), μx, 0.5 * σx.^2)
    DiffusionProcess(x, μx, σx, T, Δ)
end

# it's important to take 1e-6 to have the right tail index of Additive functional (see tests)
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
    @assert (2 * κ * xbar)/σ^2 > 1
    x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    μx = κ .* (xbar .- x)
    σx = σ .* sqrt.(x)
    DiffusionProcess(x, μx, σx)
end

generator(X::DiffusionProcess) = X.T
state_space(X::DiffusionProcess) = X.x

function ∂(X::DiffusionProcess)
    T = deepcopy(X.T)
    operator_second_order!(T, X.Δ, Zeros(length(X.x)), Ones(length(X.x)), Zeros(length(X.x)))
end

#========================================================================================

Additive functional m for diffusion
dm = μm(x) dt + σm(x) dZt
dx = μ(x) dt + σ(x) dZt

========================================================================================#

mutable struct AdditiveFunctionalDiffusion <: AdditiveFunctional
    X::DiffusionProcess
    μm::AbstractVector{<:Number}
    σm::AbstractVector{<:Number}
    ρ::Number
    δ::Number
    T::Tridiagonal
end

function AdditiveFunctionalDiffusion(X::DiffusionProcess, μm::AbstractVector{<:Number}, σm::AbstractVector{<:Number}; ρ::Number = 0.0, δ::Number = 0.0)
    length(X.x) == length(μm) || error("Vector for grid and μm should have the same size")
    length(X.x) == length(σm) || error("Vector for grid and σm should have the same size")
    AdditiveFunctionalDiffusion(X, μm, σm, ρ, δ, deepcopy(X.T))
end

function generator(M::AdditiveFunctionalDiffusion)
    ξ -> operator_second_order!(M.T, M.X.Δ, ξ .* M.μm .+ 0.5 * ξ^2 .* M.σm.^2 .- M.δ,  M.X.μx .+ ξ .* M.ρ .* M.σm .* M.X.σx, 0.5 * M.X.σx.^2)
end