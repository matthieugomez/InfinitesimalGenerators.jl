"""
    DiffusionProcess(x, μx, σx)

Returns the diffusion process

    dx_t = μ(x_t) dt + σ(x_t) dZ_t

discretized on the strictly increasing grid `x` (possibly non-uniform), where `μx` and `σx`
are the drift and volatility evaluated on the grid. The drift is discretized by upwinding
and the boundaries are reflecting, so [`generator`](@ref) always returns a valid generator
(transition-rate) matrix.
"""
struct DiffusionProcess{TX <: AbstractVector{<:Real}, Tμ <: AbstractVector{<:Real}, Tσ <: AbstractVector{<:Real}} <: ContinuousTimeMarkovProcess{1}
    x::TX
    μx::Tμ
    σx::Tσ
    function DiffusionProcess(x::TX, μx::Tμ, σx::Tσ) where {TX <: AbstractVector{<:Real}, Tμ <: AbstractVector{<:Real}, Tσ <: AbstractVector{<:Real}}
        length(x) == length(μx) == length(σx) || throw(ArgumentError("Vector for grid, drift, and volatility should have the same size"))
        length(x) >= 2 || throw(ArgumentError("State grid must contain at least two points"))
        all(x[i] < x[i + 1] for i in 1:(length(x) - 1)) || throw(ArgumentError("State grid must be strictly increasing"))
        new{TX, Tμ, Tσ}(x, μx, σx)
    end
end

state_space(X::DiffusionProcess) = (X.x,)

"""
    generator(X::DiffusionProcess)

Returns the discretized infinitesimal generator of the diffusion,

    𝔸: f ↦ μx * ∂f + 0.5 * σx^2 * ∂²f,

acting on functions with reflecting boundary conditions (∂f = 0 at the edges of the grid).
Its transpose is the discretized forward (Fokker–Planck) operator,

    𝔸': g ↦ -∂(μx * g) + 0.5 * ∂²(σx^2 * g),

acting on densities with zero-flux boundary conditions (-μx * g + 0.5 * ∂(σx^2 * g) = 0 at
the edges of the grid).
"""
function generator(X::DiffusionProcess)
    generator(X.x, X.μx, X.σx)
end

function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    n = length(x)
    𝔸 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    @inbounds for i in 1:n
        Δxp = x[min(i, n-1)+1] - x[min(i, n-1)]
        Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
        Δx = (Δxm + Δxp) / 2
        # upwinding with reflecting boundaries: outward boundary drift is dropped
        if μx[i] >= 0
            if i < n
                𝔸[i, i + 1] += μx[i] / Δxp
                𝔸[i, i] -= μx[i] / Δxp
            end
        elseif i > 1
            𝔸[i, i] += μx[i] / Δxm
            𝔸[i, i - 1] -= μx[i] / Δxm
        end
        𝔸[i, max(i - 1, 1)] += 0.5 * σx[i]^2 / (Δxm * Δx)
        𝔸[i, i] -= 0.5 * σx[i]^2 * 2 / (Δxm * Δxp)
        𝔸[i, min(i + 1, n)] += 0.5 * σx[i]^2 / (Δxp * Δx)
    end
    # ensure rows sum to zero with machine precision
    c = sum(𝔸, dims = 2)
    for i in 1:n
        𝔸[i, i] -= c[i]
    end
    return 𝔸
end

"""
    Returns the discretized version of the operator ∂

        ∂: f ↦ ∂f

    The scheme is upwind with respect to the drift (forward where μx ≥ 0,
    backward where μx < 0), matching the discretization used in `generator`.
    At interior nodes where the drift is exactly zero — where upwinding has no
    preferred direction — a central difference is used instead. (Building the
    matrix directly avoids the `Diagonal(μx) \\ generator(…)` division, which
    would produce `NaN` rows wherever μx = 0.)
"""
function ∂(X::DiffusionProcess)
    x, μx = X.x, X.μx
    n = length(x)
    D = Tridiagonal(zeros(n - 1), zeros(n), zeros(n - 1))
    @inbounds for i in 1:n
        Δxp = x[min(i, n - 1) + 1] - x[min(i, n - 1)]
        Δxm = x[max(i - 1, 1) + 1] - x[max(i - 1, 1)]
        if (μx[i] == 0) && (1 < i < n)
            # zero drift: central difference rather than 0/0
            D[i, i - 1] -= 1 / (Δxm + Δxp)
            D[i, i + 1] += 1 / (Δxm + Δxp)
        elseif μx[i] >= 0
            # forward (upwind)
            if i < n
                D[i, i + 1] += 1 / Δxp
                D[i, i]     -= 1 / Δxp
            end
        elseif i > 1
            # backward (upwind)
            D[i, i]     += 1 / Δxm
            D[i, i - 1] -= 1 / Δxm
        end
    end
    return D
end

"""
    OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, xmin, xmax, pow = 1)

Returns the Ornstein–Uhlenbeck process

    dx_t = -κ * (x_t - xbar) * dt + σ * dZ_t

discretized as a [`DiffusionProcess`](@ref) on an automatically chosen grid.

By default the grid has `length` points spanning the `p` and `1 - p` quantiles of the
stationary distribution `N(xbar, σ^2 / 2κ)`; pass `xmin` and `xmax` to override the limits.
`pow` controls the grid spacing (when `xmin > 0`, points are uniform in `x^(1/pow)`, so
`pow > 1` concentrates points near `xmin`). The default `p = 1e-10` is deliberately extreme:
reflecting boundaries distort solutions near the edges of the grid, and a wide grid pushes
that distortion where the process never goes — which matters especially for [`tail_index`](@ref).
"""
function OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, p = 1e-10, length = 100,
    xmin = quantile(Normal(xbar, σ / sqrt(2 * κ)), p), xmax = quantile(Normal(xbar, σ / sqrt(2 * κ)), 1 - p), pow = 1)
    # it's important to take low p to have the right tail index of Additive functional
    if xmin > 0
        x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    else
        x = range(xmin, stop = xmax, length = length)
    end
    DiffusionProcess(x, κ .* (xbar .- x), σ * Ones(Base.length(x)))
end

"""
    CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, xmin, xmax, pow = 2)

Returns the Cox–Ingersoll–Ross process

    dx_t = -κ * (x_t - xbar) * dt + σ * sqrt(x_t) * dZ_t

discretized as a [`DiffusionProcess`](@ref) on an automatically chosen grid. Requires the
Feller condition `2κ * xbar / σ^2 > 1`, so that 0 is not attainable.

By default the grid has `length` points spanning the `p` and `1 - p` quantiles of the
stationary Gamma distribution; pass `xmin` and `xmax` to override the limits. `pow` controls
the grid spacing (points are uniform in `x^(1/pow)`, so the default `pow = 2` concentrates
points near zero, where the volatility is smallest).
"""
function CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, α = 2 * κ * xbar / σ^2, β = σ^2 / (2 * κ), xmin = quantile(Gamma(α, β), p), xmax = quantile(Gamma(α, β), 1 - p), pow = 2)
    # check 0 is not attainable
    (2 * κ * xbar) / σ^2 > 1 || throw(ArgumentError("Feller condition not satisfied: 2κx̄/σ² must be > 1"))
    x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    DiffusionProcess(x, κ .* (xbar .- x), σ .* sqrt.(x))
end
