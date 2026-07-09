# Processes and generators

InfinitesimalGenerators.jl represents a continuous-time Markov process by a process object. Each process object knows how to construct its finite-grid generator matrix, and the rest of the package works with that generator.

## Common interface

All process types implement the same basic interface:

- [`generator(X)`](@ref generator): the generator (transition-rate) matrix of the discretized process.
- [`state_space(X)`](@ref state_space): a tuple of state-space axes.
- `size(X)`: the tensor shape of the state space.

All built-in process types subtype `ContinuousTimeMarkovProcess{N}`, where `N` is the number of tensor-product state-space axes. A scalar diffusion or finite-state chain has `N == 1`; a product of two scalar processes has `N == 2`. Arrays are flattened in Julia's column-major order when applying `generator(X)`.

`generator(X)` checks that the discretized operator is a valid Markov generator: rows must sum to zero and off-diagonal entries must be non-negative transition rates. Use `generator(X; check = :warn)` or `check = false` only when you intentionally want the raw finite-difference operator without a Markov-process interpretation.

```@example processes
using InfinitesimalGenerators
nothing # hide
```

## Finite-state chains

The simplest process is one whose generator is given directly. [`ContinuousTimeMarkovChain`](@ref)`(states, Q)` represents a continuous-time chain on `states` with generator matrix `Q` (rows sum to zero, non-negative off-diagonals), and `generator(Z)` is just `Q`. `ContinuousTimeMarkovChain(Q)`, without states, uses the indices `1:size(Q, 1)`.

```@example processes
states = [0.5, 1.5]
Q = [-0.1 0.1; 0.2 -0.2]
Z = ContinuousTimeMarkovChain(states, Q)
generator(Z)
```

## Scalar diffusions

[`DiffusionProcess`](@ref)`(x, μx, σx)` represents ``dx_t = \mu(x_t)dt + \sigma(x_t)dW_t``, given a strictly increasing grid `x` (possibly non-uniform) and the drift and volatility evaluated on it. Discretizing turns the diffusion into a finite-state chain jumping between neighboring grid points.

Drift is discretized by upwinding: forward differences where the drift is positive, backward differences where it is negative. Boundaries are reflecting, so the discretized operator is always a valid generator matrix (rows sum to zero, non-negative off-diagonals).

```@example processes
x = range(-1, 1, length = 5)
X = DiffusionProcess(x, -0.03 .* x, 0.01 .* ones(length(x)))
@assert ndims(X) == 1 # hide
generator(X)
```

In the finite-difference literature, a discretization with non-negative off-diagonal weights is called a *monotone* scheme — one of the conditions (together with consistency and stability) under which finite-difference solutions of HJB equations converge to the viscosity solution, and the reason upwinding is the standard discretization there (see [EconPDEs' discussion of upwinding](https://matthieugomez.github.io/EconPDEs.jl/dev/getting_started/#Upwinding)). "Monotone scheme" and "valid Markov generator" are the same sign condition seen from two sides, which is why the two packages produce consistent results on the same grid.

The two workhorse scalar diffusions come with constructors — [`OrnsteinUhlenbeck`](@ref) and [`CoxIngersollRoss`](@ref) — that choose the grid automatically:

```@example processes
Xou = OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.03, σ = 0.01)   # dx = -κ (x - xbar) dt + σ dW
Xcir = CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 0.07)    # dx = -κ (x - xbar) dt + σ √x dW
nothing # hide
```

By default the grid spans the `p` and `1 - p` quantiles of the stationary distribution with `length` points; pass `length`, `xmin`, `xmax`, or `pow` (grid-spacing power) to override. The default `p = 1e-10` errs on the side of a wide grid: reflecting boundaries distort solutions near the edges of the grid (see [Computing expectations (Kolmogorov backward)](expectations.md)), and a wide grid pushes that distortion where the process almost never goes — which matters especially when computing [tail indices](tail_index.md).

!!! note "Upwinding keeps the probabilistic interpretation when the state space is discretized"
    A `DiffusionProcess` is discretized so that its generator matrix has non-negative off-diagonal elements and rows summing to zero. This allows a probabilistic interpretation of the discretized operator — it can be seen as a continuous-time Markov chain, with non-negative rates of jumping between neighboring grid points. This is what upwinding is for: differencing the drift in the direction it points — forward where ``\mu > 0``, backward where ``\mu < 0`` — keeps those rates non-negative on any grid, where a centered difference would not.

## Combining processes

Scalar processes can be combined into higher-dimensional processes. These processes live on tensor-product state spaces: `size(X)` is the shape of the grid, and `state_space(X)` returns one axis per tensor dimension.

#### Products of independent processes

[`ProductProcess`](@ref) combines independent processes into their joint process. The generator is assembled from Kronecker sums, and the state space is the tensor product in argument order:

```@example processes
Y = ProductProcess(X, Z)
generator(Y)
```

#### Regime switching

[`SwitchingProcess`](@ref)`(Z, Xs)` drops the independence assumption: the dynamics depend on the autonomous modulating process `Z`. In the `i`-th state of `Z`, the process follows `Xs[i]`; all modulated processes must share the same grid.

`Z` can be any process in the package. A chain gives regime switching, while a diffusion gives continuous modulation. This is the natural representation of a solved model with a discrete state — see [Application to a consumption-saving problem](hjb.md), where wealth drifts at a policy-implied rate that switches with income:

```@example processes
xs = only(state_space(X))
Xlow  = DiffusionProcess(xs, -0.03 .* xs .- 0.01, 0.01 .* ones(length(xs)))
Xhigh = DiffusionProcess(xs, -0.03 .* xs .+ 0.01, 0.01 .* ones(length(xs)))
S = SwitchingProcess(Z, [Xlow, Xhigh])
generator(S)
```

## Correlated multivariate diffusions

[`MultivariateDiffusionProcess`](@ref) handles correlated diffusions on a tensor-product grid, with drift, variance, and covariance terms already evaluated on it (scalars or arrays), for example after solving an HJB:

Own-drift and own-variance terms always produce a valid generator under the built-in upwind scheme. The directional cross-derivative stencil used for correlated states is more restrictive: it is monotone only under a grid-scaled diagonal-dominance condition. For two states `x` and `y`, a useful local rule of thumb is

```julia
abs(covxy) <= variance.x * Δy / Δx
abs(covxy) <= variance.y * Δx / Δy
```

On equally spaced grids this becomes `abs(covxy) <= min(variance.x, variance.y)`, which is stronger than positive semidefiniteness of the covariance matrix. If `generator(M)` throws a negative-off-diagonal error, the most common fix is to rescale or refine the grids so `Δx / Δy` is closer to `sqrt(variance.x / variance.y)` in the region where the covariance is large. Equivalently, transform states so their local volatilities are more balanced, or reduce the covariance.

```@example processes
xs = range(-1, 1, length = 3)
ys = range(0, 2, length = 3)
grid = (; x = xs, y = ys)

μx = repeat(-0.1 .* xs, 1, length(ys))
μy = repeat(reshape(1 .- ys, 1, :), length(xs), 1)
σx = 0.2 .* ones(length(xs), length(ys))
σy = 0.3 .* ones(length(xs), length(ys))
covxy = 0.01 .* ones(length(xs), length(ys))

M = MultivariateDiffusionProcess(grid;
    drift = (; x = μx, y = μy),
    variance = (; x = σx .^ 2, y = σy .^ 2),
    covariance = (; xy = covxy),
)

𝔸 = generator(M)
𝔸
```
