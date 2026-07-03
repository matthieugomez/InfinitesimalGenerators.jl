# Univariate processes

Every process in the package implements the same three verbs — [`generator(X)`](@ref generator) (the generator (transition-rate) matrix of the discretized process), [`state_space(X)`](@ref state_space) (a tuple of state-space axes), and `size(X)` (the tensor shape of the state space) — and everything else ([stationary distributions, expectations, tail indices](operators.md)) is built on them. Univariate processes subtype `ContinuousTimeMarkovProcess{1}`, so `ndims(X) == 1` and `state_space(X)` and `size(X)` are one-element tuples. This page tours the univariate process types.

## Finite-state continuous-time Markov chains

The simplest process is one whose generator is given directly: [`ContinuousTimeMarkovChain`](@ref)`(states, Q)` represents a continuous-time chain on `states` with generator matrix `Q` (rows sum to zero, non-negative off-diagonals), and `generator(Z)` is just `Q`. (`ContinuousTimeMarkovChain(Q)`, without states, uses the indices `1:size(Q, 1)`.)

```@example univariate
using InfinitesimalGenerators

states = [0.5, 1.5]
Q = [-0.1 0.1; 0.2 -0.2]
Z = ContinuousTimeMarkovChain(states, Q)
generator(Z)
```

## Diffusions

[`DiffusionProcess`](@ref)`(x, μx, σx)` represents ``dx_t = \mu(x_t)dt + \sigma(x_t)dZ_t``, given a strictly increasing grid `x` (possibly non-uniform) and the drift and volatility evaluated on it. Discretizing turns the diffusion into exactly the object above — a finite-state chain jumping between neighboring grid points. The drift is discretized by upwinding — forward differences where it is positive, backward where it is negative — and the boundaries are reflecting, so the discretized operator is always a valid generator matrix (rows sum to zero, non-negative off-diagonals):

```@example univariate
x = range(-1, 1, length = 100)
X = DiffusionProcess(x, -0.03 .* x, 0.01 .* ones(100))
@assert ndims(X) == 1 # hide
generator(X)
```

In the finite-difference literature, a discretization with non-negative off-diagonal weights is called a *monotone* scheme — the property that guarantees convergence to the right (viscosity) solution of HJB equations, and the reason upwinding is the standard discretization there (see [EconPDEs' discussion of upwinding](https://matthieugomez.github.io/EconPDEs.jl/dev/getting_started/#Upwinding)). "Monotone scheme" and "valid Markov generator" are the same condition seen from two sides, which is why solutions move between the two packages exactly.

## Convenience constructors

The two workhorse processes come with constructors — [`OrnsteinUhlenbeck`](@ref) and [`CoxIngersollRoss`](@ref) — that choose the grid automatically:

```@example univariate
X = OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.03, σ = 0.01)   # dx = -κ (x - xbar) dt + σ dZ
Xcir = CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 0.07)  # dx = -κ (x - xbar) dt + σ √x dZ
nothing # hide
```

By default the grid spans the `p` and `1 - p` quantiles of the stationary distribution with `length` points; pass `length`, `xmin`, `xmax`, or `pow` (grid-spacing power) to override. The default `p = 1e-10` is deliberately extreme: reflecting boundaries distort solutions near the edges of the grid (see the [expectations tutorial](expectations.md)), and a wide grid pushes that distortion where the process never goes — which matters especially when computing [tail indices](tail_index.md).

To combine these building blocks — independent products, regime switching, correlated states — see [Multivariate processes](multivariate.md).

```@example univariate
@assert abs(sum(stationary_distribution(Z)) - 1) <= 1e-12 # hide
nothing # hide
```
