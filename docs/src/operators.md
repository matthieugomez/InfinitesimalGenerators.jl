# Operators

Once a process is represented by its generator matrix, everything else is linear algebra. This page lists the operators that apply to any `ContinuousTimeMarkovProcess{N}` — univariate or multivariate — with pointers to the tutorials that derive each one by hand.

Throughout, take a simple process:

```@example operators
using InfinitesimalGenerators

X = OrnsteinUhlenbeck(; xbar = 1.0, κ = 0.1, σ = 0.05)
xs = only(state_space(X))
nothing # hide
```

## The generator itself

`generator(X)` returns the generator (transition-rate) matrix ``\mathbb{T}`` of the discretized process — the operator ``f \mapsto \lim_{t \downarrow 0} E[f(x_t) \mid x_0 = x]/t`` acting on the flattened state space. Rows sum to zero and off-diagonals are non-negative, so the matrix can be used directly for anything the package does not provide: exponentials for transition probabilities, resolvents for present values, transposes for densities.

```@example operators
𝕋 = generator(X)
```

## Stationary distributions

`stationary_distribution(X)` solves the Kolmogorov forward equation ``\mathbb{T}' g = 0`` and returns the probability mass at each grid point, shaped like `size(X)` — see the [distribution dynamics tutorial](distributions.md) for the computation by hand:

```@example operators
g = stationary_distribution(X)
sum(g .* xs)   # ≈ xbar
```

The keyword form `stationary_distribution(X; δ = δ, ψ = ψ)` computes the stationary distribution when agents die at rate ``\delta`` and are reborn with distribution ``\psi`` — the resolvent ``(\delta I - \mathbb{T}')^{-1} \delta \psi`` — the relevant object in perpetual-youth and firm-entry models:

```@example operators
ψ = zeros(length(xs)); ψ[1] = 1.0    # everyone reborn at the bottom
gδ = stationary_distribution(X; δ = 0.05, ψ = ψ)
sum(gδ .* xs)
```

## Expectations: Feynman–Kac

`feynman_kac(X, ts; f, ψ, v, direction)` computes conditional expectations of the general form

```math
E\left[\int_0^T e^{-\int_0^t v(x_s) ds} f(x_t) \, dt + e^{-\int_0^T v(x_s) ds} \psi(x_T) \,\Big|\, x_0\right]
```

by implicit Euler time stepping — forecasts (`ψ`), flow payoffs and present values (`f`), state-dependent discounting or hazard rates (`v`). See the [expectations tutorial](expectations.md) for the computation by hand and the exact semantics of each keyword:

```@example operators
u = feynman_kac(X, 0:0.1:10; ψ = collect(xs))   # E[x_T | x_t], one column per date
maximum(abs, u[:, 1] - (1.0 .+ exp(-0.1 * 10.0) .* (xs .- 1.0)))
```

## Additive functionals: long-run CGFs and tail indices

An `AdditiveFunctionalDiffusion(X, μm, σm; ρ)` represents a cumulative quantity ``dm_t = \mu_m(x_t) dt + \sigma_m(x_t) dZ^m_t`` driven by the state — think ``m = \log w`` for a size ``w`` growing at a state-dependent rate. `cgf(m)` returns its long-run scaled cumulant generating function ``\xi \mapsto \lim_{t\to\infty} \log E[e^{\xi m_t}]/t``, computed as the principal eigenvalue of the tilted generator, and `tail_index(m; δ)` returns the Pareto exponent of the stationary distribution of ``e^m`` under death rate ``\delta``. See the [tail indices tutorial](tail_index.md):

```@example operators
m = AdditiveFunctionalDiffusion(X, collect(xs .- 1.0), 0.1 .* ones(length(xs)))
tail_index(m; δ = 0.05)
```

## Finite differences

The finite-difference operators the package is built on are also exported. In one dimension, the grid is a vector and the function values a vector of the same length; in multiple dimensions, the grid is a tuple of axis vectors such as `(xs, ys)`, dimensions are selected by number, and the function values are an array of matching size:

```@example operators
f = sin.(xs)
FirstDerivative(xs, f; direction = :forward, bc = (0, 0))
FirstDerivative(xs, f; direction = :backward, bc = (0, 0))
SecondDerivative(xs, f; bc = (0, 0))

ys = range(0, 2, length = 40)
F = [x * y for x in xs, y in ys]
grid = (xs, ys)
FirstDerivative(grid, F, 1; direction = :forward)
SecondDerivative(grid, F, 1, 1)
SecondDerivative(grid, F, 1, 2; direction = :up)
nothing # hide
```

The argument `bc` is the value of the *first derivative* at each limit of the grid. It defaults to zero, the right condition for reflecting boundaries — and the hook through which HJB boundary conditions like borrowing constraints enter (see the [HJB tutorial](hjb.md)).

```@example operators
@assert abs(sum(g) - 1) <= 1e-10 # hide
@assert abs(sum(gδ) - 1) <= 1e-10 # hide
@assert maximum(abs, u[:, 1] - (1.0 .+ exp(-0.1 * 10.0) .* (xs .- 1.0))) <= 1e-2 # hide
nothing # hide
```
