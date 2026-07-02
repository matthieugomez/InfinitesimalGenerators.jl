[![Build status](https://github.com/matthieugomez/InfinitesimalGenerators.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/InfinitesimalGenerators.jl/actions)
[![Coverage](https://codecov.io/gh/matthieugomez/InfinitesimalGenerators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/matthieugomez/InfinitesimalGenerators.jl)

This package provides tools to work with Markov processes through their finite-dimensional state spaces and infinitesimal generator matrices.

# Installation
```julia
using Pkg
Pkg.add("InfinitesimalGenerators")
```

# Markov Processes
The package allows you to compute expectations involving Markov processes.

```julia
using InfinitesimalGenerators

# Create a diffusion process (here, the Ornstein-Uhlenbeck dx = -0.03 * x * dt + 0.01 * dZ_t)
# Note that the package assumes reflecting boundaries at the limits
x = range(-1, 1, length = 100)
μx = .- 0.03 .* x
σx = 0.01 .* ones(length(x))
X = DiffusionProcess(x, μx, σx)

# Convenience constructors are also available:
# X = OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.03, σ = 0.01)
# X = CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0)

# Return its stationary distribution
g = stationary_distribution(X)

# Return the associated generator as a matrix (i.e. the operator `f -> ∂_tE[f(x_t)|x_0=x]`)
MX = generator(X)

# Use the process to compute E[∫_0^T e^{-∫_0^t v(x_s)ds}f(x_t)dt + e^{-∫_0^T v(x_s)ds}ψ(x_T) | x_0 = x]
feynman_kac(X, range(0, 100, step = 1/12); f = zeros(size(X)), ψ = ones(size(X)), v = zeros(size(X)))

# Return the grid the process is defined on
state_space(X)
```

The convenience constructors `OrnsteinUhlenbeck` and `CoxIngersollRoss` choose the grid automatically. By default it spans the `p` and `1 - p` quantiles of the stationary distribution with `length` points; pass `length`, `xmin`, `xmax`, or `pow` (grid-spacing power) to override. A small `p` is recommended when computing the tail index of an additive functional.

Any subtype of `MarkovProcess` works with `generator`, `stationary_distribution`, and `feynman_kac` as long as it defines `generator(X)` (the flat transition matrix), `state_space(X)` (the grid), and `size(X)` (the tensor shape of the state space). The package uses `UnivariateMarkovProcess` and `MultivariateMarkovProcess` subtypes for one-dimensional and tensor-shaped processes.

Finite-state Markov chains use the same interface. `Q` is a continuous-time transition-rate matrix.

```julia
z = [:low, :high]
Q = [-0.1 0.1; 0.2 -0.2]
Z = MarkovChain(z, Q)

generator(Z)
stationary_distribution(Z)
```

Use `ProductProcess` for independent Markov processes. The product follows argument order.

```julia
Y = ProductProcess(X, Z)
generator(Y)
state_space(Y)  # (state_space(X), state_space(Z))
size(Y)
length(Y)
```

Use `MultivariateDiffusionProcess` when drift and covariance terms are already
evaluated on a tensor-product grid, for example after solving an HJB.

```julia
xs = range(-1, 1, length = 50)
ys = range(0, 2, length = 40)
grid = (; x = xs, y = ys)

μx = repeat(-0.1 .* xs, 1, length(ys))
μy = repeat(reshape(1 .- ys, 1, :), length(xs), 1)
σx = 0.2 .* ones(length(xs), length(ys))
σy = 0.3 .* ones(length(xs), length(ys))
covxy = 0.01 .* ones(length(xs), length(ys))

Y = MultivariateDiffusionProcess(grid;
    drift = (; x = μx, y = μy),
    variance = (; x = σx .^ 2, y = σy .^ 2),
    covariance = (; xy = covxy),
)

generator(Y)
size(Y)
length(Y)
stationary_distribution(Y)  # array with size(Y)
feynman_kac(Y, range(0, 10, step = 1); ψ = ones(size(Y)))
```

`generator(Y)` checks that the discretized operator is a valid Markov generator:
rows must sum to zero and off-diagonal entries must be nonnegative transition
rates. The directional cross-derivative stencil used for correlated states is
monotone only under a grid-scaled diagonal-dominance condition. For two states
`x` and `y`, a useful local rule of thumb is

```julia
abs(covxy) <= variance.x * Δy / Δx
abs(covxy) <= variance.y * Δx / Δy
```

On equally spaced grids this becomes `abs(covxy) <= min(variance.x,
variance.y)`, which is stronger than positive semidefiniteness of the covariance
matrix. If `generator(Y)` throws a negative-off-diagonal error, the most common
fix is to rescale or refine the grids so `Δx / Δy` is closer to
`sqrt(variance.x / variance.y)` in the region where the covariance is large.
Equivalently, transform states so their local volatilities are more balanced, or
reduce the covariance. Use `generator(Y; check = :warn)` or `check = false` only
when you intentionally want the raw finite-difference operator without a Markov
process interpretation.

Use `SwitchingProcess` when continuous dynamics depend on the finite-state Markov chain.

```julia
Xlow = DiffusionProcess(x, μ_low, σ_low)
Xhigh = DiffusionProcess(x, μ_high, σ_high)
Y = SwitchingProcess(Z, [Xlow, Xhigh])
```

# Additive Functionals
Given a Markov process `X`, an additive functional `m` is defined by `dm = μm(x) dt + σm(x) dZm` with `corr(dZm, dZ) = ρ`.

```julia
# Create an additive functional with drift μm and volatility σm
m = AdditiveFunctionalDiffusion(X, μm, σm; ρ = 0.0)

# Return its associated generator (i.e. the operator `f -> ∂_tE[e^{m}f(x_t)|x_0=x]`)
generator(m)

# Return the long run scaled CGF of m, i.e. ξ -> lim_{t→∞} log(E[e^{ξ m_t}])/t
cgf(m)(1.0)

# Return the tail index of the stationary distribution of e^m
tail_index(m)
```

# Finite Differences
The package also allows you to compute finite differences of a function on a grid. In one dimension, the grid is a vector `x` and the function values are a vector `f` with the same length. In multiple dimensions, the grid is a `NamedTuple` such as `(; x = xs, y = ys)`, and the function values are an array `F` with size `(length(xs), length(ys))`.

```julia
using InfinitesimalGenerators
x = range(-1, 1, length = 100)
f = sin.(x)
FirstDerivative(x, f; direction = :forward, bc = (0, 0))
FirstDerivative(x, f; direction = :backward, bc = (0, 0))
SecondDerivative(x, f; bc = (0, 0))

xs = range(-1, 1, length = 50)
ys = range(0, 2, length = 40)
grid = (; x = xs, y = ys)
F = [x * y for x in xs, y in ys]
FirstDerivative(grid, F, :x; direction = :forward)
SecondDerivative(grid, F, :x, :x)
SecondDerivative(grid, F, :x, :y; direction = :up)
```
The argument `bc` refers to the value of the first derivative at each limit of the grid. This argument defaults to zero, which is the right condition when solving problems with reflecting boundaries.

# Joint Operator
For low-level operator work, `jointoperator` combines individual generators with a transition matrix. Most users can use `SwitchingProcess` instead.

```julia
Q = [-0.1 0.1; 0.2 -0.2]  # regime transition matrix
J = jointoperator([generator(X1), generator(X2)], Q)
```

## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary conditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
