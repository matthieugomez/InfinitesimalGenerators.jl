[![Build status](https://github.com/matthieugomez/InfinitesimalGenerators.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/InfinitesimalGenerators.jl/actions)

This package provides a set of tools to work with Markov Processes defined on a 1-dimensional grid.

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

# Use the generator to compute E[∫_0^T e^{-∫_0^t v(x_s)ds}f(x_t)dt + e^{-∫_0^T v(x_s)ds}ψ(x_T) | x_0 = x]
feynman_kac(MX, range(0, 100, step = 1/12); f = zeros(length(x)), ψ = ones(length(x)), v = zeros(length(x)))
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

# Derivatives
The package also allows you to compute (lazy) first and second derivatives of a function on a grid using finite difference schemes.

```julia
using InfinitesimalGenerators
x = range(-1, 1, length = 100)
f = sin.(x)
FirstDerivative(x, f; direction = :forward, bc = (0, 0))
FirstDerivative(x, f; direction = :backward, bc = (0, 0))
SecondDerivative(x, f; bc = (0, 0))
```
The argument `bc` refers to the value of the first derivative at each limit of the grid. This argument defaults to zero, which is the right condition when solving problems with reflecting boundaries.

# Joint Operator
For coupled Markov processes switching between `N` regimes, `jointoperator` combines the individual generators with a transition matrix.

```julia
Q = [-0.1 0.1; 0.2 -0.2]  # regime transition matrix
J = jointoperator([generator(X1), generator(X2)], Q)
```

## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary conditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
