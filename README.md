[![Build status](https://github.com/matthieugomez/InfinitesimalGenerators.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/InfinitesimalGenerators.jl/actions)

This package provides a set of tools to work with Markov Processes defined on a 1-dimensional Gri



# Markov Processes
The package allows you to compute compute expectations involving Markov processes

```julia
using InfinitesimalGenerators
# Create a diffusion process (here,  the Ornstein–Uhlenbeck dx = -0.03 * x * dt + 0.01 * dZ_t)
# Note that the package assumes reflecting boundaries at the limits
x = range(-1, 1, length = 100)
μx = .- 0.03 .* x
σx = 0.01 .* ones(length(x))
X = DiffusionProcess(x, μx, σx)


# Return its stationary distribution 
g = stationary_distribution(X)

# Return the associated generator as a matrix (i.e. the operator `f -> ∂_tE[f(x_t)|x_0=x]`)
MX = generator(X)

# Use the generator to compute E[\int_0^T e^{-\int_0^t f(x_s)ds}f(x_s) +  e^{-\int_0^T v(x_s)ds}ψ(x_T) | x_0 = x]
feynman_kac(MX; t = range(0, 100, step = 1/12), f = zeros(length(x)),  ψ = ones(length(x)), v = zeros(length(x)))
```

# Additive Functionals
- `M = AdditiveFunctional(X, μm, σm)` creates, given a discretized Markov Process, the Additive Functional with drift  `μm` and volatility `σm`
- `generator(M)` returns its associated generator (i.e. the operator `f -> ∂_tE[e^{m}f(x_t)|x_0=x]`)
- `cgf(m)` returns the long run scaled CGF of `m` 
- `tail_index(m)` returns the tail index of the stationary distribution of `e^m`



# Derivative
The package also allows you to compute (lazy) first and second derivatives of a function on a grid using a finite difference schemes

```julia
using InfinitesimalGenerators
x = range(-1, 1, length = 100)
f = sin.(x)
FirstDerivative(x, f, direction = :upward, bc = (0, 0))
FirstDerivative(x, f, direction = :downward, bc = (0, 0))
SecondDerivative(x, f, bc = (0, 0))
```
The argument `bc` refers to the value of the first-derivative at each limit of the grid. This argument defaults to zero, which is the right condition when solving problems with reflecting boundaries.


## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
