[![Build status](https://github.com/matthieugomez/InfinitesimalGenerators.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/InfinitesimalGenerators.jl/actions)

# Markov Processes
- `X = DiffusionProcess(x::AbstractVector, μ::AbstractVector, σ::AbstractVector)` creates the discretized Markov Process with drift `μ` and volatility `σ`, on a grid `x` with reflecting boundaries.
- `generator(X)` returns its associated generator (i.e. the operator `f -> ∂_tE[f(x_t)|x_0=x]`)
- `stationary_distribution(X)` returns its stationary distribution (i.e. the positive vector `g` such that `g * generator(X) = 0`)

# Additive Functionals
- `M = AdditiveFunctional(X, μm, σm)` creates, given a discretized Markov Process, the Additive Functional with drift  `μm` and volatility `σm`
- `generator(M, ξ)` returns its associated generator (i.e. the operator `f -> ∂_tE[e^{ξm}f(x_t)|x_0=x]`)
- `cgf(m)` returns the long run scaled CGF of `m` 
- `tail_index(m)` returns the tail index of the stationary distribution of `e^m`

## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
