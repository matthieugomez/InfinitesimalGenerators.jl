[![Build status](https://github.com/matthieugomez/InfinitesimalGenerators.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/InfinitesimalGenerators.jl/actions)

# Markov Processes
- `X = DiffusionProcess(x, μ, σ)` creates the Markov Process corresponding to a diffusive process with drift `μ(x)` and volatility `σ(x)` with reflecting boundaries.
- `generator(X)` returns its associated transition matrix
- `stationary_distribution(X)` returns its stationary distribution

# Additive Functionals
- `M = AdditiveFunctional(x, μm, σm)` creates the Additive Functional with drift  `μm(x)` and volatility `σm(x)`
- `generator(M)` creates the function `ξ -> T(ξ)` returning the tilted transition matrix (i.e. infinitesimal generator of `f -> E[e^{ξm}f(x_t)|x_0=x]`)
- `cgf(m)` returns the long run scaled CGF of `m` 
- `tail_index(m)` returns the tail index of the stationary distribution of `e^m`

## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
