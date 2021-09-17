[![Build status](https://github.com/matthieugomez/InfinitesimalGenerators.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/InfinitesimalGenerators.jl/actions)

# Markov Processes
For a Markov process `x` 
- `generator(x)` creates the transition matrix of a diffusive process with drift `μ(x)` and volatility `σ(x)` with reflecting boundaries.
- `stationary_distribution` returns the stationary distribution of `x`

# Additive Functionals
For an additive functional `m`:
- `generator(m)` creates the function `ξ -> T(ξ)` for the additive functional with drift `μm(x)` and volatility `σm(x)`
- `cgf(m)` returns the long run scaled CGF of `m` 
- `tail_index(m)` returns the tail index of the stationary distribution of `e^m`

## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
