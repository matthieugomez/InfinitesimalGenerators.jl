[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


For a Markov process defined by a matrix `T`  where `T` is the operator such that `Tf = E[df]`
- `stationary_distribution(T)` returns its stationary distribution
- `feynman_kac_backward(T,  t, ψ, f, V)` returns the solution of the PDE `u_t(x, t) + T u  - V(x, t) u + f(x, t) = 0` with `u(x, T) = ψ(x)`

For an additive functional `m` defined by a function `ξ -> T(ξ)` where `T` is the operator such that `T f= E[d(e^(ξm)f)]` 
- `cgf(f)` returns the long run scaled CGF of `m` 
- `tail_index(f)` returns the tail index of the stationary distribution of `e^m`


Moreover, 
- `generator(DiffusionProcess(x, μ, σ))` creates the transition matrix of a diffusive process with drift `μ(x)` and volatility `σ(x)` (with reflecting boundaries)
- `generator(AdditiveFunctional(DiffusionProcess(x, μ, σ), μm, σm)` creates the function ``ξ -> T(ξ)` for the additive functional with drift `μm(x)` and volatility `σm(x)`

## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
