[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


For a MarkovProcess with arithmetic drift `μx` and arithmetic volatility `σx`
- `X = MarkovProcess(x, μx, σx)` creates a `MarkovProcess`
- `stationary_distribution(X)` returns its stationary distribution
- `feynman_kac_backward(X,  t, ψ, f, V)` returns the solution of the PDE `u_t(x, t) + 𝔸 u  - V(x, t) u + f(x, t) = 0` with `u(x, T) = ψ(x)`

For a MultiplicativeFunctional `M` with geometric drift `μM(x)` and geometric volatility `σM(x)`
- `M = MultiplicativeFunctional(X, μM, σM)` creates a `MultiplicativeFunctional` 
- `cgf_longrun(M)` returns the long run scaled CGF of `log(M)` 
- `tail_index(M)` returns the tail index of its stationary distribution



## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
