[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


# General Tools

- `X = MarkovProcess(x, μx, σx)` creates a Markov Process
- `MultiplicativeFunctional(X, μM, σM)` creates an associated Multiplicative Functional

### Create Infinitesimal Generators
- `generator(MarkovProcess(X)` returns the infinitesimal generator 𝔸 associated with a Markov process `X`: <br>
	<img src="img/dx.png" height ="25%" width = "25%">: <br> <img src="img/generator.png" height ="44%" width = "44%"> <br clear="all" />
-  `generator(M)` returns the tilted infinitesimal generator 𝔸 associated with the multiplicative functional `M`: <br>
	<img src="img/dM.png" height ="33%" width = "33%">: <br> <img src="img/generator_tilted.png" height ="60%" width = "60%"> <br clear="all" />

### Work with Infinitesimal Generators
For an infinitesimal generator 𝔸:
- `principal_eigenvalue(𝔸)` returns a the principal eigenvalue of the matrix `𝔸`, its left eigenvector, and its right eigenvector
- `feynman_kac_backward(𝔸,  t, ψ, f, V)` returns the solution of the PDE `u_t(x, t) + 𝔸 u  - V(x, t) u + f(x, t) = 0` with `u(x, T) = ψ(x)`

# Convenience Functions
In addition, the package provides the following convenience functions, obtained by applying the functions above to particular generators:
- `stationary_distribution(X)` returns the stationary distribution of `x`
- `hansen_scheinkman_decomposition(M)` returns the [Hansen-Scheinkman decomposition](https://www.nber.org/papers/w12650) of `M`
- `feynman_kac(X; t, ψ, f, V, direction = :forward)`	returns <img src="img/feynman_kac.png" height ="45%" width = "45%">
- `feynman_kac_forward(M; t, ψ, direction = :forward)` returns  <img src="img/feynman_kac_tilded.png" height ="22%" width = "15%">
- `tail_index(M)` returns the tail index of the process `M`.


## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. In contrast, InfinitesimalGenerators always assumes reflecting boundaries.
- The principal eigenvalue of infinitesimal generators is found using [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)
