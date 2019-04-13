[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


# General Tools

### Create Infinitesimal Generators
- `generator(x, μx, σx)` returns the infinitesimal generator 𝔸 associated with a Markov process: <br>
	<img src="img/dx.png" height ="25%" width = "25%">: <br> <img src="img/generator.png" height ="44%" width = "44%"> <br clear="all" />
-  `generator(x, μx, σx, μM, σM)` returns the tilted infinitesimal generator 𝔸 associated with a multiplicative functional: <br>
	<img src="img/dM.png" height ="33%" width = "33%">: <br> <img src="img/generator_tilted.png" height ="40%" width = "80%"> <br clear="all" />

### Work with Infinitesimal Generators
For an infinitesimal generator 𝔸:
- `principal_eigenvalue(𝔸)` returns a the principal eigenvalue of the matrix `𝔸`, its left eigenvector, and its right eigenvector
- `feynman_kac_backward(𝔸,  t, ψ, f, V)` returns the solution of the PDE `u_t(x, t) + 𝔸 u  - V(x, t) u + f(x, t) = 0` with `u(x, T) = ψ(x)`

# Convenience Functions
In addition, the package provides the following convenience functions, obtained by applying the functions above to particular generators:
- `stationary_distribution(x, μx, σx)` returns the stationary distribution of `x`
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the [Hansen-Scheinkman decomposition](https://www.nber.org/papers/w12650) of `M`
- `feynman_kac_forward(x, μx, σx; t, ψ, f, V)`	returns <img src="img/feynman_kac.png" height ="45%" width = "45%">
- `feynman_kac_forward(x, μx, σx, μM, σM; t, ψ)` returns  <img src="img/feynman_kac_tilded.png" height ="22%" width = "15%">
- `tail_index(x, μx, σx, μM, σM)` returns the tail index of the process `M`.


## Related Packages
- [SimpleDifferentialOperators](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) contains more general tools to define operators with different boundary counditions. This package always assumes reflecting boundaries (which is the "right" boundary condition to restrict an operator defined on the whole real line on a finite grid).
- The principal eigenvalue of infinitesimal generators is found using [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)
