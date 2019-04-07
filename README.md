[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


## Markov Process
For a diffusive process 
	<img src="img/dx.png" height ="30%" width = "30%">
- `generator(x, μx, σx)` returns the infinitesimal generator 𝔸: <br> <img src="img/generator.png" height ="50%" width = "50%">

- `stationary_distribution(x, μx, σx)` returns the stationary distribution of `x`, i.e. the left principal eigenvector of  𝔸 <br> <img src="img/stationary.png" height ="35%" width = "35%">
- `feynman_kac_forward(x, μx, σx; t, ψ, f, V)`	returns <img src="img/feynman_kac.png" height ="60%" width = "60%">

## Multiplicative Functional
For an associated multiplicative functional
<img src="img/dM.png" height ="40%" width = "40%">
- `generator(x, μx, σx, μM, σM)` returns the tilted infinitesimal generator 𝔸: <br> <img src="img/generator_tilted.png" height ="70%" width = "70%">
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the [Hansen-Scheinkman decomposition](https://www.nber.org/papers/w12650) of `M`, i.e. the principal eigenvalue/eigenvectors of 𝔸.
- `feynman_kac_forward(x, μx, σx, μM, σM; t, ψ)` returns  <img src="img/feynman_kac_tilded.png" height ="22%" width = "22%">
- `impulse_response(x, μx, σx, μM, σM; t, σD)` returns  `σD(x) * (σM + σE[M_T | X_0 = x])`.

## Related Packages
- This package relies on [BandedMatrices.jl](https://github.com/JuliaMatrices/BandedMatrices.jl) to represent infinitesimal generators as matrices, as well as [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl) to solve for the principal eigenvalue/eigenvector.
- This package is related to [DiffEqOperators.jl](https://github.com/JuliaDiffEq/DiffEqOperators.jl), which contains more general tools to solve differential equations.