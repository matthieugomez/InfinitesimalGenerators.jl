[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


## Markov Process
For a diffusive process 
	<img src="img/dx.png" height ="30%" width = "30%">
- `generator(x, μx, σx)` returns the infinitesimal generator 𝔸, defined as <br> <img src="img/generator.png" height ="50%" width = "50%">

- `stationary_distribution(x, μx, σx)` returns the stationary distribution of `x`, defined as the left principal eigenvector of the infinitesimal generator  𝔸 <br> <img src="img/stationary.png" height ="35%" width = "35%">
- `feynman_kac_forward(x, μx, σx; t, ψ, f, V)`	returns <img src="img/feynman_kac.png" height ="60%" width = "60%">

## Multiplicative Functional
For an associated multiplicative functional
<img src="img/dM.png" height ="40%" width = "40%">
- `generator(x, μx, σx, μM, σM)` returns the infinitesimal generator tilted by `M` <br> <img src="img/generator_tilted.png" height ="70%" width = "70%">
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the Hansen-Scheinkman decomposition of `M`.
- `feynman_kac_forward(x, μx, σx, μM, σM; t, ψ)` returns  <img src="img/feynman_kac_tilded.png" height ="22%" width = "22%">
- `impulse_response(x, μx, σx, μM, σM; t, σD)` returns  `σD(x) * (σM + σE[M_T | X_0 = x])`.

