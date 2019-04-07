[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


## Markov Process
For a diffusive process <img src="img/dx.png">
- `generator(x, μx, σx)` returns the infinitesimal generator of `x` <img src="img/generator.png" style="height:10px;">

- `stationary_distribution(x, μx, σx)` returns the stationary distribution of `x`.
- `feynman_kac_forward(x, μx, σx; t, ψ, f, V)`	returns <img src="img/feynman_kac.png" style="height:10px;">

## Multiplicative Functional
For an associated multiplicative functional <img src="img/dM.png" style="height:10px;">
- `generator(x, μx, σx, μM, σM)` returns the infinitesimal generator of `x` tilted by `M` <img src="img/generator_tilted.png" style="height:10px;">
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the Hansen-Scheinkman decomposition of `M`.
- `feynman_kac_forward(x, μx, σx, μM, σM; t, ψ)` returns  <img src="img/feynman_kac_tilded.png" style="height:10px;">
- `impulse_response(x, μx, σx, μM, σM; t, σD)` returns  `σD(x) * (σM + σE[M_T | X_0 = x])`.

