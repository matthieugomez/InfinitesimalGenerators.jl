[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


## Markov Process
For a diffusive process `dx = μx(x)dt + σx(dZ_t)`
- `generator(x, μx, σx)` returns the infinitesimal generator of `x`.
- `stationary_distribution(x, μx, σx)` returns the stationary distribution of `x`.
- `feynman_kac_forward(x, μx, σx; t, ψ, f, V)`	returns `E[∫0^t e^{-∫0s V(x_τ)dτ}f(x_s)ds + e^{-∫0t V(x_τ)dτ}ψ(x_t)|x_0 = x]`. 

## Multiplicative Functional
For a diffusive process `dx = μx(x)dt + σx(dZ_t)` and an associated multiplicative functional `dM/M = μM(x)dt + σM(dZ_t)`
- `generator(x, μx, σx, μM, σM)` returns the infinitesimal generator of `x` tilded by `M`. 
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the Hansen-Scheinkman decomposition of `M`.
- `feynman_kac_forward(x, μx, σx, μM, σM; t, ψ)` returns  `E[M_tψ(x_t)|x_0 = x]`.
- `impulse_response(x, μx, σx, μM, σM; t, σD)` returns  `σD(x) * (σM + σE[M_T | X_0 = x])`.

