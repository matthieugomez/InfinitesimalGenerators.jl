[![Build Status](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/InfinitesimalGenerators.jl)


# Infinitesimal Generator
For a diffusive process `dx = μx(x)dt + σx(dZ_t)`
- `generator(x, μx, σx)` returns the infinitesimal generator for the process 
- `stationary_distribution(x, μx, σx)` returns the stationary distribution of the process 
- `feynman_kac_forward(x, μx, σx; t, ψ, f, V)`	returns `E[∫0^t e^{-∫0s V(x_τ)dτ}f(x_s)ds + e^{-∫0t V(x_τ)dτ}ψ(x_t)|x_0 = x]` 

# Infinitesimal Generator for Tilded Process
For a diffusive process `dx = μx(x)dt + σx(dZ_t)` and an associated geometric functional `dM/M = μM(x)dt + σM(dZ_t)`
- `generator(x, μx, σx, μM, σM)` returns the infinitesimal generator of the process `x` tilded by the geometric functional `M` 
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the hansen-scheinkman decomposition of the geometric functional `M`
- `feynman_kac_forward(x, μx, σx, μM, σM; t, ψ)` returns  `E[M_tψ(x_t)|x_0 = x]`
