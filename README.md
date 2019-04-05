[![Build Status](https://travis-ci.org/matthieugomez/ContinuousTimeMarkovOperators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/ContinuousTimeMarkovOperators.jl)

A set of tools to work with Diffusion Processes:

# Infinitesimal Generator
For a diffusive process `dx = μx(x)dt + σx(dZ_t)`
- `generator(x, μx, σx)` returns the infinetisimal generator for the process 
- `stationary_distribution(x, μx, σx)` returns the stationary distribution of the process 
- ```julia
	feynman_kac_forward(x, μx, σx; 
			t = range(0, 100, step = 1/12), ψ = ones(length(x)), f = zeros(length(x)), V = zeros(length(x)))
	```	 
	returns `E[∫0^t e^{-∫0s V(x_τ)dτ}f(x_s)ds + e^{-∫0t V(x_τ)dτ}ψ(x_t)|x_0 = x]` 

# Infinitesimal Generator for Tilded Process
For a diffusive process `dx = μx(x)dt + σx(dZ_t)` and an associated geometric functional `dM/M = μM(x)dt + σM(dZ_t)`
- `generator(x, μx, σx, μM, σM)` returns the infinitesimal generator of the process `x` tilded by the geometric functional `M` 
- `hansen_scheinkman_decomposition(x, μx, σx, μM, σM)` returns the hansen-scheinkman decomposition of the geometric functional `M`
- ```julia
	feynman_kac_forward(x, μx, σx, μM, σM; 
			t = range(0, 100, step = 1/12), ψ = ones(length(x)))
	``` 
	returns  `E[M_tψ(x_t)|x_0 = x]`
