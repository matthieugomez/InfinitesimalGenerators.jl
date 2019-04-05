[![Build Status](https://travis-ci.org/matthieugomez/ContinuousTimeMarkovOperators.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/ContinuousTimeMarkovOperators.jl)

A set of tools to work with Diffusion Processes:

# Infinitesimal Generator
For a diffusive process `dx = μx(x)dt + σx(dZ_t)`
- `generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)` returns the infinetisimal generator for the process 
- `stationary_distribution(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)` returns the stationary distribution of the process 
- `feynman_kac_forward(x, μx, σx; t::AbstractVector = range(0, 100, step = 1/12), ψ::AbstractVector, f::AbstractVector = zeros(length(x)), V::AbstractVector = zeros(length(x)))` returns the function `E[∫t^T e^{-∫ts V(x_τ)dτ}f(x_s)ds + e^{-∫tT V(x_τ)dτ}ψ(x_t)|x_0 = x]` 

# Infinitesimal Generator for Tilded Process
For a geometric functional `dM/M = μM(x)dt + σM(dZ_t)`
- `generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)` returns the infinitesimal generator of the process `x` tilded by the geometric functional `M` 
- `hansen_scheinkman_decomposition(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)` returns the hansen-scheinkman decomposition of the geometric functional `M`
- `feynman_kac_forward(x, μx, σx, μM, σM; t::AbstractVector = range(0, 100, step = 1/12), ψ = ones(length(x)))` returns the function `E[M_tψ(X_t) |X_t = x]`
