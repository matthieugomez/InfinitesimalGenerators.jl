module InfinitesimalGenerators

using LinearAlgebra, Arpack, KrylovKit, Roots, Distributions, DiffEqDiffTools, FillArrays

include("utils.jl")
include("generators.jl")
include("diffusions.jl")

@deprecate stationary_distribution(x, μx, σx, args...) stationary_distribution(MarkovDiffusion(x, μx, σx))
@deprecate cgf_longrun(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0, kwargs...) cgf_longrun(MultiplicativeFunctional(MarkovDiffusion(x, μx, σx), μM, σM; ρ = ρ, δ = δ); kwargs...)
@deprecate hansen_scheinkman(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0, kwargs...) hansen_scheinkman(MultiplicativeFunctional(MarkovDiffusion(x, μx, σx), μM, σM; ρ = ρ, δ = δ); kwargs...)
@deprecate tail_index(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0,  kwargs...) tail_index(MultiplicativeFunctional(MarkovDiffusion(x, μx, σx), μM, σM; ρ = ρ, δ = δ); kwargs...)
@deprecate feynman_kac(x, μx, σx; kwargs...) feynman_kac(MarkovDiffusion(x, μx, σx); kwargs...)
@deprecate feynman_kac(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0, kwargs...) feynman_kac(MultiplicativeFunctional(MarkovDiffusion(x, μx, σx), μM, σM; ρ = ρ, δ = δ); kwargs...)

export 
MarkovProcess,
MultiplicativeFunctional,
generator,
stationary_distribution,
feynman_kac,
cgf_longrun,
hansen_scheinkman,
tail_index,

MarkovDiffusion,
OrnsteinUhlenbeck,
CoxIngersollRoss,
MultiplicativeFunctionalDiffusion

end