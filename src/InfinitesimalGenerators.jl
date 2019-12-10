module InfinitesimalGenerators

using LinearAlgebra, Arpack, KrylovKit, Roots, Distributions, DiffEqDiffTools, FillArrays

include("utils.jl")
include("generators.jl")
include("feynman_kac.jl")

@deprecate stationary_distribution(x, μx, σx, args...) stationary_distribution(MarkovProcess(x, μx, σx))
@deprecate cgf_longrun(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0, kwargs...) cgf_longrun(MultiplicativeFunctional(MarkovProcess(x, μx, σx), μM, σM, ρ, δ); kwargs...)
@deprecate hansen_scheinkman(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0, kwargs...) hansen_scheinkman(MultiplicativeFunctional(MarkovProcess(x, μx, σx), μM, σM, ρ, δ); kwargs...)
@deprecate tail_index(x, μx, σx, μM, σM; ρ = 0.0, δ = 0.0,  kwargs...) tail_index(MultiplicativeFunctional(MarkovProcess(x, μx, σx), μM, σM, ρ, δ); kwargs...)

export 
MarkovProcess,
OrnsteinUhlenbeck,
CoxIngersollRoss,
MultiplicativeFunctional,
generator,
stationary_distribution,
feynman_kac,
hansen_scheinkman,
cgf_longrun,
tail_index



end