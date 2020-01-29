module InfinitesimalGenerators

using LinearAlgebra, Arpack, KrylovKit, Roots, Distributions, FiniteDiff, FillArrays



include("utils.jl")



abstract type MarkovProcess end

function stationary_distribution(X::MarkovProcess; kwargs...)
    stationary_distribution(generator(X); kwargs...)
end

function feynman_kac(X::MarkovProcess; kwargs...)
    feynman_kac(generator(X); kwargs...)
end

abstract type AdditiveFunctional end

function cgf(M::AdditiveFunctional; kwargs...)
    cgf(generator(M); kwargs...)
end

function tail_index(M::AdditiveFunctional; kwargs...)
    tail_index(generator(M); kwargs...)
end

function feynman_kac(M::AdditiveFunctional; kwargs...)
    ξ -> feynman_kac(generator(M)(ξ); kwargs...)
end

include("diffusions.jl")


export 
MarkovProcess,
AdditiveFunctional,
DiffusionProcess,
AdditiveFunctionalDiffusion,
generator,
stationary_distribution,
feynman_kac,
cgf,
tail_index
end