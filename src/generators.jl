
##############################################################################
##
## Markov Process
##
## must define generator(X)
## which corresponds to 𝔸f = E[df]
##
##############################################################################
abstract type MarkovProcess end

function stationary_distribution(X::MarkovProcess; kwargs...)
    stationary_distribution(generator(X); kwargs...)
end
function feynman_kac(X::MarkovProcess; kwargs...)
    feynman_kac(generator(X); kwargs...)
end

##############################################################################
##
## Multiplicative Functional
##
## must define generator(X)
## which corresponds to ξ -> 𝔸(ξ)
## where A(ξ)f = E[d(M^ξf)]

##############################################################################
abstract type MultiplicativeFunctional end

function cgf_longrun(M::MultiplicativeFunctional; kwargs...)
    cgf_longrun(generator(M); kwargs...)
end

function tail_index(M::MultiplicativeFunctional; kwargs...)
    tail_index(generator(M); kwargs...)
end

function feynman_kac(M::MultiplicativeFunctional; kwargs...)
    feynman_kac(generator(M)(1); kwargs...)
end
