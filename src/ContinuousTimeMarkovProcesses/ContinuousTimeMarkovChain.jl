"""
    ContinuousTimeMarkovChain(states, Q)
    ContinuousTimeMarkovChain(Q)

Returns a finite-state continuous-time Markov chain on the discrete state space
`states`, with generator (transition-rate) matrix `Q`.

The matrix `Q` is a generator: rows must sum to zero and off-diagonal entries
must be nonnegative transition rates. It is not a discrete-time transition
probability matrix.
"""
struct ContinuousTimeMarkovChain{TS <: AbstractVector, TQ <: AbstractMatrix{<:Real}} <: ContinuousTimeMarkovProcess{1}
    states::TS
    Q::TQ
    function ContinuousTimeMarkovChain(states::TS, Q::TQ) where {TS <: AbstractVector, TQ <: AbstractMatrix{<:Real}}
        check_generator(Q)
        length(states) == size(Q, 1) ||
            throw(DimensionMismatch("state space and generator (transition-rate) matrix should have the same length"))
        new{TS, TQ}(states, Q)
    end
end

ContinuousTimeMarkovChain(Q::AbstractMatrix{<:Real}) =
    ContinuousTimeMarkovChain(1:size(Q, 1), Q)

state_space(X::ContinuousTimeMarkovChain) = (X.states,)

generator(X::ContinuousTimeMarkovChain) = X.Q
