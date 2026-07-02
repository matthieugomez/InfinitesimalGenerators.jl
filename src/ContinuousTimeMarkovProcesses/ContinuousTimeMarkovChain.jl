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
        # validate generator (transition-rate) matrix
        size(Q, 1) == size(Q, 2) || throw(DimensionMismatch("generator (transition-rate) matrix must be square"))
        size(Q, 1) > 0 || throw(ArgumentError("generator (transition-rate) matrix cannot be empty"))
        QT = typeof(float(real(zero(eltype(Q)))))
        tol = sqrt(eps(QT)) * max(one(QT), maximum(abs.(Q)))
        maximum(abs.(sum(Q, dims = 2))) <= tol ||
            throw(ArgumentError("rows of generator (transition-rate) matrix must sum to zero"))
        all(i == j || Q[i, j] >= -tol for i in axes(Q, 1), j in axes(Q, 2)) ||
            throw(ArgumentError("off-diagonal transition rates must be nonnegative"))
        length(states) == size(Q, 1) ||
            throw(DimensionMismatch("state space and generator (transition-rate) matrix should have the same length"))
        new{TS, TQ}(states, Q)
    end
end

ContinuousTimeMarkovChain(Q::AbstractMatrix{<:Real}) =
    ContinuousTimeMarkovChain(1:size(Q, 1), Q)

state_space(X::ContinuousTimeMarkovChain) = (X.states,)

generator(X::ContinuousTimeMarkovChain) = X.Q
