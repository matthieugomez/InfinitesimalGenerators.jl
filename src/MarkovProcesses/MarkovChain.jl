"""
    MarkovChain(z, Q)
    MarkovChain(Q)

Returns a continuous-time Markov chain on the discrete state space `z`, with
transition-rate matrix `Q`.
"""
struct MarkovChain{TZ <: AbstractVector, TQ <: AbstractMatrix{<:Real}} <: UnivariateMarkovProcess
    z::TZ
    Q::TQ
    function MarkovChain(z::TZ, Q::TQ) where {TZ <: AbstractVector, TQ <: AbstractMatrix{<:Real}}
        # validate transition matrix
        size(Q, 1) == size(Q, 2) || throw(DimensionMismatch("transition matrix must be square"))
        size(Q, 1) > 0 || throw(ArgumentError("transition matrix cannot be empty"))
        T = typeof(float(real(zero(eltype(Q)))))
        tol = sqrt(eps(T)) * max(one(T), maximum(abs.(Q)))
        maximum(abs.(sum(Q, dims = 2))) <= tol ||
            throw(ArgumentError("rows of transition matrix must sum to zero"))
        all(i == j || Q[i, j] >= -tol for i in axes(Q, 1), j in axes(Q, 2)) ||
            throw(ArgumentError("off-diagonal transition rates must be nonnegative"))
        length(z) == size(Q, 1) ||
            throw(DimensionMismatch("state space and transition matrix should have the same length"))
        new{TZ, TQ}(z, Q)
    end
end

MarkovChain(Q::AbstractMatrix{<:Real}) = MarkovChain(1:size(Q, 1), Q)

state_space(X::MarkovChain) = X.z

Base.size(X::MarkovChain) = (length(X.z),)

generator(X::MarkovChain) = X.Q
