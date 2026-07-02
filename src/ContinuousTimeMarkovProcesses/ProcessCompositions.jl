"""
    jointoperator(operators, Q)

Combine one generator matrix per regime with a regime generator (transition-rate)
matrix `Q` into the block generator of the joint process. Most users can use
[`SwitchingProcess`](@ref) instead.
"""
function jointoperator(operators::AbstractVector{<:Tridiagonal}, Q::AbstractMatrix)
    # validate operator blocks and switching matrix
    N = length(operators)
    N > 0 || throw(ArgumentError("operators cannot be empty"))
    wn = size(operators[1], 1)
    all(size(o) == (wn, wn) for o in operators) || throw(DimensionMismatch("All operators must have the same square size"))
    size(Q, 1) == size(Q, 2) || throw(DimensionMismatch("Q must be square"))
    size(Q, 1) == N || throw(DimensionMismatch("Q must have one row and column per operator"))

    J = BandedBlockBandedMatrix(Zeros(wn * N, wn * N), fill(wn, N) ,fill(wn, N), (N-1, N-1), (1, 1))
    for i in 1:N
        for j in 1:N
            J[Block(i,j)] = Q[i, j] * I(wn) + (i == j) * operators[i]
        end
    end
    return J
end

function jointoperator(operators::AbstractVector{<:AbstractMatrix}, Q::AbstractMatrix)
    # validate operator blocks and switching matrix
    N = length(operators)
    N > 0 || throw(ArgumentError("operators cannot be empty"))
    wn = size(operators[1], 1)
    all(size(o) == (wn, wn) for o in operators) || throw(DimensionMismatch("All operators must have the same square size"))
    size(Q, 1) == size(Q, 2) || throw(DimensionMismatch("Q must be square"))
    size(Q, 1) == N || throw(DimensionMismatch("Q must have one row and column per operator"))

    Iwn = sparse(I, wn, wn)
    return blockdiag(sparse.(operators)...) + kron(sparse(Q), Iwn)
end

"""
    SwitchingProcess(Z, Xs)

Returns a Markov process whose dynamics switch across the states of Markov chain
`Z`. In state `only(state_space(Z))[i]`, the process follows `Xs[i]`.
"""
struct SwitchingProcess{N, TZ <: ContinuousTimeMarkovChain, TP <: AbstractVector{<:ContinuousTimeMarkovProcess}} <: ContinuousTimeMarkovProcess{N}
    Z::TZ
    processes::TP
    function SwitchingProcess(Z::TZ, processes::TP) where {TZ <: ContinuousTimeMarkovChain, TP <: AbstractVector{<:ContinuousTimeMarkovProcess}}
        # validate switching components
        length(processes) == length(Z) ||
            throw(DimensionMismatch("there should be one process per Markov-chain state"))
        length(processes) > 0 || throw(ArgumentError("processes cannot be empty"))
        shape = size(processes[1])
        all(size(X) == shape for X in processes) ||
            throw(DimensionMismatch("all switching processes should have the same state shape"))
        space = state_space(processes[1])
        all(state_space(X) == space for X in processes) ||
            throw(DimensionMismatch("all switching processes should use the same state space"))
        N = ndims(processes[1]) + ndims(Z)
        new{N, TZ, TP}(Z, processes)
    end
end

state_space(X::SwitchingProcess) = (state_space(X.processes[1])..., state_space(X.Z)...)

generator(X::SwitchingProcess) = jointoperator(generator.(X.processes), generator(X.Z))

"""
    ProductProcess(X, Y, ...)

Returns the independent product of Markov processes in argument order.
"""
struct ProductProcess{N, TP <: Tuple} <: ContinuousTimeMarkovProcess{N}
    processes::TP
    function ProductProcess(processes::TP) where {TP <: Tuple}
        length(processes) >= 2 || throw(ArgumentError("ProductProcess needs at least two processes"))
        all(X -> X isa ContinuousTimeMarkovProcess, processes) ||
            throw(ArgumentError("all ProductProcess arguments must be Markov processes"))
        N = 0
        for process in processes
            N += ndims(process)
        end
        return new{N, TP}(processes)
    end
end

ProductProcess(processes::ContinuousTimeMarkovProcess...) = ProductProcess(processes)

function state_space(X::ProductProcess)
    axes = ()
    for process in X.processes
        axes = (axes..., state_space(process)...)
    end
    return axes
end

function generator(X::ProductProcess)
    operators = sparse.(generator.(X.processes))
    n = length(X)
    T = promote_type(map(eltype, operators)...)
    A = spzeros(T, n, n)
    before = 1
    for (operator, process) in zip(operators, X.processes)
        size(operator) == (length(process), length(process)) ||
            throw(DimensionMismatch("generator size does not match process length"))
        after = div(n, before * length(process))
        A += kron(sparse(I, after, after), kron(operator, sparse(I, before, before)))
        before *= length(process)
    end
    return A
end
