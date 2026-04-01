function jointoperator(operators::AbstractVector{<:Tridiagonal}, Q::AbstractMatrix)
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
