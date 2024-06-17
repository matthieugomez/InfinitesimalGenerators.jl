function jointoperator(operators, Q::Array)
    N = length(operators)
    wn = size(operators[1], 1)
    @assert all(o isa Tridiagonal for o in operators)
    # check if all os have same size
    @assert all(size(o) == (wn, wn) for o in operators) 
    # check if the size of transition matrix is 
    # same as the number of operators 
    @assert size(Q,1) == size(Q,2) == N
    J = BandedBlockBandedMatrix(Zeros(wn * N, wn * N), fill(wn, N) ,fill(wn, N), (N-1, N-1), (1, 1))
    for i in 1:N
        for j in 1:N
            J[Block(i,j)] = Q[i, j] * I(wn) + (i == j) * operators[i]
        end
    end
    return J
end
