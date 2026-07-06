module InfinitesimalGenerators

using Distributions: Normal, Gamma, quantile
using FillArrays: Ones
using LinearAlgebra: LinearAlgebra, Adjoint, Diagonal, LAPACKException, SingularException, Transpose, Tridiagonal, I, Symmetric, diag, eigmin, factorize, kron, ldiv!
using PrecompileTools: @compile_workload
using SparseArrays: blockdiag, findnz, sparse, spzeros
using Roots: fzero

include("derivatives.jl")
include("matrix_operators/check_generator.jl")
include("matrix_operators/principal_eigenvalue.jl")
include("matrix_operators/stationary_distribution.jl")
include("matrix_operators/feynman_kac.jl")
include("ContinuousTimeMarkovProcesses.jl")
include("ContinuousTimeMarkovProcesses/ContinuousTimeMarkovChain.jl")
include("ContinuousTimeMarkovProcesses/DiffusionProcess.jl")
include("ContinuousTimeMarkovProcesses/ProcessCompositions.jl")
include("ContinuousTimeMarkovProcesses/MultivariateDiffusionProcess.jl")
include("AdditiveFunctional.jl")




export
ContinuousTimeMarkovProcess,
generator,
check_generator,
state_space,
stationary_distribution,
feynman_kac,
principal_eigenvalue,
DiffusionProcess,
MultivariateDiffusionProcess,
ContinuousTimeMarkovChain,
ProductProcess,
SwitchingProcess,
OrnsteinUhlenbeck,
CoxIngersollRoss,
AdditiveFunctional,
tilted_generator,
cgf,
cgf_eigenvector,
tail_index,
AdditiveFunctionalDiffusion,
FirstDerivative,
SecondDerivative,
jointoperator

# compile the main entry points on tiny grids so the first real call is fast
@compile_workload begin
    # univariate diffusion (Tridiagonal generator)
    X = OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 0.1, length = 10)
    stationary_distribution(X)
    feynman_kac(X, 0.0:0.5:1.0; ψ = only(state_space(X)) .^ 2)

    # chain, switching, and product processes (sparse generators)
    Z = ContinuousTimeMarkovChain([0.0, 1.0], [-0.5 0.5; 0.5 -0.5])
    stationary_distribution(SwitchingProcess(Z, [X, X]))
    stationary_distribution(ProductProcess(X, Z))

    # multivariate diffusion (sparse generator)
    xs = collect(range(-1.0, 1.0, length = 5))
    Xmv = MultivariateDiffusionProcess((; x = xs, y = xs);
        drift = (; x = -0.1 .* xs .* ones(1, 5), y = -0.1 .* ones(5) .* xs'),
        variance = (; x = 0.01, y = 0.01))
    stationary_distribution(Xmv)

    # derivative operators (lazy arrays: collect compiles the indexing paths)
    f = only(state_space(X)) .^ 2
    collect(FirstDerivative(only(state_space(X)), f; direction = :backward))
    collect(SecondDerivative(only(state_space(X)), f))
    F = xs .* xs'
    collect(FirstDerivative(Xmv.grid, F, :x))
    collect(SecondDerivative(Xmv.grid, F, 1, 1))
    collect(SecondDerivative(Xmv.grid, F, 1, 2))
end

end
