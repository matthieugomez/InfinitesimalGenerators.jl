module InfinitesimalGenerators

using Distributions: Normal, Gamma, quantile
using FillArrays: Ones, Zeros
using LinearAlgebra: Adjoint, Diagonal, Transpose, Tridiagonal, I, Symmetric, diag, eigmin, factorize, kron, ldiv!
using SparseArrays: blockdiag, findnz, sparse, spzeros
using Roots: fzero
using BlockBandedMatrices: BandedBlockBandedMatrix, Block

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
end
