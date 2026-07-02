module InfinitesimalGenerators

using Distributions: Normal, Gamma, quantile
using FillArrays: Ones, Zeros
using LinearAlgebra: Diagonal, Tridiagonal, I, Symmetric, eigmin, factorize, kron
using SparseArrays: blockdiag, findnz, sparse, spzeros
using Roots: fzero
using BlockBandedMatrices: BandedBlockBandedMatrix, Block


abstract type MarkovProcess end
abstract type UnivariateMarkovProcess <: MarkovProcess end
abstract type MultivariateMarkovProcess <: MarkovProcess end
# This type should define generator(), which returns a transition matrix 𝕋 such that
# 𝕋f = lim_{t→0} E[f(x_t)|x_0=x]/t
Base.size(X::MarkovProcess) = (length(state_space(X)),)
Base.length(X::MarkovProcess) = prod(size(X))
include("derivatives.jl")
include("MarkovProcesses/MarkovChain.jl")
include("MarkovProcesses/MarkovChainCompositions.jl")
include("MarkovProcesses/DiffusionProcess.jl")
include("MarkovProcesses/MultivariateDiffusionProcess.jl")


include("operators/principal_eigenvalue.jl")
include("operators/stationary_distribution.jl")
include("operators/feynman_kac.jl")
include("AdditiveFunctional.jl")




export 
MarkovProcess,
UnivariateMarkovProcess,
MultivariateMarkovProcess,
generator,
state_space,
stationary_distribution,
feynman_kac,
DiffusionProcess,
MultivariateDiffusionProcess,
MarkovChain,
ProductProcess,
SwitchingProcess,
OrnsteinUhlenbeck,
CoxIngersollRoss,
AdditiveFunctional,
cgf,
tail_index,
AdditiveFunctionalDiffusion,
FirstDerivative,
SecondDerivative,
jointoperator
end
