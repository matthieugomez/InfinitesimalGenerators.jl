module InfinitesimalGenerators

using Distributions: Normal, Gamma, quantile
using FillArrays: Ones, Zeros
using LinearAlgebra: Diagonal, Tridiagonal, I, Symmetric, eigmin, factorize, kron, ldiv!
using SparseArrays: blockdiag, findnz, sparse, spzeros
using Roots: fzero
using BlockBandedMatrices: BandedBlockBandedMatrix, Block

"""
    ContinuousTimeMarkovProcess{N}

Abstract supertype for continuous-time Markov processes with `N` tensor-product
state-space axes.

Concrete subtypes must define `state_space(X)` and `generator(X)`.
`state_space(X)` must return a tuple of state-space axes, including for
one-dimensional processes. `size(X)` is derived from those axes, and `ndims(X)`
returns the type parameter `N`.
"""
abstract type ContinuousTimeMarkovProcess{N} end
# This type should define generator(), which returns a generator matrix 𝕋 such that
# 𝕋f = lim_{t→0} E[f(x_t)|x_0=x]/t

"""
    state_space(X::ContinuousTimeMarkovProcess)

Return the state-space axes of `X` as a tuple.

The return value is always a tuple, including for one-dimensional processes:
`state_space(X)[1]` is the first axis, and `only(state_space(X))` is the grid of
a one-dimensional process.
"""
function state_space end

"""
    size(X::ContinuousTimeMarkovProcess)

Return the tensor shape of arrays defined on the state space of `X`.

This is the length of each state-space axis. Operators such as
`stationary_distribution(X)` and `feynman_kac(X, ...)` use this shape for
process-level inputs and outputs.
"""
Base.size(X::ContinuousTimeMarkovProcess) = map(length, state_space(X))

"""
    length(X::ContinuousTimeMarkovProcess)

Return the total number of grid points or finite states in `X`, equal to
`prod(size(X))`. This is the length of flattened state vectors and the row/column
dimension of `generator(X)`.
"""
Base.length(X::ContinuousTimeMarkovProcess) = prod(size(X))

"""
    ndims(X::ContinuousTimeMarkovProcess)

Return the number of tensor-product state-space axes of `X`.
"""
Base.ndims(::ContinuousTimeMarkovProcess{N}) where {N} = N
include("derivatives.jl")
include("ContinuousTimeMarkovProcesses/ContinuousTimeMarkovChain.jl")
include("ContinuousTimeMarkovProcesses/DiffusionProcess.jl")
include("ContinuousTimeMarkovProcesses/ProcessCompositions.jl")
include("ContinuousTimeMarkovProcesses/MultivariateDiffusionProcess.jl")


include("operators/principal_eigenvalue.jl")
include("operators/stationary_distribution.jl")
include("operators/feynman_kac.jl")
include("AdditiveFunctional.jl")




export 
ContinuousTimeMarkovProcess,
generator,
state_space,
stationary_distribution,
feynman_kac,
DiffusionProcess,
MultivariateDiffusionProcess,
ContinuousTimeMarkovChain,
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
