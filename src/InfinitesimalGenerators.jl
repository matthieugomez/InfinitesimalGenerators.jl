module InfinitesimalGenerators

using Arpack: Arpack
using Distributions: Normal, Gamma, quantile
using FillArrays: Ones, Zeros
using KrylovKit: KrylovKit
using LinearAlgebra: Diagonal, Tridiagonal, I, factorize, ldiv!, diag
using Roots: fzero
using BlockBandedMatrices: BandedBlockBandedMatrix, Block





include("MarkovProcess.jl")
include("AdditiveFunctional.jl")
include("feynman_kac.jl")
include("principal_eigenvalue.jl")
include("derivatives.jl")
include("jointoperator.jl")




export 
MarkovProcess,
generator,
state_space,
stationary_distribution,
feynman_kac,
DiffusionProcess,
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