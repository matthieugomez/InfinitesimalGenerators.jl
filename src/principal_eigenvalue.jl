"""
Compute the principal eigenvector and eigenvalue of a linear operator 𝕋, where 𝕋 is a Metzler matrix (i.e. off-diagonal components are nonnegative)

Denote a = -minimum(Diagonal(V)), which implies 𝕋 + a * I has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for aI + 𝕋, which is real, and the correspondind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, and so this also correspoinds to the eigenvalue with largest real part of 𝕋, which happens to be real.
Denote η(𝕋) the eigenvalue with largest real part of a matrix and ρ(𝕋) the eigenvalue with largest modulus. We have

    η(𝕋) = ρ(𝕋 + a * I) - a

Moreover, the associated eigenvector is real and strictly positive.


Note that, when 𝕋 is generator, its rows sum to zero. This means that eigenvalue with largest real part is 
    η(𝕋) = 0
In other words, all eigenvalues of 𝕋 have real part <= 0. This means that 𝕋 is a singular M matrix.
(another proof is to say that, for any s, sI - 𝕋 is a non-singular M-Matrix for any s> 0, since there exists x = e such that (sI - 𝕋) * x > 0). 

"""
function principal_eigenvalue(𝕋; r0 = ones(size(𝕋, 1)))
    η, r = 0.0, r0
    a = - minimum(diag(𝕋))
    try
        vals, vecs = Arpack.eigs(𝕋 + a * I; v0 = collect(r0), nev = 1, which = :LM)
        η = vals[1]
        r = vecs[:, 1]
    catch
        vals, vecs = KrylovKit.eigsolve(𝕋 + a * I, collect(r0), 1, :LM; maxiter = size(𝕋, 1))
        η = vals[1]
        r = vecs[1]
    end
    abs(imag(η)) <= eps() || @warn "Principal Eigenvalue has an imaginary part"
    maximum(abs.(imag.(r))) <= eps() || @warn "Principal Eigenvector has an imaginary part"
    real(η) - a, abs.(r)
end

