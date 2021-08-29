"""
Compute the principal eigenvector and eigenvalue of a linear operator 𝕋, where 𝕋 is a Metzler matrix (i.e. off-diagonal components are nonnegative) (or, equilvalently, -𝕋 is a Z-matrix).

Denote a = minimum(Diagonal(V)). We have that  𝕋 - a * I has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for aI + 𝕋, which is real, and the correspondind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, which means that one can look for the eigenvalue with largest real part 

If, moreover, -𝕋 is a M-matrix, then all its eigenvalues have positive real part. Therefore, all the eigenvalues of 𝕋 have negative real part. Therefore, the eigenvalue with largest real part is also the eigenvalue with smallest magnitude.
"""
function principal_eigenvalue(𝕋::Matrix; r0 = ones(size(T, 1)))
    e = eigen(𝕋)
    λs = e.values
    vs = e.vectors
    i0 = argmax(real.(λs))
    η = λs[i0]
    r = vs[:, i0]
    for i in 1:length(λs)
        if λs[i] ≈ λs[i0]
            if all(real.(vs[:, i]) .>= - eps()) & all(abs.(imag.(vs[:, i])) .<= eps())
                r = vs[:, i]
            end
        end
    end
    abs(imag(η)) <= eps() || @warn "Principal Eigenvalue has an imaginary part"
    maximum(abs.(imag.(r))) <= eps() || @warn "Principal Eigenvector has an imaginary part"
    real(η), abs.(r)
end

function principal_eigenvalue(𝕋; r0 = ones(size(𝕋, 1)))
    η, r = 0.0, r0
    a = minimum(diag(𝕋))
    try
        vals, vecs = Arpack.eigs(𝕋 - a * I; v0 = collect(r0), nev = 1, which = :LM)
        η = vals[1]
        r = vecs[:, 1]
    catch
        vals, vecs = KrylovKit.eigsolve(𝕋 - a * I, collect(r0), 1, :LM, maxiter = size(𝕋, 1))
        η = vals[1]
        r = vecs[1]
    end
    abs(imag(η)) <= eps() || @warn "Principal Eigenvalue has an imaginary part"
    maximum(abs.(imag.(r))) <= eps() || @warn "Principal Eigenvector has an imaginary part"
    real(η) + a, abs.(r)
end

