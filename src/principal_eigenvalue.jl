"""
Compute the principal eigenvector and eigenvalue of a linear operator 𝕋, where 𝕋 is a Metzler matrix (i.e. off-diagonal components are nonnegative) (or, equilvalently, -𝕋 is a Z-matrix).

In this case, there exists a positive a such that aI + 𝕋 has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for aI + 𝕋, which is real, and the correspongind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, which means that one can look for the eigenvalue with largest real part 

If, moreover, -𝕋 is a M-matrix, then all its eigenvalues have positive real part. Therefore, all the eigenvalues of 𝕋 have negative real part. Therefore, the eigenvalue with largest real part is also the eigenvalue with smallest magnitude.
"""
function principal_eigenvalue(𝕋::Matrix; eigenvector = :right, r0 = ones(size(T, 1)))
    l, η, r = nothing, nothing, nothing
    if eigenvector ∈ (:left, :both)
        e = eigen(adjoint(𝕋))
        λs = e.values
        vs = e.vectors
        i0 = argmax(real.(λs))
        η = λs[i0]
        l = vs[:, i0]
        for i in 1:length(λs)
            if λs[i] ≈ λs[i0]
                if all(real.(vs[:, i]) .>= - eps()) & all(abs.(imag.(vs[:, i])) .<= eps())
                    l = vs[:, i]
                end
            end
        end
    end
   if eigenvector ∈ (:right, :both)
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
    end
    return clean_eigenvector_left(l), clean_eigenvalue(η), clean_eigenvector_right(r)
end

function principal_eigenvalue(𝕋; eigenvector = :right, r0 = ones(size(𝕋, 1)))
    l, η, r = nothing, nothing, nothing
    V = minimum(diag(𝕋))
    if eigenvector ∈ (:left, :both)
        try
            vals, vecs = Arpack.eigs(adjoint(𝕋 - V * I); nev = 1, which = :LM)
            η = vals[1]
            l = vecs[:, 1]
        catch
            vals, vecs = KrylovKit.eigsolve(adjoint(𝕋 - V * I), collect(r0), 1, :LM, maxiter = size(𝕋, 1))
            l = vecs[1]
            η = vals[1]
        end
    end
    if eigenvector ∈ (:right, :both)
        try
            vals, vecs = Arpack.eigs(𝕋 - V * I; v0 = r0, nev = 1, which = :LM)
            η = vals[1]
            r = vecs[:, 1]
        catch
            vals, vecs = KrylovKit.eigsolve(𝕋 - V * I, collect(r0), 1, :LM, maxiter = size(𝕋, 1))
            η = vals[1]
            r = vecs[1]
        end
    end
    clean_eigenvector_left(l), clean_eigenvalue(η + V), clean_eigenvector_right(r)
end


clean_eigenvalue(η::Union{Nothing, Real}) = η
function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    return real(η)
end

clean_eigenvector_left(::Nothing) = nothing
clean_eigenvector_left(l::AbstractVector) = abs.(l) ./ sum(abs.(l))

clean_eigenvector_right(::Nothing) = nothing
clean_eigenvector_right(r::AbstractVector) = abs.(r)
