


#========================================================================================

Compute the principal eigenvector and eigenvalue of 𝔸
By definition, it is the one associated with a positive eigenvector.
In particular, it must be real.

B = -𝔸 is a Z matrix (all off diagonal are negative). Therefore, there exists a positive s such that sI + A has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for sI + A, which is real, and the correspongind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, which means that I can look for the eigenvalue with largest real part 



If, moreover, B, is a M-matrix, then all its eigenvalues have positive real part. Therefore, all the eigenvalues of A have negative real part. Therefore, the eigenvalue with largest real part is also the eigenvalue with smallest magnitude.

========================================================================================#
function principal_eigenvalue(𝔸::AbstractMatrix; which = :SM, eigenvector = :right, r0 = ones(size(𝔸, 1)))
    l, η, r = nothing, nothing, nothing
    if which == :SM
        vals, vecs = Arpack.eigs(adjoint(𝔸); nev = 1, which = :SM)
        η = vals[1]
        l = vecs[:, 1]
        if eigenvector ∈ (:right, :both)
                vals, vecs = Arpack.eigs(𝔸; v0 = r0, nev = 1, which = :SM)
                η = vals[1]
                r = vecs[:, 1]
        end
    elseif which == :LR
        # Arpack LR tends to fail if the LR is close to zero, which is the typical case if we want to compute tail index
        # Arpack SM is much faster, but it does not always give the right eigenvector (either because LR ≠ SM (happens when the eigenvalue is very positive)
        # Even when it gives the right eigenvalue, it can return a complex eigenvector
        vals, vecs, info = KrylovKit.eigsolve(adjoint(𝔸), r0, 1, :LR, maxiter = size(𝔸, 1))
        info.converged == 0 &&  @warn "KrylovKit did not converge"
        η = vals[1]
        l = vecs[1]
        if eigenvector ∈ (:right, :both)
            vals, vecs, info = KrylovKit.eigsolve(𝔸, 1, :LR, maxiter = size(𝔸, 1))
            info.converged == 0 &&  @warn "KrylovKit did not converge"
            η = vals[1]
            r = vecs[1]
        end
    end
    l = clean_eigenvector_left(l)
    return l, clean_eigenvalue(η), clean_eigenvector_right(l, r)
end

clean_eigenvalue(η::Union{Nothing, Real}) = η
function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    real(η)
end

clean_eigenvector_left(::Nothing) = nothing
clean_eigenvector_left(l::AbstractVector) = abs.(l) ./ sum(abs.(l))


# correct normalization is \int r l = 1
clean_eigenvector_right(l, ::Nothing) = nothing
clean_eigenvector_right(l, r::AbstractVector) = abs.(r) ./ sum(l .* abs.(r))






# f is a function that for each ξ gives the generator matrix
# find_root return ζ such that the principal eigenvalue of f(ζ) is zero
function find_root(f::Function; which = :SM, xatol = 1e-2, verbose = false, r0 = ones(size(f(1.0), 1)), kwargs...)
    out = 0.0
    if which == :SM
        try
            # SM is so much faster. So try if it works.
            f = ξ -> begin
                out = principal_eigenvalue(f(ξ); which = :SM, r0 = r0)
                eltype(out[3]) <: Float64 && copyto!(r0, out[3])
                verbose && @show (:SM, ξ, out[2])
                return out[2]
            end
            D = ξ -> FiniteDiff.finite_difference_derivative(f, ξ)
            out = find_zero((f, D), 1.0, Roots.Newton(); xatol = xatol, kwargs...)
            out2 = principal_eigenvalue(f(out); which = :LR, r0 = r0)[2]
            if abs(out2) > 1e-2 
                @warn "Algorithm looking for SM eigenvalue = 0 converged to ζ = $out. However, the :LR eigenvalue for this ζ is  $out2"
                throw("there is an error")
            end
        catch
            which = :LR
        end
    end
    if which == :LR
        f = ξ -> begin
            out = principal_eigenvalue(f(ξ); which = :LR, r0 = r0)
            eltype(out[3]) <: Float64 && copyto!(r0, out[3])
            verbose && @show (:LR, ξ, out[2])
            return out[2]
        end
        D = ξ -> FiniteDiff.finite_difference_derivative(f, ξ)
        try
            out = find_zero((f, D), 1.0, Roots.Newton(); xatol = xatol, kwargs...)
        catch
            out = find_zero((f, D), (1e-2, 10.0); xatol = xatol, kwargs...)
        end
    end
    return out
end

##############################################################################
##
## Feynman Kac
##
##############################################################################

"""
With direction = :backward
Solve the PDE backward in time
u(x, T) = ψ(x)
0 = u_t + 𝔸u_t - V(x, t)u +  f(x, t)


With direction = :forward
Solve the PDE forward in time
u(x, 0) = ψ(x)
u_t = 𝔸u - V(x)u + f(x)
"""
function feynman_kac(𝔸::AbstractMatrix; 
    t::AbstractVector = range(0, 100, step = 1/12), 
    ψ::AbstractVector = ones(size(𝔸, 1)), 
    f::Union{AbstractVector, AbstractMatrix} = zeros(size(𝔸, 1)), 
    V::Union{AbstractVector, AbstractMatrix} = zeros(size(𝔸, 1)),
    direction= :backward)
    if direction == :backward
        u = zeros(size(𝔸, 1), length(t))
        u[:, end] = ψ
        if isa(f, AbstractVector) && isa(V, AbstractVector)
            if isa(t, AbstractRange)
                dt = step(t)
                𝔹 = factorize(I + (Diagonal(V) - 𝔸) * dt)
                for i in (length(t)-1):(-1):1
                    ψ = ldiv!(𝔹, u[:, i+1] .+ f .* dt)
                    u[:, i] = ψ
                end
            else
                for i in (length(t)-1):(-1):1
                    dt = t[i+1] - t[i]
                    𝔹 = I + (Diagonal(V) - 𝔸) * dt
                    u[:, i] = 𝔹 \ (u[:, i+1] .+ f .* dt)
                end
            end
        elseif isa(f, AbstractMatrix) && isa(V, AbstractMatrix)
            for i in (length(t)-1):(-1):1
                dt = t[i+1] - t[i]
                𝔹 = I + (Diagonal(view(V, :, i)) - 𝔸) * dt
                u[:, i] = 𝔹 \ (u[:, i+1] .+ f[:, i] .* dt)
            end
        else
            error("f and V must be Vectors or Matrices")
        end
        return u
    elseif direction == :forward
        u = feynman_kac(𝔸; t = - reverse(t), ψ = ψ, f = f, V = V, direction = :backward)
        return u[:,end:-1:1]
    else
        error("Direction must be :backward or :forward")
    end
end

