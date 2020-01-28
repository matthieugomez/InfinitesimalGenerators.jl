
#========================================================================================

Compute the principal eigenvector and eigenvalue of 𝔸
By definition, it is the one associated with a positive eigenvector.
In particular, it must be real.

B = -𝔸 is a Z matrix (all off diagonal are negative). Therefore, there exists a positive s such that sI + A has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for sI + A, which is real, and the correspongind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, which means that I can look for the eigenvalue with largest real part 



If, moreover, B, is a M-matrix, then all its eigenvalues have positive real part. Therefore, all the eigenvalues of A have negative real part. Therefore, the eigenvalue with largest real part is also the eigenvalue with smallest magnitude.

========================================================================================#
function principal_eigenvalue(𝔸::Matrix; which = :SM, eigenvector = :right, r0 = ones(size(𝔸, 1)))
    l, η, r = nothing, nothing, nothing
    if eigenvector ∈ (:left, :both)
        e = eigen(𝔸')
        λs = e.values
        vs = e.vectors
        if which == :SM
            i0 = argmin(abs.(λs))
        elseif which == :LR
            i0 = argmax(real.(λs))
        end
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
        e = eigen(𝔸)
        λs = e.values
        vs = e.vectors
        if which == :SM
            i0 = argmin(abs.(λs))
        elseif which == :LR
            i0 = argmax(real.(λs))
        end
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


function principal_eigenvalue(𝔸::AbstractMatrix; which = :SM, eigenvector = :right, r0 = ones(size(𝔸, 1)))
    l, η, r = nothing, nothing, nothing
    if which == :SM
        if eigenvector ∈ (:left, :both)
            vals, vecs = Arpack.eigs(adjoint(𝔸); nev = 1, which = :SM)
            η = vals[1]
            l = vecs[:, 1]
        end
        if eigenvector ∈ (:right, :both)
            vals, vecs = Arpack.eigs(𝔸; v0 = r0, nev = 1, which = :SM)
            η = vals[1]
            r = vecs[:, 1]
        end
    elseif which == :LR
        # Arpack LR tends to fail if the LR is close to zero, which is the typical case when computing tail index
        # Arpack SM is much faster, but (i) it does not always give the right eigenvector (either because LR ≠ SM (happens when the eigenvalue is very positive) (ii) even when it gives the right eigenvalue, it can return a complex eigenvector
        if eigenvector ∈ (:left, :both)
            vals, vecs, info = KrylovKit.eigsolve(adjoint(𝔸), r0, 1, :LR, maxiter = size(𝔸, 1))
            info.converged == 0 &&  @warn "KrylovKit did not converge"
            η = vals[1]
            l = vecs[1]
        end
        if eigenvector ∈ (:right, :both)
            vals, vecs, info = KrylovKit.eigsolve(𝔸, 1, :LR, maxiter = size(𝔸, 1))
            info.converged == 0 &&  @warn "KrylovKit did not converge"
            η = vals[1]
            r = vecs[1]
        end
    end
    return clean_eigenvector_left(l), clean_eigenvalue(η), clean_eigenvector_right(r)
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
clean_eigenvector_right(::Nothing) = nothing
clean_eigenvector_right(r::AbstractVector) = abs.(r)




##############################################################################
##
## Compute Distributions
##
##############################################################################

function stationary_distribution(T::AbstractMatrix)
    g, η, _ = principal_eigenvalue(T; which = :SM, eigenvector = :left)
    abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    return g
end

# Death rate δ and reinjection ψ
function stationary_distribution(T::AbstractMatrix, δ::Number, ψ::AbstractVector{<:Number})
    clean_eigenvector_left((δ * I - T') \ (δ * ψ))
end


function tail_index(μ::Number, σ::Number; δ::Number = 0)
    if σ > 0
        (1 - 2 * μ / σ^2 + sqrt((1- 2 * μ / σ^2)^2 + 8 * δ / σ^2)) / 2
    else
        δ / μ
    end
end
# f is a function that for each ξ gives an AbstractMatrix
# find_root return ζ such that the principal eigenvalue of f(ζ) is zero
function tail_index(@nospecialize(f::Function); xatol = 1e-2, verbose = false, r0 = ones(size(f(1.0), 1)), kwargs...)
    ζ, r = nothing, nothing
    g = ξ -> begin
       out = principal_eigenvalue(f(ξ); which = :LR, r0 = r0)
       copyto!(r0, out[3])
       verbose && @show (:LR, ξ, out[2])
       return out[2]
    end
    ζ = fzero(g, (1e-2, 10.0); xatol = xatol, kwargs...)
    return ζ
end

function cgf_longrun(f::Function; which = :LR, eigenvector = :right, r0 = Ones(size(f(1), 1)))
    ξ -> principal_eigenvalue(f(ξ); which = which, eigenvector = eigenvector, r0 = r0)
end
##############################################################################
##
## Feynman Kac
##
##############################################################################

"""
With direction = :backward
Solve the PDE backward in time
u(x, t[end]) = ψ(x)
0 = u_t + 𝔸u_t - V(x, t)u +  f(x, t)


With direction = :forward
Solve the PDE forward in time
u(x, t[1]) = ψ(x)
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

