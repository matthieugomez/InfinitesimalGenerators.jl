#========================================================================================

Compute the foward looking operator

v_0 * g - ∂(v1 * g) + 0.5 * ∂∂(v2 * g)

========================================================================================#
function build_operator(x, v0, v1, v2)
    Δ = EconPDEs.make_Δ(x)
    n = length(x)
    T = BandedMatrix(Zeros(n, n), (1, 1))
    build_operator!(T, Δ, v0, v1, v2)
end

function build_operator!(T, Δ, v0, v1, v2)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(T, 0.0)
    # construct matrix T. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    for i in 1:n
        if v1[i] >= 0
            T[min(i + 1, n), i] += v1[i] * invΔxp[i]
            T[i, i] -= v1[i] * invΔxp[i]
        else
            T[i, i] += v1[i] * invΔxm[i]
            T[max(i - 1, 1), i] -= v1[i] * invΔxm[i]
        end
        T[max(i - 1, 1), i] += v2[i] * invΔxm[i] * invΔx[i]
        T[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        T[min(i + 1, n), i] += v2[i] * invΔxp[i] * invΔx[i]
        # Make sure each column sums to zero. Important in some cases: for isntance, otherwise cannot find sdf decomposition in GP model
        T[i, i] += v0[i] - sum(view(T, :, i))
    end
    return T
end


#========================================================================================

Compute the principal eigenvector and eigenvalue of an operator

========================================================================================#
clean_eigenvalue(η::Union{Nothing, Real}) = η

function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    real(η)
end
clean_density(::Nothing) = nothing
clean_density(v::Vector) = abs.(v) ./ sum(abs.(v))

clean_f(v::Vector) = abs.(v)
clean_f(::Nothing) = nothing

function principal_eigenvalue(T; method = :krylov, eigenvector = :right)
    η = nothing
    if method == :krylov
        vl, η, vr = principal_eigenvalue_krylov(T; eigenvector = eigenvector)
        if η == nothing
            @warn "Krylov Methods Failed"
        end
    end
    if η == nothing
        vl, η, vr = principal_eigenvalue_BLAS(convert(Matrix{Float64}, T); eigenvector = eigenvector)
    end
    return clean_f(vl), clean_eigenvalue(η), clean_density(vr)
end

function principal_eigenvalue_krylov(T; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        vals, vecs, info = KrylovKit.eigsolve(T, 1, :LR)
        if info.converged > 0
            η = vals[1]
            vr = vecs[1]
        end
    end
    if eigenvector ∈ (:left, :both)
        vals, vecs, info = KrylovKit.eigsolve(T', 1, :LR)
        if info.converged > 0
            η = vals[1]
            vl = vecs[1]
        end
    end 
    return vl, η, vr
end


function principal_eigenvalue_BLAS(T; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        e = eigen(T)
        _, out = findmax(real.(e.values))
        η = e.values[out]
        vr = e.vectors[:, out]
    end
    if eigenvector ∈ (:left, :both)
        e = eigen(copy(T'))
        _, out = findmax(real.(e.values))
        η = e.values[out]
        vl = e.vectors[:, out]
    end 
    return vl, η, vr
end
