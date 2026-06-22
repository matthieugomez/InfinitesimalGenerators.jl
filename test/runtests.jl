using InfinitesimalGenerators, Test, Statistics, LinearAlgebra, Expokit

# Shared Ornstein–Uhlenbeck process used across several test sets
xbar = 0.0
κ = 0.1
σ = 0.02
X = OrnsteinUhlenbeck(; xbar = xbar, κ = κ, σ = σ, length = 1000)


@testset "Feynman-Kac" begin
    ψ = X.x.^2
    ts = range(0, stop = 100, step = 1/10)
    u = feynman_kac(generator(X), ts; ψ = ψ, direction = :forward)
    @test maximum(abs, u[:, 50] .- expmv(ts[50], generator(X), ψ)) <= 1e-3
    @test maximum(abs, u[:, 200] .- expmv(ts[200], generator(X), ψ)) <= 1e-3
    @test maximum(abs, u[:, end] .- expmv(ts[end], generator(X), ψ)) <= 1e-5
    @test maximum(abs, feynman_kac(generator(X), ts; ψ = ψ, direction = :forward) .- feynman_kac(generator(X), ts; ψ = ψ, direction = :forward)) <= 1e-5

    # scalar generator with time-varying f and v
    T_scalar = zeros(1, 1)
    ts_scalar = 0:1:2
    f_scalar = reshape([1.0, 2.0, 3.0], 1, :)
    v_scalar = reshape([0.0, 1.0, 3.0], 1, :)
    u_scalar = feynman_kac(T_scalar, ts_scalar; f = f_scalar, ψ = [1.0], v = v_scalar, direction = :forward)
    u_expected = zeros(1, length(ts_scalar))
    u_expected[:, 1] .= 1.0
    for i in 1:(length(ts_scalar) - 1)
        dt = ts_scalar[i + 1] - ts_scalar[i]
        B = I + Diagonal(v_scalar[:, i + 1]) * dt
        u_expected[:, i + 1] = B \ (u_expected[:, i] .+ f_scalar[:, i + 1] .* dt)
    end
    @test u_scalar ≈ u_expected
end


@testset "feynman_kac element type" begin
    # output element type should follow the inputs, not be hard-coded to Float64
    Xf = OrnsteinUhlenbeck(; κ = 0.1, σ = 0.02, length = 50)
    M = generator(Xf)
    𝕋 = Tridiagonal(Float32.(M.dl), Float32.(M.d), Float32.(M.du))
    ts = Float32.(0:0.1:10)   # concrete Float32 grid (a Float32 range has Float64 eltype on Julia 1.6)
    u = feynman_kac(𝕋, ts; ψ = Float32.(Xf.x .^ 2), direction = :forward)
    @test eltype(u) == Float32
    @test all(isfinite, u)
end


@testset "stationary_distribution" begin
    g_discounted = stationary_distribution(X; δ = 1e-2)
    @test all(isfinite, g_discounted)
    @test sum(g_discounted) ≈ 1.0 atol = 1e-12
    @test all(g_discounted .>= 0.0)
    @test_throws ArgumentError stationary_distribution(X; δ = 1e-2, ψ = zeros(length(X.x)))
    @test_throws ArgumentError DiffusionProcess([0.0], [0.0], [1.0])
    @test_throws ArgumentError DiffusionProcess([0.0, 0.0, 1.0], zeros(3), ones(3))
    @test_throws ArgumentError DiffusionProcess([0.0, 1.0, 0.5], zeros(3), ones(3))
end


@testset "Multiplicative functional: cgf" begin
    # dM/M = x dt
    m = AdditiveFunctionalDiffusion(X, X.x, zeros(length(X.x)))
    η, r = cgf(m)(1)
    @test η ≈ xbar + 0.5 * σ^2 / κ^2 atol = 1e-2
    r_analytic = exp.(X.x ./ κ)
    @test norm(r ./ sum(r) .- r_analytic ./ sum(r_analytic)) <= 2 * 1e-3
    ts = range(0, stop = 200, step = 1/10)
    u = feynman_kac(generator(m), ts; direction = :forward, ψ = ones(size(generator(m), 1)))
    @test log.(stationary_distribution(X)' * u[:, end]) ./ ts[end] ≈ η atol = 1e-2
end


@testset "tail_index and speed" begin
    μm = -0.06
    m = AdditiveFunctionalDiffusion(X, X.x .+ μm, zeros(length(X.x)))
    ζ = tail_index(m)
    @test μm * ζ + 0.5 * ζ^2 * (σ^2 / κ^2) ≈ 0.0 atol = 1e-2
    η, r = cgf(m; eigenvector = :right)(ζ)
    η, l = cgf(m; eigenvector = :left)(ζ)
    f = exp.(ζ .* X.x ./ κ)
    @test norm(f ./ sum(f) .- r ./ sum(r)) <= 1e-2
    ψ_reaching = r .* l ./ sum(r .* l)
    speed = sum(ψ_reaching .* m.μm)
    @test speed ≈ μm + ζ * (σ^2 / κ^2) atol = 1e-2
end


@testset "left/right eigenvectors (correlated)" begin
    m = AdditiveFunctionalDiffusion(X, X.x, 0.01 * ones(length(X.x)); ρ = 1)
    η, r = cgf(m; eigenvector = :right)(1)
    η, l = cgf(m; eigenvector = :left)(1)
    ψ_tilde = stationary_distribution(DiffusionProcess(X.x, X.μx .+ m.ρ .* m.σm .* X.σx, X.σx))
    @test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3
end


@testset "Multiplicative functional (ρ = 0)" begin
    μm = -0.01
    σm = 0.1
    m = AdditiveFunctionalDiffusion(X, μm .+ X.x, σm .* ones(length(X.x)))
    ζ = tail_index(m)
    ζ_analytic = 2 * (-μm) / (σm^2 + (σ / κ)^2)
    @test ζ ≈ ζ_analytic atol = 1e-2
    η, r = cgf(m; eigenvector = :right)(ζ)
    η, l = cgf(m; eigenvector = :left)(ζ)
    @test η ≈ 0.0 atol = 1e-4
    ψ = stationary_distribution(X)
    @test (r .* ψ) ./ sum(r .* ψ) ≈ l rtol = 1e-3
end


@testset "Twisted process stationary distribution" begin
    # the modified process μ + σ² ∂ ln(r) has stationary distribution r²ψ
    Xl = OrnsteinUhlenbeck(; κ = κ, σ = σ, length = 1000)
    μm = -0.01
    σm = 0.1
    m = AdditiveFunctionalDiffusion(Xl, μm .+ Xl.x .- 0.02, σm .* ones(length(Xl.x)))
    ψ = stationary_distribution(Xl)
    ζ = tail_index(m)
    η, r = cgf(m; eigenvector = :right)(ζ)
    η, l = cgf(m; eigenvector = :left)(ζ)
    ψ_cond = stationary_distribution(DiffusionProcess(Xl.x, Xl.μx .+ Xl.σx.^2 .* (InfinitesimalGenerators.∂(Xl) * log.(r)), Xl.σx))
    @test (r.^2 .* ψ) ./ sum(r.^2 .* ψ) ≈ ψ_cond rtol = 1e-1
end


@testset "Multiplicative functional (ρ = 1)" begin
    Xl = OrnsteinUhlenbeck(; κ = κ, σ = σ, length = 1000)
    μm = -0.01
    σm = 0.1
    m0 = AdditiveFunctionalDiffusion(Xl, μm .+ Xl.x .- 0.02, σm .* ones(length(Xl.x)))
    m = AdditiveFunctionalDiffusion(Xl, m0.μm, m0.σm; ρ = 1.0)
    ζ = tail_index(m)
    η, r = cgf(m; eigenvector = :right)(ζ)
    η, l = cgf(m; eigenvector = :left)(ζ)
    @test η ≈ 0.0 atol = 1e-3
    ψ_tilde = stationary_distribution(DiffusionProcess(Xl.x, Xl.μx .+ ζ .* m.σm .* m.ρ .* Xl.σx, Xl.σx))
    @test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3
end


@testset "Cox-Ingersoll-Ross" begin
    gbar = 0.03
    σ_cir = 0.01
    Xc = CoxIngersollRoss(xbar = gbar, κ = κ, σ = σ_cir)
    m = AdditiveFunctionalDiffusion(Xc, Xc.x, zeros(length(Xc.x)))
    η_analytic = gbar * κ^2 / σ_cir^2 * (1 - sqrt(1 - 2 * σ_cir^2 / κ^2))
    @test cgf(m)(1.0)[1] ≈ η_analytic rtol = 1e-2
end


@testset "FirstDerivative and SecondDerivative" begin
    x = range(0.0, stop = 1.0, length = 1000)
    y = x.^2
    dy = FirstDerivative(x, y; direction = :forward)
    @test length(dy) == length(x)
    @test dy[500] ≈ 2 * x[500] atol = 1e-2
    dy_back = FirstDerivative(x, y; direction = :backward)
    @test dy_back[500] ≈ 2 * x[500] atol = 1e-2
    # boundary conditions: forward derivative at last point returns bc
    @test dy[end] == 0.0
    # backward derivative at first point returns bc
    @test dy_back[1] == 0.0

    d2y = SecondDerivative(x, y)
    @test length(d2y) == length(x)
    @test d2y[500] ≈ 2.0 atol = 1e-2
end


@testset "∂ with zero-drift node" begin
    # grid placed so that x = 0 (and hence μx = 0) is exactly an interior node
    x = range(-1.0, 1.0, length = 101)
    μx = -0.1 .* x
    σx = 0.02 .* ones(length(x))
    Xz = DiffusionProcess(x, μx, σx)
    @test μx[51] == 0                       # zero drift at the middle node
    D = InfinitesimalGenerators.∂(Xz)
    @test all(isfinite, D)                  # Diagonal(μx) \ generator would give NaN here
    df = D * (x .^ 2)
    @test df[51] ≈ 0.0 atol = 1e-12         # central difference: d(x²)/dx = 0 at x = 0
    @test df[30] ≈ 2 * x[30] atol = 5e-2    # upwind away from the zero-drift node
end


@testset "jointoperator" begin
    X1 = OrnsteinUhlenbeck(; κ = 0.1, σ = 0.02, length = 50)
    X2 = OrnsteinUhlenbeck(; κ = 0.2, σ = 0.03, length = 50)
    T1 = generator(X1)
    T2 = generator(X2)
    Q = [-0.1 0.1; 0.2 -0.2]
    J = jointoperator([T1, T2], Q)
    @test size(J) == (100, 100)
    # rows should sum to zero (generator property)
    @test maximum(abs.(sum(Matrix(J), dims = 2))) < 1e-10
    # principal eigenvalue of a generator should be zero
    Jdense = Matrix(J)
    η, r = InfinitesimalGenerators.principal_eigenvalue(Jdense)
    @test η ≈ 0.0 atol = 1e-8
    @test all(r .> 0)
    @test_throws DimensionMismatch jointoperator([T1], Q)
    @test_throws DimensionMismatch jointoperator([T1, generator(OrnsteinUhlenbeck(; κ = 0.1, σ = 0.02, length = 60))], Q)
    @test_throws DimensionMismatch jointoperator([T1, T2], [-0.1 0.1 0.0; 0.2 -0.2 0.0])
end
