using InfinitesimalGenerators, Test, Statistics, LinearAlgebra, SparseArrays, Expokit
using BlockBandedMatrices: BandedBlockBandedMatrix

# Shared Ornstein–Uhlenbeck process used across several test sets
xbar = 0.0
κ = 0.1
σ = 0.02
X = OrnsteinUhlenbeck(; xbar = xbar, κ = κ, σ = σ, length = 1000)

struct CustomMarkovProcess <: ContinuousTimeMarkovProcess{1} end

InfinitesimalGenerators.state_space(::CustomMarkovProcess) = (1:2,)
InfinitesimalGenerators.generator(::CustomMarkovProcess) = [-1.0 1.0; 2.0 -2.0]


@testset "Feynman-Kac" begin
    ψ = X.x.^2
    ts = range(0, stop = 100, step = 1/10)
    u = feynman_kac(generator(X), ts; ψ = ψ, direction = :forward)
    @test feynman_kac(X, ts; ψ = ψ, direction = :forward) ≈ u
    @test maximum(abs, u[:, 50] .- expmv(ts[50], generator(X), ψ)) <= 1e-3
    @test maximum(abs, u[:, 200] .- expmv(ts[200], generator(X), ψ)) <= 1e-3
    @test maximum(abs, u[:, end] .- expmv(ts[end], generator(X), ψ)) <= 1e-5
    @test maximum(abs, feynman_kac(generator(X), ts; ψ = ψ, direction = :forward) .- feynman_kac(generator(X), ts; ψ = ψ, direction = :forward)) <= 1e-5

    𝔸_sparse = sparse([-1.0 1.0; 1.0 -1.0])
    u_sparse = feynman_kac(𝔸_sparse, 0.0:1.0:2.0; f = zeros(2), ψ = ones(2))
    @test u_sparse ≈ ones(2, 3)

    # scalar generator with time-varying f and v
    𝔸_scalar = zeros(1, 1)
    ts_scalar = 0:1:2
    f_scalar = reshape([1.0, 2.0, 3.0], 1, :)
    v_scalar = reshape([0.0, 1.0, 3.0], 1, :)
    u_scalar = feynman_kac(𝔸_scalar, ts_scalar; f = f_scalar, ψ = [1.0], v = v_scalar, direction = :forward)
    u_expected = zeros(1, length(ts_scalar))
    u_expected[:, 1] .= 1.0
    for i in 1:(length(ts_scalar) - 1)
        dt = ts_scalar[i + 1] - ts_scalar[i]
        B = I + Diagonal(v_scalar[:, i + 1]) * dt
        u_expected[:, i + 1] = B \ (u_expected[:, i] .+ f_scalar[:, i + 1] .* dt)
    end
    @test u_scalar ≈ u_expected

    u_matrix_f = feynman_kac(𝔸_scalar, ts_scalar; f = f_scalar, ψ = [0.0], v = [0.0])
    @test u_matrix_f ≈ feynman_kac(𝔸_scalar, ts_scalar; f = f_scalar, ψ = [0.0], v = zeros(1, length(ts_scalar)))
    u_matrix_v = feynman_kac(𝔸_scalar, ts_scalar; f = [1.0], ψ = [0.0], v = v_scalar)
    @test u_matrix_v ≈ feynman_kac(𝔸_scalar, ts_scalar; f = ones(1, length(ts_scalar)), ψ = [0.0], v = v_scalar)
end


@testset "feynman_kac element type" begin
    # output element type should follow the inputs, not be hard-coded to Float64
    Xf = OrnsteinUhlenbeck(; κ = 0.1, σ = 0.02, length = 50)
    𝔸64 = generator(Xf)
    𝔸 = Tridiagonal(Float32.(𝔸64.dl), Float32.(𝔸64.d), Float32.(𝔸64.du))
    ts = Float32.(0:0.1:10)   # concrete Float32 grid (a Float32 range has Float64 eltype on Julia 1.6)
    u = feynman_kac(𝔸, ts; ψ = Float32.(Xf.x .^ 2), direction = :forward)
    @test eltype(u) == Float32
    @test all(isfinite, u)
end


@testset "stationary_distribution" begin
    g_discounted = stationary_distribution(X; δ = 1e-2)
    @test size(X) == (length(X.x),)
    @test length(X) == length(X.x)
    @test all(isfinite, g_discounted)
    @test sum(g_discounted) ≈ 1.0 atol = 1e-12
    @test all(g_discounted .>= 0.0)
    @test_throws ArgumentError stationary_distribution(X; δ = 1e-2, rebirth = zeros(length(X)))
    ψ_rebirth = zeros(length(X)); ψ_rebirth[1] = 1.0
    @test stationary_distribution(X; δ = 1e-2, rebirth = ψ_rebirth) ≈
          stationary_distribution(X; δ = 1e-2, ψ = ψ_rebirth)   # deprecated keyword still works
    @test_throws ArgumentError stationary_distribution(X; δ = 1e-2, rebirth = ψ_rebirth, ψ = ψ_rebirth)
    @test_throws ArgumentError DiffusionProcess([0.0], [0.0], [1.0])
    @test_throws ArgumentError DiffusionProcess([0.0, 0.0, 1.0], zeros(3), ones(3))
    @test_throws ArgumentError DiffusionProcess([0.0, 1.0, 0.5], zeros(3), ones(3))
    X_positive = OrnsteinUhlenbeck(; xbar = 1.0, κ = 0.1, σ = 0.02, p = 0.01, length = 10, pow = 2)
    @test all(only(state_space(X_positive)) .> 0.0)
    @test issorted(only(state_space(X_positive)); lt = <)

    C = CustomMarkovProcess()
    @test size(C) == (2,)
    @test length(C) == 2
    @test stationary_distribution(C) ≈ [2 / 3, 1 / 3]
end


@testset "ContinuousTimeMarkovChain, ProductProcess, and SwitchingProcess" begin
    Q = [-0.1 0.1; 0.2 -0.2]
    Z = ContinuousTimeMarkovChain([:low, :high], Q)
    @test Z.states == [:low, :high]
    @test state_space(Z) == ([:low, :high],)
    @test size(Z) == (2,)
    @test length(Z) == 2
    @test Z isa ContinuousTimeMarkovProcess{1}
    @test ndims(Z) == 1
    @test Z isa ContinuousTimeMarkovChain
    @test generator(Z) == Q
    @test stationary_distribution(Z) ≈ [2 / 3, 1 / 3]
    @test ContinuousTimeMarkovChain(Q).Q == Q
    @test_logs (:warn, r"row sum") ContinuousTimeMarkovChain([-0.1 0.2; 0.1 -0.2])
    @test_logs (:warn, r"Metzler") ContinuousTimeMarkovChain([-0.1 0.1; -0.2 0.2])

    x_small = range(-0.2, stop = 0.2, length = 40)
    Xbase = DiffusionProcess(x_small, -0.1 .* x_small, 0.02 .* ones(length(x_small)))
    @test state_space(Xbase) == (x_small,)
    @test Xbase isa ContinuousTimeMarkovProcess{1}
    @test ndims(Xbase) == 1
    Y = ProductProcess(Xbase, Z)
    G = generator(Y)
    @test state_space(Y) == (x_small, [:low, :high])
    @test length.(state_space(Y)) == (length(x_small), 2)
    @test size(Y) == (length(x_small), 2)
    @test length(Y) == 2 * length(x_small)
    @test Y isa ContinuousTimeMarkovProcess{2}
    @test ndims(Y) == 2
    @test G isa SparseMatrixCSC
    @test size(G) == (2 * length(x_small), 2 * length(x_small))
    @test maximum(abs.(sum(Matrix(G), dims = 2))) < 1e-10
    @test Matrix(G) ≈ Matrix(generator(SwitchingProcess(Z, fill(Xbase, length(Z)))))

    ψx = stationary_distribution(Xbase)
    πz = stationary_distribution(Z)
    ψY = stationary_distribution(Y)
    @test size(ψY) == size(Y)
    @test ψY ≈ ψx * πz' rtol = 1e-5

    Xlow = DiffusionProcess(x_small, -0.1 .* x_small, 0.02 .* ones(length(x_small)))
    Xhigh = DiffusionProcess(x_small, 0.03 .- 0.2 .* x_small, 0.03 .* ones(length(x_small)))
    Yswitch = SwitchingProcess(Z, [Xlow, Xhigh])
    @test state_space(Yswitch) == (x_small, [:low, :high])
    @test length.(state_space(Yswitch)) == (length(x_small), 2)
    @test size(Yswitch) == (length(x_small), 2)
    @test length(Yswitch) == 2 * length(x_small)
    @test Yswitch isa ContinuousTimeMarkovProcess{2}
    @test ndims(Yswitch) == 2
    Gswitch = generator(Yswitch)
    @test Gswitch isa BandedBlockBandedMatrix
    @test Matrix(Gswitch) ≈ Matrix(jointoperator([generator(Xlow), generator(Xhigh)], Q))
    @test Matrix(jointoperator([generator(Xlow), generator(Xhigh)], Q)) ≈
          Matrix(jointoperator(sparse.([generator(Xlow), generator(Xhigh)]), Q))
    @test_throws DimensionMismatch SwitchingProcess(Z, [Xlow])
    @test_throws DimensionMismatch SwitchingProcess(Z, [Xlow, DiffusionProcess(range(-0.2, stop = 0.2, length = 41), zeros(41), ones(41))])
end

@testset "check_generator" begin
    # every generator the package builds passes, whatever its matrix type
    A = generator(X)
    @test check_generator(A) === A                                       # Tridiagonal, chainable
    Zc = ContinuousTimeMarkovChain([0.5, 1.5], [-0.1 0.1; 0.2 -0.2])
    @test check_generator(generator(ProductProcess(X, Zc))) isa SparseMatrixCSC
    @test check_generator(generator(SwitchingProcess(Zc, [X, X]))) isa BandedBlockBandedMatrix
    @test check_generator(Matrix(generator(Zc))) isa Matrix

    # violations warn (with the size of the violation); structural problems still throw
    @test_logs (:warn, r"row sum") check_generator([-0.1 0.2; 0.1 -0.2])
    @test_logs (:warn, r"Metzler") check_generator([0.1 -0.1; 0.2 -0.2])
    @test_throws ArgumentError check_generator(zeros(2, 1))            # non-square
    @test_logs (:warn, r"Metzler") check_generator(sparse([1.0 -1.0; -1.0 1.0]))
    # tolerance override
    @test_logs (:warn, r"row sum") check_generator([-1e-6 1e-6; 2e-6 -1e-6]; atol = 1e-9)
    @test check_generator([-1e-6 1e-6; 2e-6 -1e-6]; atol = 1e-5) isa Matrix

    # principal_eigenvalue checks Metzlerity, its Perron–Frobenius hypothesis
    @test_logs (:warn, r"Metzler") principal_eigenvalue([0.1 -0.1; 0.2 -0.2])
    @test principal_eigenvalue(A')[1] ≈ 0.0 atol = 1e-8   # transposed generators pass
end

@testset "SwitchingProcess with a diffusion modulator" begin
    # modulation by an autonomous diffusion == a MultivariateDiffusionProcess with independent innovations
    as = collect(range(0.0, 1.0, length = 15))
    zs = collect(range(-1.0, 1.0, length = 20))
    Zd = DiffusionProcess(zs, -0.2 .* zs, 0.05 .* ones(length(zs)))
    S = SwitchingProcess(Zd, [DiffusionProcess(as, (0.1 + 0.05 * z) .- 0.1 .* as, 0.02 .* ones(length(as))) for z in zs])
    @test size(S) == (length(as), length(zs))
    Xmv = MultivariateDiffusionProcess((; a = as, z = zs);
        drift = (; a = [0.1 + 0.05 * z - 0.1 * a for a in as, z in zs],
                   z = repeat(reshape(-0.2 .* zs, 1, :), length(as), 1)),
        variance = (; a = 0.02^2, z = 0.05^2))
    @test Matrix(generator(S)) ≈ Matrix(generator(Xmv)) atol = 1e-10
    @test stationary_distribution(S) ≈ stationary_distribution(Xmv) atol = 1e-8
end

@testset "MultivariateDiffusionProcess" begin
    xs = collect(range(-1.0, 1.0, length = 8))
    ys = collect(range(-2.0, 2.0, length = 7))
    μx = -0.1 .* xs
    μy = -0.2 .* ys
    σx = 0.3 .* ones(length(xs))
    σy = 0.4 .* ones(length(ys))

    Xx = DiffusionProcess(xs, μx, σx)
    Xy = DiffusionProcess(ys, μy, σy)

    grid = (; x = xs, y = ys)
    drift = (; x = repeat(μx, 1, length(ys)),
               y = repeat(reshape(μy, 1, :), length(xs), 1))
    variance = (; x = repeat(σx .^ 2, 1, length(ys)),
                  y = repeat(reshape(σy .^ 2, 1, :), length(xs), 1))

    Xxy = MultivariateDiffusionProcess(grid; drift = drift, variance = variance)
    Gxy = generator(Xxy)
    @test state_space(Xxy) == (xs, ys)
    @test size(Xxy) == (length(xs), length(ys))
    @test length(Xxy) == length(xs) * length(ys)
    @test Xxy isa ContinuousTimeMarkovProcess{2}
    @test ndims(Xxy) == 2
    Gxy_expected = kron(Matrix(I, length(ys), length(ys)), Matrix(generator(Xx))) +
                   kron(Matrix(generator(Xy)), Matrix(I, length(xs), length(xs)))
    @test Matrix(Gxy) ≈ Gxy_expected

    Pxy = ProductProcess(Xx, Xy)
    @test state_space(Pxy) == (xs, ys)
    @test size(Pxy) == size(Xxy)
    @test Matrix(generator(Pxy)) ≈ Gxy_expected
    Zxy = ContinuousTimeMarkovChain([:low, :high], [-0.1 0.1; 0.2 -0.2])
    Pxyz = ProductProcess(Xxy, Zxy)
    @test state_space(Pxyz) == (xs, ys, [:low, :high])
    @test size(Pxyz) == (length(xs), length(ys), 2)
    Sxyz = SwitchingProcess(Zxy, [Xxy, Xxy])
    @test state_space(Sxyz) == (xs, ys, [:low, :high])
    @test size(Sxyz) == (length(xs), length(ys), 2)

    ψxy = stationary_distribution(Xxy)
    ψ_expected = stationary_distribution(Xx) * stationary_distribution(Xy)'
    @test size(ψxy) == size(Xxy)
    @test sum(ψxy) ≈ 1.0 atol = 1e-12
    @test ψxy ≈ ψ_expected rtol = 1e-8 atol = 1e-10
    @test stationary_distribution(Pxy) ≈ ψ_expected rtol = 1e-8 atol = 1e-10

    ts_xy = 0.0:0.25:1.0
    ψ_terminal = ones(size(Xxy))
    u_xy = feynman_kac(Xxy, ts_xy; ψ = ψ_terminal, direction = :forward)
    u_xy_flat = feynman_kac(Gxy, ts_xy; ψ = vec(ψ_terminal), direction = :forward)
    @test size(u_xy) == (size(Xxy)..., length(ts_xy))
    @test reshape(u_xy, length(Xxy), length(ts_xy)) ≈ u_xy_flat
    u_xy_backward = feynman_kac(Xxy, ts_xy; ψ = ψ_terminal)
    u_xy_backward_flat = feynman_kac(Gxy, ts_xy; ψ = vec(ψ_terminal))
    @test reshape(u_xy_backward, length(Xxy), length(ts_xy)) ≈ u_xy_backward_flat

    covxy = 0.05 .* ones(length(xs), length(ys))
    Xcorr = MultivariateDiffusionProcess(grid; drift = drift, variance = variance,
        covariance = (; xy = covxy))
    Gcorr = generator(Xcorr)
    @test maximum(abs.(sum(Gcorr, dims = 2))) < 1e-10
    @test minimum([Gcorr[i, j] for i in axes(Gcorr, 1), j in axes(Gcorr, 2) if i != j]) >= -1e-12
    Xcorr_negative = MultivariateDiffusionProcess(grid; drift = drift, variance = variance,
        covariance = (; xy = -covxy))
    Gcorr_negative = generator(Xcorr_negative)
    @test maximum(abs.(sum(Gcorr_negative, dims = 2))) < 1e-10
    @test minimum([Gcorr_negative[i, j] for i in axes(Gcorr_negative, 1), j in axes(Gcorr_negative, 2) if i != j]) >= -1e-12

    bad_drift = (; x = zeros(length(xs), length(ys)),
                   y = zeros(length(xs), length(ys)))
    bad_variance = (; x = 0.01 .* ones(length(xs), length(ys)),
                      y = ones(length(xs), length(ys)))
    bad_covxy = 0.09 .* ones(length(xs), length(ys))
    Xnonmonotone = MultivariateDiffusionProcess(grid; drift = bad_drift,
        variance = bad_variance, covariance = (; xy = bad_covxy))
    err = try
        generator(Xnonmonotone)
        nothing
    catch err
        err
    end
    @test err isa ArgumentError
    @test occursin("negative off-diagonal", sprint(showerror, err))
    @test occursin("Δx / Δy", sprint(showerror, err))
    Gnonmonotone = generator(Xnonmonotone; check = false)
    Gnonmonotone_warn = @test_logs (:warn, r"negative off-diagonal") generator(Xnonmonotone; check = :warn)
    @test Gnonmonotone_warn == Gnonmonotone
    @test minimum([Gnonmonotone[i, j] for i in axes(Gnonmonotone, 1), j in axes(Gnonmonotone, 2) if i != j]) < 0
    @test_throws ArgumentError generator(Xnonmonotone; check = :invalid)

    @test_throws ArgumentError MultivariateDiffusionProcess(grid; drift = drift,
        variance = variance, covariance = (; x = ones(length(xs), length(ys))))
    @test_throws ArgumentError MultivariateDiffusionProcess(grid; drift = drift,
        variance = variance, covariance = (; xy = 2.0 .* ones(length(xs), length(ys))))
    @test_throws ArgumentError MultivariateDiffusionProcess(Dict(:x => xs, :y => ys);
        drift = drift, variance = variance)
    @test_throws ArgumentError MultivariateDiffusionProcess(grid;
        drift = Dict(:x => drift.x, :y => drift.y), variance = variance)
    @test_throws ArgumentError MultivariateDiffusionProcess(grid;
        drift = drift, variance = Dict(:x => variance.x, :y => variance.y))
    @test_throws ArgumentError MultivariateDiffusionProcess(grid, (; bad = drift.x), variance, nothing)
    ambiguous_grid = (; a = [0.0, 1.0], bc = [0.0, 1.0], ab = [0.0, 1.0], c = [0.0, 1.0])
    ambiguous_drift = (; a = 0.0, bc = 0.0, ab = 0.0, c = 0.0)
    ambiguous_variance = (; a = 1.0, bc = 1.0, ab = 1.0, c = 1.0)
    @test_throws ArgumentError MultivariateDiffusionProcess(ambiguous_grid;
        drift = ambiguous_drift, variance = ambiguous_variance)

    μ_boundary = [-1.0; zeros(length(xs) - 2); 1.0]
    σ_boundary = 0.1 .* ones(length(xs))
    X_boundary = DiffusionProcess(xs, μ_boundary, σ_boundary)
    X_boundary_nd = MultivariateDiffusionProcess((; x = xs);
        drift = (; x = μ_boundary), variance = (; x = σ_boundary .^ 2))
    @test Matrix(generator(X_boundary)) ≈ Matrix(generator(X_boundary_nd))
    G_boundary = Matrix(generator(X_boundary))
    @test minimum([G_boundary[i, j] for i in axes(G_boundary, 1), j in axes(G_boundary, 2) if i != j]) >= -1e-12
end


@testset "Multiplicative functional: cgf" begin
    # dM/M = x dt
    m = AdditiveFunctional(X, X.x, zeros(length(X.x)))
    η, r = cgf_eigenvector(m, 1)
    @test cgf(m, 1) ≈ η
    @test η ≈ xbar + 0.5 * σ^2 / κ^2 atol = 1e-2
    r_analytic = exp.(X.x ./ κ)
    @test norm(r ./ sum(r) .- r_analytic ./ sum(r_analytic)) <= 2 * 1e-3
    ts = range(0, stop = 200, step = 1/10)
    u = feynman_kac(tilted_generator(m, 1), ts; direction = :forward, ψ = ones(length(m.X)))
    @test log.(stationary_distribution(X)' * u[:, end]) ./ ts[end] ≈ η atol = 1e-2
    # deprecated closure form still works
    η_dep, r_dep = cgf(m)(1)
    @test η_dep ≈ η
    @test_throws ArgumentError cgf(m; eigenvector = :middle)(1)
end


@testset "tail_index and speed" begin
    μm = -0.06
    m = AdditiveFunctional(X, X.x .+ μm, zeros(length(X.x)))
    ζ = tail_index(m)
    @test μm * ζ + 0.5 * ζ^2 * (σ^2 / κ^2) ≈ 0.0 atol = 1e-2
    η, r = cgf_eigenvector(m, ζ)
    η, l = cgf_eigenvector(m, ζ, :left)
    f = exp.(ζ .* X.x ./ κ)
    @test norm(f ./ sum(f) .- r ./ sum(r)) <= 1e-2
    ψ_reaching = r .* l ./ sum(r .* l)
    speed = sum(ψ_reaching .* m.drift)
    @test speed ≈ μm + ζ * (σ^2 / κ^2) atol = 1e-2
end


@testset "left/right eigenvectors (correlated)" begin
    m = AdditiveFunctional(X, X.x, 0.01 * ones(length(X.x)); ρ = 1)
    η, r = cgf_eigenvector(m, 1)
    η, l = cgf_eigenvector(m, 1, :left)
    ψ_tilde = stationary_distribution(DiffusionProcess(X.x, X.μx .+ m.covariance, X.σx))
    @test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3
end


@testset "Multiplicative functional (ρ = 0)" begin
    μm = -0.01
    σm = 0.1
    m = AdditiveFunctional(X, μm .+ X.x, σm .* ones(length(X.x)))
    ζ = tail_index(m)
    ζ_analytic = 2 * (-μm) / (σm^2 + (σ / κ)^2)
    @test ζ ≈ ζ_analytic atol = 1e-2
    η, r = cgf_eigenvector(m, ζ)
    η, l = cgf_eigenvector(m, ζ, :left)
    @test η ≈ 0.0 atol = 1e-4
    ψ = stationary_distribution(X)
    @test (r .* ψ) ./ sum(r .* ψ) ≈ l rtol = 1e-3
end


@testset "Twisted process stationary distribution" begin
    # the modified process μ + σ² ∂ ln(r) has stationary distribution r²ψ
    Xl = OrnsteinUhlenbeck(; κ = κ, σ = σ, length = 1000)
    μm = -0.01
    σm = 0.1
    m = AdditiveFunctional(Xl, μm .+ Xl.x .- 0.02, σm .* ones(length(Xl.x)))
    ψ = stationary_distribution(Xl)
    ζ = tail_index(m)
    η, r = cgf_eigenvector(m, ζ)
    ψ_cond = stationary_distribution(DiffusionProcess(Xl.x, Xl.μx .+ Xl.σx.^2 .* (InfinitesimalGenerators.∂(Xl) * log.(r)), Xl.σx))
    @test (r.^2 .* ψ) ./ sum(r.^2 .* ψ) ≈ ψ_cond rtol = 1e-1
end


@testset "Multiplicative functional (ρ = 1)" begin
    Xl = OrnsteinUhlenbeck(; κ = κ, σ = σ, length = 1000)
    μm = -0.01
    σm = 0.1
    m = AdditiveFunctional(Xl, μm .+ Xl.x .- 0.02, σm .* ones(length(Xl.x)); ρ = 1.0)
    ζ = tail_index(m)
    η, r = cgf_eigenvector(m, ζ)
    η, l = cgf_eigenvector(m, ζ, :left)
    @test η ≈ 0.0 atol = 1e-3
    ψ_tilde = stationary_distribution(DiffusionProcess(Xl.x, Xl.μx .+ ζ .* m.covariance, Xl.σx))
    @test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3
end


@testset "Cox-Ingersoll-Ross" begin
    gbar = 0.03
    σ_cir = 0.01
    Xc = CoxIngersollRoss(xbar = gbar, κ = κ, σ = σ_cir)
    m = AdditiveFunctional(Xc, Xc.x, zeros(length(Xc.x)))
    η_analytic = gbar * κ^2 / σ_cir^2 * (1 - sqrt(1 - 2 * σ_cir^2 / κ^2))
    @test cgf(m, 1.0) ≈ η_analytic rtol = 1e-2
end


@testset "generic AdditiveFunctional" begin
    Xs = OrnsteinUhlenbeck(; xbar = xbar, κ = κ, σ = σ, length = 200)
    m_diff = AdditiveFunctional(Xs, Xs.x, zeros(length(Xs.x)))
    @test tilted_generator(m_diff, 0) ≈ generator(Xs)

    # keyword and positional forms agree; scalars broadcast
    m_kw = AdditiveFunctional(Xs; drift = Xs.x, variance = 0.01^2)
    m_pos = AdditiveFunctional(Xs, Xs.x, 0.01 .* ones(length(Xs.x)))
    @test tilted_generator(m_kw, 1.3) ≈ tilted_generator(m_pos, 1.3)

    # the legacy AdditiveFunctionalDiffusion type agrees, including ρ ≠ 0
    m_ρ = AdditiveFunctional(Xs, Xs.x, 0.01 .* ones(length(Xs.x)); ρ = 0.5)
    m_afd = AdditiveFunctionalDiffusion(Xs, collect(Xs.x), 0.01 .* ones(length(Xs.x)); ρ = 0.5)
    @test m_ρ.covariance ≈ 0.5 .* 0.01 .* Xs.σx
    @test tilted_generator(m_ρ, 2.0) ≈ tilted_generator(m_afd, 2.0)
    @test cgf(m_ρ, 1.0) ≈ cgf(m_afd, 1.0) atol = 1e-10

    # the generic path (raw generator wrapped in a chain) agrees with the diffusion path
    Zx = ContinuousTimeMarkovChain(collect(Xs.x), Matrix(generator(Xs)))
    m_generic = AdditiveFunctional(Zx, collect(Xs.x), zeros(length(Xs.x)))
    @test cgf(m_generic, 1.0) ≈ cgf(m_diff, 1.0) atol = 1e-10
    @test tail_index(AdditiveFunctional(Zx, collect(Xs.x) .- 0.06, zeros(length(Xs.x)))) ≈
          tail_index(AdditiveFunctional(Xs, Xs.x .- 0.06, zeros(length(Xs.x)))) atol = 1e-3

    # two-state Markov multiplicative process: tail index satisfies Λ(ζ) = δ
    Zg = ContinuousTimeMarkovChain([0.0, 0.06], [-0.1 0.1; 0.5 -0.5])
    ν = 0.1
    δ_death = 0.05
    mz = AdditiveFunctional(Zg, [0.0, 0.06], [ν, ν])
    ζz = tail_index(mz; δ = δ_death)
    Λz(ξ) = maximum(real, eigvals(Matrix(generator(Zg)) + Diagonal(ξ .* [0.0, 0.06] .+ 0.5 .* ξ .^ 2 .* ν^2)))
    @test Λz(ζz) ≈ δ_death atol = 1e-4

    # state-shaped arrays for a multivariate process; state-independent growth has a closed form
    Y2 = ProductProcess(Xs, ContinuousTimeMarkovChain([0.5, 1.5], [-0.1 0.1; 0.2 -0.2]))
    mY = AdditiveFunctional(Y2, zeros(size(Y2)), fill(ν, size(Y2)))
    @test tail_index(mY; δ = δ_death) ≈ tail_index(ν^2 / 2, ν; δ = δ_death) atol = 1e-3

    # a covariance loading on one axis of a multivariate diffusion matches the univariate functional
    xs2 = collect(range(-0.5, 0.5, length = 30))
    ys2 = collect(range(-0.5, 0.5, length = 10))
    Xmv = MultivariateDiffusionProcess((; x = xs2, y = ys2);
        drift = (; x = repeat(-κ .* xs2, 1, length(ys2)), y = repeat(reshape(-0.2 .* ys2, 1, :), length(xs2), 1)),
        variance = (; x = σ^2, y = 0.03^2))
    X1 = DiffusionProcess(xs2, -κ .* xs2, σ .* ones(length(xs2)))
    σm2, c2 = 0.1, 0.5 * 0.1 * σ
    m_mv = AdditiveFunctional(Xmv; drift = repeat(xs2, 1, length(ys2)), variance = σm2^2, covariance = (; x = c2))
    m_1d = AdditiveFunctional(X1; drift = xs2, variance = σm2^2, covariance = c2)
    @test cgf(m_mv, 1.0) ≈ cgf(m_1d, 1.0) rtol = 1e-6

    # validation
    @test_throws ArgumentError AdditiveFunctional(Zg; drift = 0.0, variance = ν^2, covariance = [0.1, 0.1])
    @test_throws ArgumentError AdditiveFunctional(Xs; drift = 0.0, variance = 0.01^2, covariance = 1.0)
    @test_throws ArgumentError AdditiveFunctional(Xmv; drift = 0.0, variance = σm2^2, covariance = (; x = 1.0))
    @test_throws ArgumentError AdditiveFunctional(Xmv; drift = 0.0, variance = σm2^2, covariance = (; z = 0.0))
    @test_throws ArgumentError AdditiveFunctional(Xs; drift = 0.0, variance = -1.0)
    @test_throws DimensionMismatch AdditiveFunctional(Zg, [0.0], [ν, ν])
    @test_throws DimensionMismatch AdditiveFunctional(Y2, zeros(reverse(size(Y2))), fill(ν, size(Y2)))
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
    @test d2y[1] ≈ 1.0 atol = 1e-12
    @test d2y[end] ≈ -1997.0 atol = 1e-8

    xs = range(-1.0, stop = 1.0, length = 31)
    ys = range(-2.0, stop = 2.0, length = 41)
    grid = (xs, ys)
    named_grid = (; x = xs, y = ys)
    f = [x^2 + y^3 + x * y for x in xs, y in ys]

    fx = FirstDerivative(grid, f, 1; direction = :forward)
    fy = FirstDerivative(grid, f, 2; direction = :backward)
    @test size(fx) == size(f)
    @test fx[15, 20] ≈ 2 * xs[15] + ys[20] atol = 1e-1
    @test fy[15, 20] ≈ 3 * ys[20]^2 + xs[15] atol = 2e-1

    fxx = SecondDerivative(grid, f, 1, 1)
    fyy = SecondDerivative(grid, f, 2)
    fxy_up = SecondDerivative(grid, f, 1, 2; direction = :up)
    fxy_down = SecondDerivative(grid, f, 1, 2; direction = :down)
    @test fxx[15, 20] ≈ 2.0 atol = 1e-10
    @test fyy[15, 20] ≈ 6 * ys[20] atol = 1e-10
    @test fxy_up[15, 20] ≈ 1.0 atol = 1e-10
    @test fxy_down[15, 20] ≈ 1.0 atol = 1e-10
    @test FirstDerivative(named_grid, f, :x; direction = :forward) ≈ fx
    @test SecondDerivative(named_grid, f, :x, :y; direction = :up) ≈ fxy_up
    @test_throws ArgumentError FirstDerivative(grid, f, 3)
    @test_throws ArgumentError FirstDerivative(named_grid, f, :z)
    @test_throws DimensionMismatch SecondDerivative(grid, f[1:end-1, :], 1, 1)
    @test_throws ArgumentError SecondDerivative(grid, f, 1, 2; direction = :sideways)

    # multi-dimensional operators are lazy, like the one-dimensional ones
    @test !(fx isa Array)
    @test fx isa AbstractMatrix
    @test collect(fx) == fx

    # :up/:down and :forward/:backward are synonyms everywhere
    @test FirstDerivative(x, y; direction = :up) == dy
    @test FirstDerivative(x, y; direction = :down) == dy_back
    @test FirstDerivative(grid, f, 1; direction = :up) == fx
    @test SecondDerivative(grid, f, 1, 2; direction = :forward) == fxy_up
    @test SecondDerivative(grid, f, 1, 2; direction = :backward) == fxy_down
    @test_throws ArgumentError FirstDerivative(x, y; direction = :sideways)

    # cross derivatives do not take boundary conditions
    @test_throws ArgumentError SecondDerivative(grid, f, 1, 2; bc = (1, 0))
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
    𝔸1 = generator(X1)
    𝔸2 = generator(X2)
    Q = [-0.1 0.1; 0.2 -0.2]
    J = jointoperator([𝔸1, 𝔸2], Q)
    @test size(J) == (100, 100)
    # rows should sum to zero (generator property)
    @test maximum(abs.(sum(Matrix(J), dims = 2))) < 1e-10
    # principal eigenvalue of a generator should be zero
    Jdense = Matrix(J)
    η, r = InfinitesimalGenerators.principal_eigenvalue(Jdense)
    @test η ≈ 0.0 atol = 1e-8
    @test all(r .> 0)
    A = [1.0 0.2; 0.3 0.5]
    η_stalled, r_stalled = @test_logs (:warn, r"Inverse iteration") InfinitesimalGenerators.principal_eigenvalue(A; η0 = 1.23, maxiter = 0)
    @test η_stalled == 1.23
    @test norm(r_stalled) ≈ 1.0
    @test_throws DimensionMismatch jointoperator([𝔸1], Q)
    @test_throws DimensionMismatch jointoperator([𝔸1, generator(OrnsteinUhlenbeck(; κ = 0.1, σ = 0.02, length = 60))], Q)
    @test_throws DimensionMismatch jointoperator([𝔸1, 𝔸2], [-0.1 0.1 0.0; 0.2 -0.2 0.0])
end
