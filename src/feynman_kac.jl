
#========================================================================================

Compute u(x_t, t) = E[∫t^T e^{-∫ts V(x_τ, τ)dτ}f(x_s, s)ds + e^{-∫tT V(x_τ, τ)dτ}ψ(x_T)|x_t = x]
using
0 = (u_{t+1} - u_{t})/dt + 𝔸u_t + f
that is
(I - 𝔸dt)u_t =  u_{t+1} + f dt
========================================================================================#

function feynman_kac_backward(x, μx, σx; ψ::AbstractVector, t::AbstractVector = range(0, 100, step = 1/12), f::T = zeros(length(x)), V::T = zeros(length(x))) where {T <: Union{AbstractVector, AbstractMatrix}}
	u = zeros(length(x), length(t))
	u[:, length(t)] = ψ
	Δ = make_Δ(x)
	𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
	if (T <: AbstractVector)
		dt = t[2] - t[1]
		𝔹 = factorize(I - build_operator!(𝔸, Δ, V .* dt, μx .* dt, 0.5 .* σx.^2 .* dt))
		for i in (length(t)-1):(-1):1
			ψ = ldiv!(𝔹, u[:, i+1] .+ f .* dt)
			u[:, i] = ψ
		end
	elseif T <: AbstractVector
		for i in (length(t)-1):(-1):1
			dt = t[i+1] - t[i]
			𝔹 = I - build_operator!(𝔸, Δ, V .* dt, μx .* dt, 0.5 .* σx.^2 .* dt)
			ψ = 𝔹 \  (u[:, i+1] .+ f .* dt)
			u[:, i] = ψ
		end
	else
		for i in (length(t)-1):(-1):1
			dt = t[i+1] - t[i]
			𝔹 = (I - build_operator!(𝔸, Δ, V[:, i] .* dt, μx .* dt, 0.5 .* σx.^2 .* dt))
			ψ = 𝔹 \ (u[:, i+1] .+ f[:, i] .* dt)
			u[:, i] = ψ
		end
	end
	return u
end

#========================================================================================

Compute u(x_t, T)= E[∫t^T e^{-∫ts V(x_τ)dτ}f(x_s)ds + e^{-∫tTV(x_τ)dτ} ψ(x_T)|x_t = x]

========================================================================================#

function feynman_kac_forward(x, μx, σx; ψ::AbstractVector, t::AbstractVector = range(0, 100, step = 1/12), f::AbstractVector = zeros(length(x)), V::AbstractVector = zeros(length(x)))
	u = feynman_kac_backward(x, μx, σx; ψ = ψ, t = .- reverse(t), f = f, V = V)
	return u[:,end:-1:1]
end
