using Random, Distributions, Optim
using JLD2
include("utils_sim.jl")
include("utils_est.jl")

#### Code to test functions

# Set seed
Random.seed!(1117)

# Parameters
n = 1000000
p = 0.5
T_bar = 6
ϕ = [1, 2]
ν_pars = [2, 2]
opt = "inc"

# DGP quantities and true parameters
λ, γ, μ, μ_nrm = dgp_quants(T_bar, ν_pars, ϕ, opt)
Θ = vcat(γ, λ, μ_nrm[2:end])
nrm = ϕ[1] * μ[1]

# Model moments
g_model = model_moms(γ, λ, μ_nrm, nrm)

# Simulate data
g = sim_data(n, T_bar, p, ϕ, ν_pars, opt)

# Check error in model vs data moments
println("Model vs data moments: ", sum((g - g_model).^2))

# Estimate parameters when true moments given
@time Θ_hat = gmm(g, nrm)

# Print error in estimation
error = sum((Θ_hat - Θ).^2)
println("Estimates vs true parameters: ", sum((Θ_hat - Θ).^2))


