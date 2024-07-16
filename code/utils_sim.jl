
##########################################################
# Data generating process
##########################################################

# Weibull Hazard
function weibull(x; opt="con")
    if opt == "con"
        a, b = 3, 1
    elseif opt == "inc"
        a, b = 3.15, 1.25
    elseif opt == "dec"
        a, b = 2, 0.75
    end
    return (b/a) * (x / a) ^ (b - 1) 
end

# Structural hazard
function λ_func(T_bar, opt)
    λ = zeros(T_bar)
    for t in 1:T_bar
        λ[t] = weibull(t, opt=opt)
    end
    return λ
end

# Function that outputs K moments of Beta Distribution
function beta_moms(K, pars)
    a, b = pars
    μ = zeros(K)
    μ[1] = a / (a + b)
    for k in 2:K
        μ[k] = μ[k - 1] * (a + k - 1) / (a + b + k - 1)
    end
    μ_nrm = [μ[k] / (μ[1] ^ (k)) for k in 1:K]
    return μ, μ_nrm
end

# DGP Quants
function dgp_quants(T_bar, ν_pars, ϕ, opt)
    λ = λ_func(T_bar, opt)
    γ = ϕ[2] / ϕ[1]
    μ, μ_nrm = beta_moms(T_bar, ν_pars)
    return λ, γ, μ, μ_nrm
end

##########################################################
# Function to siμlate data & get data moments
##########################################################

function sim_data(n, T_bar, p, ϕ, ν_pars, opt="con")

    """
    n: number of individuals
    T_bar: number of periods
    p: probability of Z = 1
    ϕ: [ϕ_0, ϕ_1]
    ν_pars: [a, b] for Beta Distribution
    opt: "con", "inc", "dec" for structural hazard
    """

    # Generate structural hazard
    J = length(ϕ)
    λ = λ_func(T_bar, opt)

    # Generate Z & nu
    Z = rand(Binomial(1, p), n)
    ν = rand(Beta(ν_pars[1], ν_pars[2]), n)
        
    # Individual exit probabilities
    exit_probs = zeros(n, T_bar)
    for i in 1:n
        exit_probs[i, :] = λ .* ϕ[Z[i]+1] .* ν[i]  
    end

    # If any exit_prob > 0 or <1, print error
    if any(exit_probs .< 0) || any(exit_probs .> 1)
        println("Error: Exit probability out of bounds")
    end

    # Exits
    exited = zeros(n)
    crit = rand(n, T_bar)
    exit_ = exit_probs .> crit
    for t in 1:T_bar
        exited .+= exit_[:, t]
        exit_[exited .!= 0, (t+1):end] .= false
    end

    # Calculate density
    g = zeros(T_bar, J)
    for t in 1:T_bar
        for z in 1:J
            g[t, z] = mean(exit_[Z .== z-1, t])
        end
    end

    return g
end

##########################################################