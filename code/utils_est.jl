
##########################################################
# Function that gives model moments
##########################################################

function model_moms(γ, λ, μ_nrm, nrm)

    """
    γ: ϕ_1 / ϕ_0
    λ: structural hazard
    μ_nrm: normalized moments (μ_k / μ_1^k)
    nrm: ϕ_0 * μ_1 (normalization constant)
    """

    T_bar = length(λ)
    c = zeros(T_bar, T_bar, 2)
    c[:, 1, :] .= 1 
    for t in 2:T_bar
        for k in 2:T_bar
            c[t, k, :] .= c[t-1, k, :] - λ[t-1] * c[t-1, k-1, :]
        end
    end
    for k in 1:T_bar
        c[:, k, 1] .= (nrm ^ (k)) * c[:, k, 1]
        c[:, k, 2] .= ((γ * nrm)^(k)) * c[:, k, 2]
    end
    g = λ .* hcat(c[:, :, 1] * μ_nrm, c[:, :, 2] * μ_nrm)
    return g

end

##########################################################
# Other helper functions for GMM estimation
##########################################################

# Function to unstack parameters
function unstacker(T_bar, Θ)
    γ = Θ[1]
    λ = Θ[2:T_bar+1]
    μ_nrm = vcat(1, Θ[T_bar+2:end])
    return γ, λ, μ_nrm
end

# Define objective function
function objfun(Θ, g, nrm)
    γ, λ, μ_nrm = unstacker(T_bar, Θ)
    g_model = model_moms(γ, λ, μ_nrm, nrm)
    return sum((g .- g_model).^2)
end

# Numeric gradient
function num_grad(x, g, nrm)
    epsilon = 1e-8
    grad = zeros(length(x))
    for i in 1:length(x)
        epsvec = zeros(length(x))
        epsvec[i] = epsilon
        objplus = objfun(x + epsvec, g, nrm)
        objminus = objfun(x - epsvec, g, nrm)
        grad[i] = (objplus - objminus) / (2 * epsilon)
    end
    return grad
end

##########################################################
# GMM estimation
##########################################################

function gmm(g, nrm, bounds = true, usegrad = false)

    # Initial guess
    T_bar, J = size(g)
    num_vars = 2T_bar + J - 2
    x0 = 0.5 * ones(num_vars)
    x0[T_bar+2:end] .= 1.5

    # Optimization options
    options = Optim.Options(g_tol = 1e-32, iterations = 100000,
                                f_tol = 1e-32, x_tol = 1e-32,
                                show_trace = false)  
    optimizer = Optim.LBFGS()

    # Objective function and gradient
    obj_fun = (x) -> objfun(x, g, nrm)
    grad_fun! = (grad, x) -> (grad[:] = num_grad(x, g, nrm))

    # Run optimization without bounds
    if bounds == false
        if usegrad
            result = optimize(obj_fun, grad_fun!, x0, optimizer, options)
        else
            result = optimize(obj_fun, x0, optimizer, options)
        end
    end
    
    # Run optimization with bounds
    if bounds == true
        lb = zeros(num_vars)
        lb[T_bar+2:end] .= 1
        ub = fill(Inf, num_vars)
        if usegrad
            result = optimize(obj_fun, grad_fun!, lb, ub, x0, 
                Fminbox(optimizer), options)
        else
            result = optimize(obj_fun, lb, ub, x0, 
                Fminbox(optimizer), options)
        end
    end

    # Return estimated parameters
    x_hat = result.minimizer
    return x_hat

end

##########################################################