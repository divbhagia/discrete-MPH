
#######################################################
# Function to display rounded numbers
#######################################################

function pretty(v, digs=3)
    new_v = Vector{String}(undef, length(v))
    for i in eachindex(v)
        fmtd_num = string(round(v[i], digits=digs))
        decimal_pos = findfirst('.', fmtd_num)
        while length(fmtd_num) - decimal_pos < digs 
            fmtd_num *= "0"
        end
        new_v[i] = fmtd_num
    end
    return new_v
end

#######################################################
# Function to run simulation
#######################################################

function run_sim(iters, n, T_bar, p, ϕ, ν_pars, opt, folder)

    # DGP quantities and true parameters
    λ, γ, μ, μ_nrm = dgp_quants(T_bar, ν_pars, ϕ, opt)
    Θ = vcat(γ, λ, μ_nrm[2:end])
    nrm = ϕ[1] * μ[1]

    # Run siμlations and save thta_hats to file
    Θ_hats = zeros(length(Θ), iters)
    for i in 1:iters
        if i % 20 == 0
            println(100 * i / iters, "% simulations done")
        end
        g = sim_data(n, T_bar, p, ϕ, ν_pars, opt)
        Θ_hats[:, i] = gmm(g, nrm)
    end
    @save "$folder/$(opt)-n$n-ests.jld" Θ_hats 

    
    # Bias, St.Dev, RMSE in a dataframe for all parameters
    bias = mean(Θ_hats, dims=2) .- Θ
    stdev = std(Θ_hats, dims=2)
    rmse = sqrt.(bias.^2 .+ stdev.^2)
    df = DataFrame(bias=pretty(bias), stdev=pretty(stdev), rmse=pretty(rmse))
    @save "$folder/$(opt)-n$n-res.jld" df

    # Structural hazard dataframe
    λ̂ = unstacker(T_bar, mean(Θ_hats, dims=2))[2]
    σ = unstacker(T_bar, stdev)[2]
    df = DataFrame(λ=λ, λ̂=λ̂, σ=σ)
    @save "$folder/$(opt)-n$n-lmbda.jld" df

end

#######################################################
# Function to concate dataframes and save to file
#######################################################

function concat_save(df1, df2, df3, names, filepath)

    # Create table
    empty = fill("", size(df1, 1))
    df = hcat(hcat(empty, df1), hcat(empty, df2), makeunique=true)
    df = hcat(df, hcat(empty, df3), makeunique=true) 
    df = hcat(names, df, makeunique=true)

    # Convert to LaTeX
    latex_code = latexify(df; env=:table, latex=false)
    lines = split(latex_code, '\n')
    latex_code_updated = join(lines[3:end-2], '\n')

    # Write to file & return dataframe
    open(filepath, "w") do file
        write(file, latex_code_updated)
    end
    return df

end

#######################################################
# Save simulation parameters to a file
#######################################################

function save_simdetails(folder, iters, T_bar, p, ϕ, ν_pars, n_vec, opts)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH.MM.SS")
    open("$folder/simdetails.txt", "w") do io
        println(io, "Timestamp: $timestamp")
        println(io, "Number of iterations: $iters")
        println(io, "T_bar: $T_bar")
        println(io, "π: $p")
        println(io, "ϕ: $ϕ")
        println(io, "ν parameters: $ν_pars")
        println(io, "Sample sizes: $n_vec")
        println(io, "Options: $opts")
    end
end

#######################################################