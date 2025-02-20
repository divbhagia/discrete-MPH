using Random, Distributions, Optim
using DataFrames, Latexify
using JLD2, Dates, Plots, Measures, Printf
include("utils_sim.jl")
include("utils_est.jl")
include("utils_oth.jl")

# Install packages using packages.jl if missing packages

#######################################################
# Parameters
#######################################################

# Set seed
Random.seed!(1117)

# Parameters
iters = 5000
T_bar = 6
p = 0.5
ϕ = [1, 2]
ν_pars = [2, 2]
n_vec = [1500, 5000, 10000]
n1, n2, n3 = n_vec
opts = ["inc", "dec", "con"]
run_sim_again_flag = true

# Folder paths
outfolder = "output"
simfolder = "$outfolder/sims"

#######################################################
# Run simulations and save output
#######################################################

if run_sim_again_flag

    # Overwrite existing folders
    isdir(outfolder) && rm(outfolder; recursive=true)
    mkpath(simfolder)

    # Save simulation details
    save_simdetails(outfolder, iters, T_bar, p, ϕ, ν_pars, n_vec, opts)

    # Run simulations
    for opt in opts
        for n in n_vec
            println("Running simulations for $opt with n=$n")
            run_sim(iters, n, T_bar, p, ϕ, ν_pars, opt, simfolder)
        end
    end

end

#######################################################
# Create table with results
#######################################################

# Pull saved output to a dictionary
restabs = Dict()
times = Dict()
for opt in opts
    for n in n_vec
        key = "$(opt)-n$n"
        resfile = jldopen("$simfolder/$key-res.jld", "r")
        restabs["$key"] = read.(Ref(resfile), keys(resfile))[1]
        timefile = jldopen("$simfolder/$key-times.jld", "r")
        times["$key"] = read.(Ref(timefile), keys(timefile))[1]
    end
end

# Parameter names for table
names = ["\\(\\gamma\\)"]
for i in 1:T_bar
    push!(names, "\\(\\lambda_$i\\)")
end
for i in 2:T_bar
    push!(names, "\\(\\tilde{\\mu}_$i\\)")
end

# Create a table for each option with three sample sizes
for opt in opts
    tab = concat_save(restabs["$(opt)-n$n1"], 
                    restabs["$(opt)-n$n2"], 
                    restabs["$(opt)-n$n3"], 
                    names, "$outfolder/tab_$(opt).tex")
end

# Create tables for computation times
for opt in opts
    file = open("$outfolder/tab_$(opt)_times.tex", "w")
    t1 = round(times["$(opt)-n$n1"], digits=2)
    t2 = round(times["$(opt)-n$n2"], digits=2)
    t3 = round(times["$(opt)-n$n3"], digits=2)
    write(file, "\\multicolumn{1}{l}{Compute time} & & \\multicolumn{3}{c}{$t1} & & \\multicolumn{3}{c}{$t2} & & \\multicolumn{3}{c}{$t3} \\\\ \n")
    close(file)
end

#######################################################
# Plots
#######################################################

# Pull saved output to dictionaries
λ, λ̂, σ,  = Dict(), Dict(), Dict()
for opt in opts
    for n in n_vec
        key = "$(opt)-n$n"
        lfile = jldopen("$simfolder/$key-lmbda.jld", "r")
        λ["$opt"] = read.(Ref(lfile), keys(lfile))[1][!, "λ"]
        λ̂["$key"] = read.(Ref(lfile), keys(lfile))[1][!, "λ̂"]
        σ["$key"] = read.(Ref(lfile), keys(lfile))[1][!, "σ"]
    end
end

# Create and store all plots in a dictionary
plots = Dict()
for opt in opts
    for n in n_vec
        key = "$(opt)-n$n"
        λ_i, λ̂_i, σ_i = λ[opt], λ̂[key], σ[key]
        plots[key] = plot(λ̂_i, ribbon=1.96*σ_i, color="#A6A29D")
        plot!(plots[key], λ̂_i, color=:black, linewidth=1.25)
        plot!(plots[key], λ_i, color=:red, linewidth=1.25, linestyle=:dot)
        plot!(plots[key], title="n=$n", titlefontsize=5, legend=false)
        plot!(plots[key], xtickfontsize=4, ytickfontsize=4)
    end
end

# Plot for increasing hazard
plot_inc = plot(plots["inc-n$n1"], plots["inc-n$n2"], plots["inc-n$n3"],
                layout=(1, 3), size=(250, 70), margin=-1.5mm,
                ylim=(0.125, 0.8), yticks=0:0.1:1)
savefig(plot_inc, "$outfolder/plot_inc.pdf")

# Plot for decreasing hazard
plot_dec = plot(plots["dec-n$n1"], plots["dec-n$n2"], plots["dec-n$n3"],
                layout=(1, 3), size=(250, 70), margin=-1.5mm, 
                ylim=(0.125, 0.55), yticks=0:0.1:1)
savefig(plot_dec, "$outfolder/plot_dec.pdf")

# Plot for constant hazard
plot_con = plot(plots["con-n$n1"], plots["con-n$n2"], plots["con-n$n3"],
                layout=(1, 3), size=(250, 70), margin=-1.5mm, 
                ylim=(0.15, 0.55), yticks=0:0.1:1)
savefig(plot_con, "$outfolder/plot_con.pdf")

# Save legend separately
leg = plot((1:2)', framestyle=:none, 
        legend_column=2, color = [:black :red],
        linestyle = [:solid :dot], linewidth=1.5,
        labels=[" Estimate   " " True"], legendfontsize=5,
        fg_legend = :transparent, 
        margins = -2mm,
        legend = :top,  size=(250, 20))
savefig(leg, "$outfolder/plot_legend.pdf")

#######################################################