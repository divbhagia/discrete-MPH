This repository contains the code for the paper "Discrete-Time Mixed Proportional Hazard Model" by Div Bhagia. 

The analysis for the paper was conducted using Julia 1.10.4. To replicate the results, ensure you have [Julia](https://julialang.org/downloads/) installed on your machine, then execute the following commands in the terminal:

```shell
git clone https://github.com/divbhagia/DiscreteMPH.git
cd DiscreteMPH
make
```

The `make` command will install the necessary packages, execute the code, and compile the document `draft/figures_and_tables.tex`, which contains the final output for the paper. The code also creates a new folder, `output`, which will contain the plots and tables included in the final document, and all simulation results will be saved in `output/sims`.