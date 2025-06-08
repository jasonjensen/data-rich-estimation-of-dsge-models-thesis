# Data-Rich Estimation of DSGE models

This repo contains code for reproducing the content of my MA thesis in economics *Explorations in Data-Rich Estimation of DSGE Models*.

## Organization

* **analysis**: contains scripts for reproducing the tables and figures in the paper.
* **control**: contains scripts for running the different estimation variants.
* **data**: contains both input and output data. Empty in the github repo. Folder contents can be accessed [here](http://data-rich-ma-thesis-jensen.s3-website.ca-central-1.amazonaws.com/index.html
    ).
* **graphs**: contains plots produced by the scripts in the analysis folder. Empty in the github repo.
* **models**: contains the MacroModelling.jl linearized SW and SWFF models.
* **src**: contains functions called by the other parts of the code.

## Example usage:
1. Clone the repo.
2. Download [data/inputs/fred_data.jl](https://data-rich-ma-thesis-jensen.s3.ca-central-1.amazonaws.com/data%2Finput%2Ffred_data.jld2) and place it at that path in the repo folder.
3. Start julia in the folder and run:
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
4. run `include("setup.jl")`. 
5. Run one of the estimations, i.e. `include("control/run_sw_reg_estimation_2007Q2.jl")`


