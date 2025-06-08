using Pkg; Pkg.activate(".")
import Turing # should be loaded before MacroModelling
using Plots
using StatsPlots # should be loaded before MacroModelling
using MacroModelling
# using Zygote
# import Turing, Pigeons
# import Turing: NUTS, sample, logpdf, AutoZygote
import Optim, LineSearches
using Random, DataFrames, CSV, AxisKeys
# CSV, DataFrames, MCMCChains, AxisKeys
import DynamicPPL
using JLD2, FredData
using TimeSeriesEcon, Dates
using Base.Threads
using LinearAlgebra
using Plots
using FiniteDiff
using ProgressBars
using TimerOutputs
using StatsBase
using GLM
import ModelBaseEcon
import StateSpaceEcon
import LsqFit
using Optim
using Distributions
import MAT

BLAS.set_num_threads(16)
const to = TimerOutput()

DE = TimeSeriesEcon.DataEcon

# models
@eval TimeSeriesEcon begin
    function _fconvert_lower(F_to::Type{<:Union{Monthly,Quarterly{N},Quarterly,HalfYearly,HalfYearly{N},Yearly,Yearly{N},Weekly,Weekly{N}}}, t::TSeries{BDaily}, aggregator::Function; ref=:end, skip_all_nans::Bool=false, skip_holidays::Bool=false, holidays_map::Union{Nothing, TSeries{BDaily}}=nothing) where {N}
        rng_from = rangeof(t)
        
        dates = [Dates.Date(val) for val in rng_from]
        if holidays_map !== nothing
            dates = dates[holidays_map[rng_from].values]
        elseif skip_holidays
            holidays_map = copy(TimeSeriesEcon.getoption(:bdaily_holidays_map))
            dates = dates[holidays_map[rng_from].values]
        end
    
        trim = aggregator ∈ (first, last) ? ref : :both
        if skip_all_nans == true
            nanmap = .!isnan.(t)
            if holidays_map == nothing
                holidays_map = copy(TimeSeriesEcon.getoption(:bdaily_holidays_map))
            end
            if holidays_map === nothing
                holidays_map = TSeries(first(rng_from)-600, trues(length(rng_from)+1200))
            end
            holidays_map[rng_from] .= holidays_map[rng_from] .& nanmap
        end
        fi, li, trunc_start, trunc_end = _fconvert_using_dates_parts(F_to, rng_from, trim=trim, holidays_map=holidays_map)
        out_index = _get_out_indices(F_to, dates)
        
        ret = aggregator == mean ? Vector{Float64}() : Vector{eltype(t)}()
        for target in unique(out_index)
            target_range = collect(rng_from)[out_index .== target]
            vals = cleanedvalues(t[target_range[begin]:target_range[end]]; skip_all_nans=skip_all_nans, skip_holidays=skip_holidays, holidays_map=holidays_map)
            push!(ret, aggregator(vals))
        end

        return copyto!(TSeries(eltype(ret), fi+trunc_start:li-trunc_end), ret[begin+trunc_start:end-trunc_end])
    end
end


# functions
include("$(pwd())/src/utils.jl")
include("$(pwd())/src/data.jl")
include("$(pwd())/src/data_dicts.jl")
include("$(pwd())/src/draws.jl")
include("$(pwd())/src/kalman.jl")
include("$(pwd())/src/model_configs.jl")
include("$(pwd())/src/metropolis_within_gibbs.jl")
include("$(pwd())/src/storage.jl")
include("$(pwd())/src/analysis_functions.jl")
