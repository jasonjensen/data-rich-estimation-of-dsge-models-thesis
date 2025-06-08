"""
    transform_fred_data(key, data_dict::Dict, return_rng=nothing, pop_data=nothing; assume_no_split=false, treatment_rng=return_rng)

Fetch, combine (if necessary), and transform economic time series data, typically from FRED.

This function retrieves raw time series data based on a `key` from `data_dict`.
It handles potential series combinations (e.g., ratios, differences) specified by
operators like '/', '-', '+' in `fred_codes[key]`. The raw series is converted to
quarterly frequency. Finally, a transformation (e.g., log, detrend, difference)
specified by `fred_transformations[key]` is applied.

# Arguments
- `key::String`: The primary key identifying the series or combination of series to transform. This key is used to look up codes in `fred_codes` and transformation types in `fred_transformations`.
- `data_dict::Dict`: A dictionary where keys are FRED series codes (or components of a combined series) and values are objects containing `.data` (a DataFrame-like structure with `date` and `value` columns).
- `return_rng=nothing`: An optional `MIT` range. If provided, the final transformed series and `pop_data` will be truncated to this range.
- `pop_data=nothing`: An optional `TSeries` object representing population data, used for per-capita transformations. Must align with `return_rng` if provided.

# Keyword Arguments
- `assume_no_split=false`: If `true`, assumes `fred_codes[key]` is a single series code and does not attempt to split it by operators.
- `treatment_rng=return_rng`: An `MIT` range used for calculating means or trends during transformations (e.g., demeaning over this specific period). Defaults to `return_rng`.

# Returns
- `ts::TSeries`: The transformed quarterly time series.

# Details
- **Frequency Conversion**: All input series are converted to `Quarterly` frequency. Daily series (`isbdaily`) are aggregated using `fconvert(Quarterly, ..., skip_all_nans=true)`.
- **Transformations (`fred_transformations` mapping):**
    - `0`: Demean.
    - `1`: Log and demean.
    - `2`: Log per capita (for specific keys like "RGDP"), then linear detrend and demean, scaled by 100. For other keys, log, linear detrend, demean, and scale by 100.
    - `3`: Log difference, scaled by 100, and demeaned.
    - `4`: Log, linear detrend, demean, and scale by 100.
    - `5`: Per capita, linear detrend, and demean.
    - `6`: No transformation (returns the quarterly series as is, within `return_rng`).
    - `7`: Demean (originally intended for FedFunds, now just demeans).
    - `8`: Log per capita difference (for specific keys like "RGDP"), scaled by 100, demeaned. For other keys, log difference, scaled by 100, demeaned.
- **Dependencies**: Relies on global `fred_codes` and `fred_transformations` dictionaries. Requires `TimeSeriesEcon` for `TSeries` operations and frequency conversions.
"""
function transform_fred_data(key, data_dict::Dict, return_rng=nothing, pop_data=nothing; assume_no_split=false, treatment_rng=return_rng)
    # TimeSeriesEcon.clear_holidays_map()
    # TimeSeriesEcon.set_holidays_map("US", "NY")
    # if key == "WTI"
        TimeSeriesEcon.set_holidays_map("US", "NY")
    # end
    
    if !assume_no_split
        fred_series_key = fred_codes[key]
        op =    contains(fred_series_key, '/') ? '/' :
                contains(fred_series_key, '-') ? '-' :
                contains(fred_series_key, '+') ? '+' :
                '.'
    else
        op = '.'
    end
    
    if op == '.'
        # convert to TSeries
        freq = infer_frequency(data_dict[key].data)
        ts_high_freq = TSeries(dat2mit(freq, first(data_dict[key].data.date)), data_dict[key].data.value)
    else
        key1, key2 = split(fred_series_key, op)
        freq1 = infer_frequency(data_dict[key1].data)
        freq2 = infer_frequency(data_dict[key2].data)
        
        ts_high_freq1 = TSeries(dat2mit(freq1, first(data_dict[key1].data.date)), data_dict[key1].data.value)
        ts_high_freq2 = TSeries(dat2mit(freq2, first(data_dict[key2].data.date)), data_dict[key2].data.value)
        if freq1 !== freq2
            if freq1 > freq2
                ts_high_freq2 = fconvert(freq1, ts_high_freq2, ref=:end)
            else
                ts_high_freq1 = fconvert(freq2, ts_high_freq1, ref=:end)
            end
        end
        
        ts_high_freq =  op == '/' ? ts_high_freq1 ./ ts_high_freq2 :
                        op == '-' ? ts_high_freq1 .- ts_high_freq2 :
                        op == '+' ? ts_high_freq1 .+ ts_high_freq2 :
                        ts_high_freq1
    end

    # convert to Quarterly
    if isquarterly(ts_high_freq)
        ts_full = ts_high_freq
    elseif isbdaily(ts_high_freq)
        if key == "WTI"
            ts_high_freq[ts_high_freq.firstdate - 1] = first(values(ts_high_freq))
        end
        ts_full = fconvert(Quarterly, ts_high_freq; skip_all_nans=true)
    else
        ts_full = fconvert(Quarterly, ts_high_freq)
    end

    transformation = fred_transformations[key]
    ts = copy(ts_full)
    if return_rng !== nothing
        ts = ts[return_rng]
        pop_data = pop_data[return_rng]
    end
    # @show key, transformation

    if transformation == 0
        # demean
        # ts = ts[rng]
        ts = ts .- mean(ts[treatment_rng])
    elseif transformation == 1
        # log and demean
        ts = log.(ts)
        ts = ts .- mean(ts[treatment_rng])
    elseif transformation == 2
        if key ∈ ["RGDP", "RCONS", "RINV"]
            ts = ts * 1_000_000 # billions of dollars / thousands for pop
            # Note: All per capita variables are calculated using the adult population series (CNP16OV).
            # Linear detrended log per capital
            ts = log.(ts ./ pop_data)
        else 
            ts = log.(ts)
        end
        # TODO: may be assuming ranges have the same start date
        time = collect(0.0:length(ts[treatment_rng])-1)
        time_full = collect(0.0:length(ts)-1)
        X = hcat(ones(length(ts[treatment_rng])), time)
        betas = inv(X'*X)*X'*values(ts[treatment_rng])
        ts.values = ts.values .- time_full*betas[2]
        ts = ts .- mean(ts[treatment_rng]) # might not need to demean?
        ts = ts .* 100 # non-core series will be normalized, so this doesn't matter
    elseif transformation == 3
        # log differenced and demeaned
        ts = diff(log.(ts_full))

        if return_rng !== nothing
            ts = ts[return_rng]
        end
        ts = ts .* 100
        ts = ts .- mean(ts[treatment_rng])
    elseif transformation == 4
        # detrended log
        ts = log.(ts)
        time = collect(0.0:length(ts[treatment_rng])-1)
        time_full = collect(0.0:length(ts)-1)
        X = hcat(ones(length(ts[treatment_rng])), time)
        betas = inv(X'*X)*X'*values(ts[treatment_rng])
        ts.values = ts.values .- time_full*betas[2]
        ts = ts .- mean(ts[treatment_rng]) # might not need to demean?
        ts = ts .* 100
    elseif transformation == 5
        # detrended per capita level
        ts = ts ./ pop_data
        time = collect(0.0:length(ts[treatment_rng])-1)
        time_full = collect(0.0:length(ts)-1)
        X = hcat(ones(length(ts[treatment_rng])), time)
        betas = inv(X'*X)*X'*values(ts[treatment_rng])
        ts.values = ts.values .- time_full*betas[2]
        ts = ts .- mean(ts[treatment_rng]) # might not need to demean?
    elseif transformation == 6
        return ts
    elseif transformation == 7
        # detrend and demean
        # for FedFunds rate
        # update!: just demean and divide by 4
        ts = ts .- mean(ts[treatment_rng])
        # time = collect(0.0:length(ts)-1)
        # X = hcat(ones(length(ts)), time)
        # betas = inv(X'*X)*X'*values(ts)
        # ts.values = ts.values .- time*betas[2]
        # ts = ts .- mean(ts) # might not need to demean?
        ts = ts #./ 4
    elseif transformation == 8
        
        if key ∈ ["RGDP", "RCONS", "RINV"]
            ts = ts_full * 1000000 # billions of dollars / thousands for pop
            # Note: All per capita variables are calculated using the adult population series (CNP16OV).
            # Linear detrended log per capital
            ts = log.(ts_full ./ pop_data)
        else 
            ts = log.(ts_full)
        end
        ts = diff(ts)
        ts = ts .- mean(ts[treatment_rng]) # might not need to demean?
        ts = ts .* 100
        ts = ts[return_rng]
    end

    return ts
end

"""
    normalize_noncore_series!(data::Workspace, groups::Workspace)

Normalize (demean and standardize) non-core time series within a data Workspace.

This function iterates through the series in `data.X`. If a series is not found
in the list of `core` series defined in `groups.core`, it is transformed by
subtracting its mean and dividing by its standard deviation. The modification
is done in-place on `data.X`.

# Arguments
- `data::Workspace`: A Workspace object containing the data. Expected to have `data.X` (a DataFrame or KeyedArray of time series data) and potentially other fields.
- `groups::Workspace`: A Workspace object that defines core series. Expected to have `groups.core`, which should be an iterable (e.g., a Vector or Dict values) where each element is itself an iterable of core series names (strings).

# Returns
- Nothing. The `data.X` field is modified in-place.
"""
function normalize_noncore_series!(data::Workspace, groups::Workspace)
    all_core_series = Vector{String}()
    for g ∈ values(groups.core)
        for s ∈ g
            push!(all_core_series, s)
        end
    end
    for n ∈ names(data.X)
        if n ∉ all_core_series
            data.X[!, n] = data.X[!, n] .- mean(data.X[!, n])
            data.X[!, n] = data.X[!, n] ./ sqrt(StatsBase.var(data.X[!, n]))
        end
    end
end

"""
    get_data_obj(m, df, core=false, excluded_series=[])

Prepare a data Workspace object for a given model and DataFrame.

This function constructs a `Workspace` containing data relevant to the model `m`.
It selects series from the input `df` based on model type ("SW" or "SWFF") and
whether only `core` series are requested. Non-core series are normalized.
The function also sets up mappings between data series names and model symbols.

# Arguments
- `m`: The model object (`MacroModelling.ℳ`). Its `model_name` attribute determines the set of core series and variable groups. It also provides `m.var` (state variables) and `m.exo` (shock variables).
- `df::DataFrame`: The input DataFrame containing all available time series data.
- `core=false`: If `true`, only the predefined core series for the model `m` are selected. If `false`, all series from `df` that are also in `fred_codes` (and not in `excluded_series`) are used.
- `excluded_series=[]`: A list of series names (strings) to explicitly exclude from the dataset.

# Returns
- `data::Workspace`: A Workspace object with the following fields:
    - `:core_series::Vector{String}`: List of core series names.
    - `:core_series_modelsyms::Vector{Symbol}`: Corresponding model symbols for core series.
    - `:X::DataFrame`: The selected and potentially normalized data series.
    - `:state_vars::Vector{Symbol}`: State variables from the model.
    - `:shock_vars::Vector{Symbol}`: Shock variables from the model.
    - `:X_core::DataFrame`: A DataFrame containing only the core series, with columns renamed to model symbols.
    - `:nObs::Int`: Number of observations (rows) in `data.X`.
    - `:nStates::Int`: Number of state variables.
    - `:nVars::Int`: Number of selected data series (columns) in `data.X`.

# Details
- Core series and their corresponding model symbols are predefined based on `m.model_name`.
- Non-core series are normalized (demeaned and standardized) using `normalize_noncore_series!`.
- Relies on global `fred_codes`, `fred_groups_sw`, and `fred_groups` (the latter two selected based on `m.model_name`).
"""
function get_data_obj(m, df, core=false, excluded_series=[])
    core_series = ["RGDP", "PGDP", "RCONS", "RINV", "RWAGE", "HOURS", "FedFunds"]
    core_series_modelsyms = [:Y, :π, :C, :I, :w, :L, :R]
    groups = fred_groups_sw
    if m.model_name == "SWFF"
        core_series = ["RGDP", "PGDP", "RCONS", "RINV", "RWAGE", "HOURS", "FedFunds", "SFYBAAC"]
        core_series_modelsyms = [:Y, :π, :C, :I, :w, :L, :R, :S]
        groups = fred_groups
    end
    
    
    selected_data_series=  intersect(collect(names(df)),keys(fred_codes));
    selected_data_series = filter(x -> x ∉ excluded_series, selected_data_series)
    if core
        selected_data_series=  core_series;
    end

    
    data = Workspace(
        :core_series => core_series,
        :core_series_modelsyms => core_series_modelsyms,
        :X => df[1:end, selected_data_series],
        :state_vars => m.var,
        :shock_vars => m.exo
    );
    # TODO: don't normalize pseudo-core
    normalize_noncore_series!(data,groups);
    data.X_core = data.X[!, data.core_series];
    rename!(data.X_core, [str => sym for (str, sym) ∈ zip(data.core_series, data.core_series_modelsyms)]);
    data.nObs = size(data.X,1);
    data.nStates = length(data.state_vars);
    data.nVars = size(data.X,2);
    return data
end


"""
    load_and_treat_data(return_rng, treatment_rng=return_rng, fred_groups::Workspace=..., data_path::String=...)

Loads and processes multiple economic time series from a JLD2 file.

This function acts as a wrapper that loads a raw dataset, iterates through a list of series from the global `fred_codes`, 
and applies transformations using `transform_fred_data`. It reorders the final series to place `core` variables first and 
returns the data as both a DataFrame and an MVTSeries.

# Arguments
- `return_rng`: The `MIT` range for the final output data.
- `treatment_rng=return_rng`: The `MIT` range used for calculating transformation statistics (e.g., demeaning). Defaults to `return_rng`.
- `fred_groups::Workspace`: A Workspace defining series groups, used to order `core` series first.
- `data_path::String`: Path to the input JLD2 data file.

# Returns
- `fred_df::DataFrame`: A DataFrame containing the transformed time series.
- `fred_mvts::MVTSeries`: An `MVTSeries` object with the same transformed data.
"""
function load_and_treat_data(return_rng, treatment_rng=return_rng, fred_groups::Workspace=Workspace(; :core => Workspace()), data_path::String="$(pwd())/data/input/fred_data.jld2")
    fred_data = load(data_path)
    pop_data = transform_fred_data("POP", fred_data; assume_no_split = true)
    fred_db = Workspace()
    skipped_keys = []
    keys_to_skip = ["SP500","DJIA", "NAPMI", "NAPM_NEW_ORDERS", "NAPM_INVENTORIES", "NAPM_SUP_DEL", "PMP", "MAN_NEW_ORDERS", "INVENT"]
    for key in keys(fred_codes)
        if key ∈ keys_to_skip
            push!(skipped_keys, key)
        else
            try
                fred_db[Symbol(key)] = transform_fred_data(key, fred_data, return_rng, pop_data, treatment_rng=treatment_rng)
            catch e
                println(e)
                @warn "$(typeof(e)) for $(key) with range $(rangeof(e.a))"
            end
        end
    end
    
    # reorganize
    added_series = Set{Symbol}()
    ordered_fred_db = Workspace()
    for key in keys(fred_groups.core)
        for series ∈ fred_groups.core[key]
            ordered_fred_db[Symbol(series)] = copy(fred_db[Symbol(series)])
            push!(added_series, Symbol(series))
        end
    end
    for key ∈ keys(fred_db)
        if key ∉ added_series
            ordered_fred_db[key] = copy(fred_db[key])
        end
    end


    local fred_mvts = MVTSeries(;pairs(ordered_fred_db)...)[return_rng]
    if !isempty(skipped_keys)
        @warn "Skipped data series: $skipped_keys"
    end
 
    local fred_df = DataFrame(fred_mvts.values, collect(colnames(fred_mvts)))

    return fred_df, fred_mvts
end