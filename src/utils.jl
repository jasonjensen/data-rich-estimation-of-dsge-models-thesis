

"""
    dat2mit(f::Type{<:Frequency}, d::Date) -> MIT

Convert a `Date` object to an `MIT` (moment in time) object of a specified frequency `f`.

This function uses `TimeSeriesEcon.fconvert` to perform the conversion.
It first converts the `Date` `d` to a `BDaily` `MIT` object using `daily(d)`
and then converts this `BDaily` `MIT` to the target frequency `f`.

# Arguments
- `f::Type{<:Frequency}`: The target `TimeSeriesEcon` frequency type (e.g., `Quarterly`, `Monthly`, `BDaily`).
- `d::Date`: The `Date` object to convert.

# Returns
- `MIT`: An `MIT` object representing the same moment in time as `d`, but at frequency `f`.

# Example
```julia
using TimeSeriesEcon, Dates
my_date = Date(2023, 10, 26)
monthly_mit = dat2mit(Monthly, my_date) # Converts to 2023M10
quarterly_mit = dat2mit(Quarterly, my_date) # Converts to 2023Q4
```
"""
dat2mit(f, d::Date) = fconvert(f, daily(d))

"""
    infer_frequency(df::DataFrame) -> Type{<:Frequency}

Infer the `TimeSeriesEcon` frequency of time series data stored in a `DataFrame`.

The function inspects the `date` column of the input `DataFrame` to determine the
most likely frequency of the data. It calculates the number of calendar days between
consecutive dates (or a slightly longer span for daily/business daily) to make this inference.

Heuristics used:
- If days between first two dates <= 5:
    - If DataFrame has < 8 rows: Assumes `BDaily`.
    - Else, checks days between first and 8th date:
        - If 7 days: Assumes `Daily`.
        - Else: Assumes `BDaily`.
- Else if days between first two dates <= 10: Assumes `Weekly`. (Note: Specific weekly type like `Weekly começando domingo` is not distinguished).
- Else if days between first two dates <= 40: Assumes `Monthly`.
- Else if days between first two dates <= 100: Assumes `Quarterly`.
- Else if days between first two dates <= 200: Assumes `HalfYearly`.
- Else if days between first two dates <= 400: Assumes `Yearly`.

If none of these conditions are met, the function might not return a frequency, or its behavior for higher frequencies is undefined (as indicated by "TODO: handle more frequencies").

# Arguments
- `df::DataFrame`: A `DataFrame` containing a `date` column with `Date` objects.
                   It's assumed the dates are sorted chronologically.

# Returns
- `Type{<:Frequency}`: The inferred `TimeSeriesEcon` frequency type (e.g., `Daily`, `Monthly`, `Quarterly`).

# Notes
- The function relies on the first few dates to make its inference.
- It has TODOs for handling more specific weekly frequencies and potentially other less common frequencies.
"""
function infer_frequency(df::DataFrame)
    start_date = df[1, :date]
    next_date = df[2, :date]
    nDays = length(start_date:Day(1):next_date)
    # @show nDays
    if nDays <= 5
        size(df,1) < 8 && return BDaily
        seventh_date = df[8, :date]
        nDays7 = length(start_date:Day(1):seventh_date)
        nDays7 == 7 && return Daily
        return BDaily
    end
    nDays <= 10 && return Weekly
    # TODO: which weekly?
    nDays <= 40 && return Monthly
    nDays <= 100 && return Quarterly
    nDays <= 200 && return HalfYearly
    nDays <= 400 && return Yearly

    # TODO: handle more frequencies
end

"""
    TimeSeriesEcon.MVTSeries(x::KeyedArray, firstdate=1U)

Convert the state output from the `generate_states` function to an MVTSeries.
"""
TimeSeriesEcon.MVTSeries(x::KeyedArray, firstdate=1U) = MVTSeries(firstdate, x.Variable, Matrix(x'))

"""
    make_psd!(A::SubArray)
    make_psd(A::SubArray) -> Symmetric
    make_psd(A::Matrix) -> Symmetric
    make_psd(A_sym::Symmetric) -> Symmetric

Ensure a matrix is positive semi-definite (PSD).

This function takes a square matrix `A` (or a `SubArray` view of a matrix) and
modifies it (for `make_psd!`) or returns a new PSD matrix (for `make_psd`)
by adjusting its eigenvalues.

The process involves:
1. Ensuring the input matrix is symmetric. For `make_psd`, a symmetric copy is made if needed.
2. Performing an eigen decomposition of the symmetric matrix.
3. Replacing any negative eigenvalues with zero.
4. Reconstructing the matrix using the original eigenvectors and the modified (non-negative) eigenvalues.
5. Forcing symmetry on the reconstructed matrix by averaging it with its transpose.
6. Adding a small positive constant (e.g., `1e-6` or `1e-9`) to the diagonal elements (i.e., adding `ridge*I`)
   to ensure strict positive definiteness and improve numerical stability, effectively making it positive definite.

The `make_psd!` function modifies the input `SubArray` in-place.
The `make_psd` functions return a new `Symmetric` PSD matrix.

# Arguments
- `A::SubArray`: A subarray view of a matrix (for `make_psd!` and one version of `make_psd`).
- `A::Matrix`: A matrix (for one version of `make_psd`).
- `A_sym::Symmetric`: A symmetric matrix (for one version of `make_psd`).

# Returns
- For `make_psd!`: Modifies `A` in-place.
- For `make_psd`: Returns a new `Symmetric` matrix that is positive semi-definite (practically positive definite due to the ridge).

# Notes
- The addition of `ridge*I` makes the matrix strictly positive definite, not just PSD.
- This is a common technique to handle numerical issues with covariance matrices that should be PSD but might not be due to floating-point inaccuracies.
"""
function make_psd!(A::SubArray)
    @timeit to "make_psd" begin
        A_sym = Symmetric(copy(A))         # Ensure symmetry
        F = eigen(A_sym)             # Eigen decomposition
        # display(F)
        λ = max.(F.values, 0)        # Replace negative eigenvalues with 0
        A_psd = F.vectors * Diagonal(λ) * F.vectors'
        A_psd = (A_psd + A_psd') / 2  # Force symmetry post-reconstruction
        A_psd = A_psd + 1e-9*I # add perturbation
    end
    copyto!(A, A_psd)
end
function make_psd(A::SubArray)
    @timeit to "make_psd" begin
        A_sym = Symmetric(A)         # Ensure symmetry
        F = eigen(A_sym)             # Eigen decomposition
        # display(F)
        λ = max.(F.values, 0)        # Replace negative eigenvalues with 0
        A_psd = F.vectors * Diagonal(λ) * F.vectors'
        A_psd = (A_psd + A_psd') / 2  # Force symmetry post-reconstruction
        A_psd = A_psd + 1e-6*I # add perturbation
    end
    return Symmetric(A_psd)
end
function make_psd(A::Matrix)
    @timeit to "make_psd" begin
        A_sym = Symmetric(A)         # Ensure symmetry
        F = eigen(A_sym)             # Eigen decomposition
        # display(F)
        λ = max.(F.values, 0)        # Replace negative eigenvalues with 0
        A_psd = F.vectors * Diagonal(λ) * F.vectors'
        A_psd = (A_psd + A_psd') / 2  # Force symmetry post-reconstruction
        A_psd = A_psd + 1e-6*I # add perturbation
    end
    return Symmetric(A_psd)
end
function make_psd(A_sym::Symmetric)
    @timeit to "make_psd" begin
        F = eigen(A_sym)             # Eigen decomposition
        # display(F)
        λ = max.(F.values, 0)        # Replace negative eigenvalues with 0
        A_psd = F.vectors * Diagonal(λ) * F.vectors'
        A_psd = (A_psd + A_psd') / 2  # Force symmetry post-reconstruction
        A_psd = A_psd + 1e-6*I # add perturbation
    end
    return Symmetric(A_psd)
end



parseparams(m, θ) = Workspace(k => v for (k,v) ∈ zip(m.parameters, θ))


"""
    plot_irfs(m::MacroModelling.ℳ, θ_baseline, θ_new, shk, l=20; normalize=false, negative=false)
    plot_irfs(m::MacroModelling.ℳ, θ_vals::Vector, labels::Vector, shk, l=20; normalize=false, negative=false)

Generate and plot impulse response functions (IRFs) for a model under different parameter sets.

This function simulates the model's response to a one-standard-deviation shock (`shk`)
for a specified number of periods (`l`). It can compare IRFs from two parameter sets
(`θ_baseline`, `θ_new`) or plot IRFs for multiple parameter sets (`θ_vals`) with corresponding `labels`.

The process for each parameter set involves:
1. Solving the model using the given parameters to get the state-space representation (G, H matrices).
2. Simulating the model forward for `l` periods, starting from steady state, with a one-unit shock to `shk` in the first period.
   - If `normalize=true`, the shock size is `1 / σ_shk`, where `σ_shk` is the standard deviation of the shock `shk` (parameter `σ{shk_name_without_ϵ}`).
   - If `negative=true`, the shock is applied with a negative sign.
3. Plotting the responses of all non-shock state variables.

The first signature plots IRFs for a baseline and a new parameter set on the same axes for each variable.
The second signature plots IRFs for multiple parameter sets, each with its own color, on the same axes.
A separate legend plot is also generated in the second signature.

# Arguments
- `m::MacroModelling.ℳ`: The model object.
- `θ_baseline`: The baseline parameter vector (for the first signature).
- `θ_new`: The new parameter vector to compare against the baseline (for the first signature).
- `θ_vals::Vector`: A vector of parameter vectors (for the second signature).
- `labels::Vector`: A vector of strings corresponding to each parameter set in `θ_vals`, used for plot labels (for the second signature).
- `shk::Symbol`: The symbol of the shock to impulse (e.g., `:ϵr`).
- `l=20`: The length of the impulse response horizon in periods.
- `normalize=false`: If `true`, normalizes the shock size to one standard deviation.
- `negative=false`: If `true`, applies a negative shock.

# Returns
- `all_plots::Plots.Plot`: A combined plot object containing a grid of IRFs for each variable,
                           and a legend plot if the second signature is used. The suptitle indicates the shock.

# Errors
- Throws an error if the model does not solve for any of the provided parameter sets.
"""
function plot_irfs(m::MacroModelling.ℳ, θ_baseline, θ_new, shk, l=20; normalize=false, negative=false)
    TT, SS_and_pars, 𝐒_baseline, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ_baseline, m)
    TT, SS_and_pars, 𝐒_new, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ_new, m)
    if !solved
        error("Model not solving!")
    end

    G_keyed_baseline, H_keyed_baseline, SS = get_sparse_solution(m, 𝐒_baseline)
    G_baseline = Matrix(G_keyed_baseline)
    H_baseline = Matrix(H_keyed_baseline)

    G_keyed_new, H_keyed_new, SS = get_sparse_solution(m, 𝐒_new)
    G_new = Matrix(G_keyed_new)
    H_new = Matrix(H_keyed_new)


    vars = filter(v -> !in(v, m.exo), G_keyed_baseline.Variables)

    init_states = zeros(size(G_baseline,1))
    init_shks = zeros(size(H_baseline,2))
    shk_idx = findfirst(s -> s == shk, m.exo)
    init_shks[shk_idx] = 1.0
    

    states_baseline = zeros(size(G_baseline,1), l)
    states_new = zeros(size(G_new,1), l)
    for i = 1:l
        if i == 1
            states_baseline[:,i] = G_baseline*init_states + H_baseline*init_shks
            states_new[:,i] = G_new*init_states + H_new*init_shks
        else
            states_baseline[:,i] = G_baseline*states_baseline[:,i-1]
            states_new[:,i] = G_new*states_new[:,i-1]
        end
    end

    shks = Symbol.(replace.(string.(m.exo), "ϵ"=>"ε"))
    vars = filter(v -> v ∉ shks, G_keyed_baseline.Variables)

    q_plots = Vector{Any}()
    for var ∈ vars
        var_idx = findfirst(v -> v == var, G_keyed_baseline.Variables)
        p = Plots.plot(states_baseline[var_idx,:], label = string(var), title=string(var))
        Plots.plot!(p, states_new[var_idx,:], label = "$var, new", title=string(var))
        hline!(p, [0], legend=false)
        push!(q_plots, p)
    end

    all_plots = Plots.plot(q_plots..., layout=10, size=(1000,800), suptitle="Response to a 1 std.dev. $(shk)-shock")
    return all_plots
end
function plot_irfs(m::MacroModelling.ℳ, θ_vals, labels, shk, l=20; normalize=false, negative=false)
    default_colours = get_color_palette(:auto, plot_color(:white))
    q_plots = Vector{Any}()
    for (i, θ, label) ∈ zip(1:length(θ_vals), θ_vals, labels)
        TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
        if !solved
            error("Model $(i) not solving!")
        end
        G_keyed, H_keyed, SS = get_sparse_solution(m, 𝐒)
        G = Matrix(G_keyed)
        H = Matrix(H_keyed)
        init_states = zeros(size(G,1))
        init_shks = zeros(size(H,2))
        shk_idx = findfirst(s -> s == shk, m.exo)
        init_shks[shk_idx] = 1.0
        if normalize
            σ = Symbol(replace(string.(shk), "ϵ"=>"σ"))
            println(σ)
            init_shks[shk_idx] = 1.0 / θ[findfirst(x -> x == σ, m.parameters)]
        end
        if negative
            init_shks[shk_idx] = -init_shks[shk_idx]
        end

        states = zeros(size(G,1), l)
        for j = 1:l
            if j == 1
                states[:,j] = G*init_states + H*init_shks
            else
                states[:,j] = G*states[:,j-1]
            end
        end
        shks = Symbol.(replace.(string.(m.exo), "ϵ"=>"ε"))
        vars = filter(v -> v ∉ shks, G_keyed.Variables)

        for (v, var) ∈ enumerate(vars)
            var_idx = findfirst(v -> v == var, G_keyed.Variables)
            if i == 1
                p = Plots.plot(states[var_idx,:], label = string(var), title=string(var), color=default_colours[i])
                hline!(p, [0], legend=false, color="black")
                push!(q_plots, p)
            else
                Plots.plot!(q_plots[v], states[var_idx,:], label = "$var, new", title=string(var), color=default_colours[i])
            end
        end
        if i == 1
            p = Plots.plot(states[1,:], label = label, title="Legend", color=default_colours[i])
            # hline!(p, [0], color="black", label="zero")
            push!(q_plots, p)
        else
            Plots.plot!(q_plots[end], states[1,:], label = label, title="Legend", color=default_colours[i])
        end
        
    end
   
    all_plots = Plots.plot(q_plots..., layout=length(q_plots), size=(1000,800), suptitle="Response to a 1 std.dev. $(shk)-shock",)
    return all_plots
end

"""
    states2mvts(S::KeyedArray, last_period) -> MVTSeries

Convert a `KeyedArray` of state trajectories into an `MVTSeries`.

This function takes a `KeyedArray` `S`, where rows typically represent state variables
(indexed by `S.Variable`) and columns represent time periods, and transforms it into
an `MVTSeries` object from the `TimeSeriesEcon` package.

The time dimension of the output `MVTSeries` is determined by `last_period` and the
number of columns (time periods) in `S`. The first period of the `MVTSeries` will be
`last_period - size(S,2) + 1`. The columns of the `MVTSeries` will correspond to the
state variables `S.Variable`, and the values will be the transpose of the matrix data in `S`.

# Arguments
- `S::KeyedArray`: A KeyedArray representing state trajectories. Expected to have dimensions
                   `(num_variables, num_periods)` and row keys `S.Variable`.
- `last_period`: The time index (e.g., an `MIT` instance) corresponding to the last
                 column (last time period) in the `S` array.

# Returns
- `mvts::MVTSeries`: An `MVTSeries` object where columns are state variables and rows
                     are time periods, with the time range inferred from `last_period`
                     and the dimensions of `S`.
"""
function states2mvts(S, last_period)
    first_period = last_period - size(S,2) + 1
    mvts = MVTSeries(first_period, S.Variable, Matrix(S)')
end



function plot_irfs_sse(m::MacroModelling.ℳ, msse, θ, shk, l=20)
    TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    if !solved
        error("Model not solving!")
    end

    msse = deepcopy(sw07)
    params = parseparams(m, θ)
    for key ∈ keys(params)
        msse.parameters[key] = params[key]
    end
    StateSpaceEcon.clear_sstate!(msse)
    StateSpaceEcon.sssolve!(msse)

    pl = StateSpaceEcon.Plan(msse, 1U:(1U+l-1))
    sim_data = StateSpaceEcon.steadystatedata(msse, pl)

    sim_data[1U, shk] = 1.0

    data_sim = StateSpaceEcon.simulate(msse, pl, sim_data, deviation=true, anticipate=true)

    shks = Symbol.(replace.(string.(msse.shocks), "ϵ"=>"ε"))
    vars = filter(v -> v ∉ shks, Symbol.(msse.variables))

    q_plots = Vector{Any}()
    for var ∈ vars
        p = Plots.plot(data_sim[var], label = string(var), title=string(var), trange=1U:(1U+l-1))
        hline!(p, [0], legend=false)
        push!(q_plots, p)
    end

    all_plots = Plots.plot(q_plots..., layout=10, size=(1000,800), suptitle="Response to a 1 std.dev. $(shk)-shock")
    return all_plots
end

# plot_irfs_sse(m, sw07, m.parameter_values, :ϵr)
# plot_irfs(m, m.parameter_values, :ϵr)