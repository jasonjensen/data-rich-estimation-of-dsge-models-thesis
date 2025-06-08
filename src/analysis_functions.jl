# plotting parameters
"""
    plot_param_progress(chain::Workspace, m::MacroModelling.ℳ, targets::Dict{Symbol, Vector{Real}}; plot_step=100, suptitle="Evolution of parameters", burnin=75000, plot_size=(1000,1000))

Plot the evolution of model parameters over the course of the estimation chain.

# Arguments
- `chain::Workspace`: Workspace containing the estimation results.
- `m::MacroModelling.ℳ`: Model object.
- `targets::Dict{Symbol, Vector{Real}}`: Dictionary of target values for each parameter.
- `plot_step::Int`: Step size for plotting (default: 100).
- `suptitle::String`: Title for the entire plot (default: "Evolution of parameters").
- `burnin::Int`: Number of initial samples to discard (default: 75000).
- `plot_size::Tuple{Int, Int}`: Size of the plot (default: (1000,1000)).

# Returns
- `p_param_progress::Plots.Plot`: Plot object containing the parameter evolution plots.
"""
function plot_param_progress(chain::Workspace, m::MacroModelling.ℳ, targets::Dict{Symbol, Vector{Real}}; plot_step=100, suptitle="Evolution of parameters", burnin=75000, plot_size=(1000,1000))
    default_colours = get_color_palette(:auto, plot_color(:white))
    plots_vector = Vector{Any}()
    num_params = 29
    if m.model_name == "SWFF"
        num_params = 31
    end

    for i ∈ 1:num_params
        p = vspan([size(chain.𝛉,1) / 2 / plot_step, size(chain.𝛉,1) / plot_step], alpha=0.3, color="grey")
        title = m.parameters[i]
        if title == :Φ
            title = "S''"
        end
        param_mean = mean(chain.𝛉[burnin+1:end, i])
        param_stddev = sqrt(StatsBase.var(chain.𝛉[burnin+1:end, i]))
        param_mean_line = vcat(repeat([NaN], Int(burnin/plot_step)), repeat([param_mean], Int((size(chain.𝛉,1)-burnin)/plot_step)))
        param_stddev_line = vcat(repeat([NaN], Int(burnin/plot_step)), repeat([param_stddev], Int((size(chain.𝛉,1)-burnin)/plot_step)))
        # @show param_mean_line
        p = Plots.plot(chain.𝛉[1:plot_step:end, i], title=title, titlefont=font(10,"Bookman Light"), color=default_colours[1])
        Plots.plot!(p, param_mean_line, color="#2F4F4F", width=3, ribbon=param_stddev_line)
        # println(targets[m.parameters[i]][1])
        Plots.hline!(p, [targets[m.parameters[i]][1]], color=default_colours[2], width=3)
        # Plots.vline!(p, [1200])
        # println([size(chain.𝛉,1) / 2 / plot_step, size(chain.𝛉,1) / plot_step])
        push!(plots_vector, p)
    end

    p_param_progress = Plots.plot(plots_vector..., 
        layout=length(plots_vector), size=plot_size, legend=false,
        suptitle=suptitle,  plot_titlefontfamily="Bookman Light")
    return p_param_progress
end

# generating forecasts
"""
    zlb_objective(α::Vector, G::Matrix, H::Matrix, v::Vector, S::Vector, idx::Int64, zlb_level::Float64)

Objective function for optimizing the shock multiplier in the presence of a zero lower bound (ZLB) constraint.

# Arguments
- `α::Vector`: Vector of shock multipliers.
- `G::Matrix`: State transition matrix.
- `H::Matrix`: Shock impact matrix.
- `v::Vector`: Vector of shocks.
- `S::Vector`: Vector of state variables.
- `idx::Int64`: Index of the shock to be optimized.
- `zlb_level::Float64`: ZLB level.

# Returns
- `R::Float64`: Objective value.
"""
function zlb_objective(α::Vector, G::Matrix, H::Matrix, v::Vector, S::Vector, idx::Int64, zlb_level::Float64)
    shadow = zeros(8)
    shadow[idx] = 1.0  # Julia uses 1-based indexing
    S_upd_temp = G*S + H*(v + α[1]*shadow)
    R = abs(S_upd_temp[7] - zlb_level) # subtracting a negative number
    return R
end

# Optimization setup
"""
    optimize_zlb(G, H, v, S; α₀=[0.0], method=Brent(), shk_idx=6, zlb_level=0.0)

Optimize the shock multiplier in the presence of a zero lower bound (ZLB) constraint.

# Arguments
- `G`: State transition matrix.
- `H`: Shock impact matrix.
- `v`: Vector of shocks.
- `S`: Vector of state variables.
- `α₀`: Initial guess for the shock multiplier (default: [0.0]).
- `method`: Optimization method (default: Brent()).
- `shk_idx`: Index of the shock to be optimized (default: 6).
- `zlb_level`: ZLB level (default: 0.0).

# Returns
- `optimal_α::Float64`: Optimized shock multiplier.
- `min_R::Float64`: Minimum objective value.
"""
function optimize_zlb(G, H, v, S; α₀=[0.0], method=Brent(), shk_idx=6, zlb_level=0.0)
    # Create objective function with fixed parameters
    obj(α) = zlb_objective(α, G, H, v, S, shk_idx, zlb_level)
    
    # Non-negativity constraints (α ≥ 0)
    lx = zeros(length(α₀))  # Lower bounds
    ux = fill(Inf, length(α₀))  # Upper bounds
    
    # Create constraints object
    constraints = TwiceDifferentiableConstraints(lx, ux)


    result = optimize(obj, constraints, α₀, IPNewton(), Optim.Options(iterations=1000),
                      autodiff=:forward)
                  
    
    return Optim.minimizer(result)[1], Optim.minimum(result)
end

"""
    forecast_by_endhist!(forecasts, target_series, draw_num, m, df, endhist, θ, R, 𝛙, 𝚲, n_shock_draws, hrz; speedup=false, data_X_start=1986Q2, core=false, zlb_level=0.0, excluded_series=[])

Generate forecasts for a given end history and update the forecasts matrix.

# Arguments
- `forecasts`: Matrix to store the forecasts.
- `target_series`: Vector of target series indices.
- `draw_num`: Draw number.
- `m`: Model object.
- `df`: Data frame containing the data.
- `endhist`: End history date.
- `θ`: Vector of model parameters.
- `R`: Measurement error covariance matrix.
- `𝛙`: Measurement equation matrix.
- `𝚲`: Observation equation matrix.
- `n_shock_draws`: Number of shock draws.
- `hrz`: Forecast horizon.
- `speedup`: Flag to enable speedup (default: false).
- `data_X_start`: Start date for the data (default: 1986Q2).
- `core`: Flag to use core data (default: false).
- `zlb_level`: ZLB level (default: 0.0).
- `excluded_series`: Vector of excluded series indices (default: []).
"""
function forecast_by_endhist!(forecasts, target_series, draw_num, m, df, endhist, θ, R, 𝛙, 𝚲, n_shock_draws, hrz; speedup=false, data_X_start=1986Q2, core=false, zlb_level=0.0, excluded_series=[])
    TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    if !solved
        error("model doesn't solve")
    end
    G_keyed, H_keyed, SS = get_sparse_solution(m, 𝐒)
    G = Matrix(G_keyed)
    H = Matrix(H_keyed)
    
    # constrain data
    df_constrained = df[1:length(data_X_start:endhist),:]
    
    samp_data = get_data_obj(m, df_constrained, core, excluded_series)

    S, _ = generate_states(m, θ, 𝐒, samp_data, 𝚲, 𝛙, R; backwards_pass=false, speedup=speedup)  
    R_idx = findfirst(x->x==:R, S.Variable)
    
    state_errors_dist = MvNormal(zeros(length(m.exo)), I(length(m.exo)))
    measurement_errors_dist = MvNormal(zeros(size(R,1)), R)

    base_index = (draw_num - 1)*n_shock_draws
    Threads.@threads for n = 1:n_shock_draws
        #something?
        measure_lag = Vector(samp_data.X[end,:]) - 𝚲* S[:,end]; 

        # draw shocks
        v=rand(state_errors_dist, hrz+2)

        # draw measurement shocks?
        e=rand(measurement_errors_dist, hrz+2)

        # initialize forecast states
        states_forecast = copy(S(:, 1:hrz+2))
        states_forecast .= 0
        states_forecast[:,1:2] = S[:,end-1:end]

        data_forecast = copy(samp_data.X[1:hrz+2, :])
        data_forecast .= 0
        data_forecast[1:2,:] = samp_data.X[end-1:end, :]
        # zlb_level = -1.66
        # zero_lower flag
        # zero_lower=0;
        for i=3:hrz+2
            pred_states = rekey(Matrix(G) * states_forecast(time=i-1), 1 => G_keyed.Variables)
            states_forecast[:,i] = pred_states + H * v[:,i]

            # States are equal to the transformed previous states plus shocks
            # Sim_states_temp=G*Sim_states(:,i-1,jj)+H*v(i-1,:)';
            # println(states_forecast[5, i]) #R
            if states_forecast[R_idx, i] < zlb_level
                # println("zlb")
                shk_idx = findfirst(s -> s == :ϵr, m.exo)
                optimal_α, min_R = optimize_zlb(G, H, v[:,i], Vector(pred_states), α₀=[1.0], method=BFGS(), shk_idx=shk_idx, zlb_level=zlb_level)
                # @show optimal_α
                # @show min_R
                # println(optimal_α)

                # _,_,α,_,_,_,_,_ = csminwel(zlb, ones(1),ones(1,1), (G, H, v[:,i], pred_states); crit=0.005, max_iters=100, verbose=false)
                v[shk_idx,i] = v[shk_idx,i] + optimal_α
                states_forecast[:,i] = pred_states + H * v[:,i]
            end
            # special zlb logic
            # if Sim_states_temp(7)<-1.29
            #     z=@zlb;
            #     [~,~,alpha,~,~,~,~,~] = evalc('csminwel(z, 1,1,[],.005,100,{G, H, v(i-1,:), Sim_states(:,i-1,jj)})');
            # # alpha is about 0.63, they're adding it to the MP shock. I.e. a tightening
            #     v(i-1,6)=v(i-1,6)+alpha;
            #     Sim_states(:,i,jj)=G*Sim_states(:,i-1,jj)+H*v(i-1,:)';
            #     zero_lower=zero_lower+1;
            # else
            #     Sim_states(:,i,jj)=Sim_states_temp;
            # end
            
            # measurement error
            measure = 𝛙*measure_lag + e[:,i];

            # predicted data
            data_forecast[i,:] = 𝚲 * states_forecast[:,i] + measure;
            # Sim_X(:,i,jj)=DELTA_mean*Sim_states(:,i,jj)+measure;

            # next period's lagged measurement error is this period's measurment error
            measure_lag=measure;
        end
        # println(data_forecast[!,"RGDP"])
        # println(samp_data[!,"RGDP"])
        forecasts[Symbol(endhist)][base_index + n, :] = data_forecast[!,target_series]
    end

end

"""
    gen_forecasts_for_range(m, target_series, res, df, burnin, n_parameter_draws, n_shock_draws, rng, hrz; speedup=true, data_X_start=1986Q2, core = false, zlb_level=0.0, excluded_series=[])

Generate forecasts for a range of end histories and parameter draws.

# Arguments
- `m`: Model object.
- `target_series`: Vector of target series indices.
- `res`: Estimation results.
- `df`: Data frame containing the data.
- `burnin`: Number of initial samples to discard.
- `n_parameter_draws`: Number of parameter draws.
- `n_shock_draws`: Number of shock draws.
- `rng`: Range of end histories.
- `hrz`: Forecast horizon.
- `speedup`: Flag to enable speedup (default: true).
- `data_X_start`: Start date for the data (default: 1986Q2).
- `core`: Flag to use core data (default: false).
- `zlb_level`: ZLB level (default: 0.0).
- `excluded_series`: Vector of excluded series indices (default: []).

# Returns
- `forecasts`: Workspace containing the forecasts.
"""
function gen_forecasts_for_range(m, target_series, res, df, burnin, n_parameter_draws, n_shock_draws, rng, hrz; speedup=true, data_X_start=1986Q2, core = false, zlb_level=0.0, excluded_series=[])
    forecasts = Workspace()
    draw_num = 1
    failure_count = 0
    success_count = 0
    while draw_num < n_parameter_draws
        if draw_num == 1
            for endhist ∈ rng
                forecasts[Symbol(endhist)] = zeros(n_parameter_draws*n_shock_draws,hrz+2)
            end
        end
        draw_index = rand(burnin+1:size(res.𝛉,1))
        # draw_index = rand(7500:10000)
        # draw_index = rand(1:100)
        # _θ = res.𝛉[draw_index, :]
        _θ = res.𝛉[draw_index, :]
        _R = Diagonal(res.R[draw_index, :])
        _𝛙 = Diagonal(res.𝛙[draw_index, :])
        _𝚲 = copy(res.𝚲[draw_index])

        # if draw_num < 5
        #     walk_model_to_θ(m, _θ)
        # end
        # if m.model_name == "SWFF"
        #     # _θ = copy(SWFF.parameter_values)
        #     _θ = [[gelfer_means_swff_dfm_estimation[v][1] for v ∈ m.parameters[1:31]]..., m.parameter_values[32:end]...]
        #     if speedup
        #         @info "overwriting"
        #         theirΛ = vars["DELTA_estimates"][draw_index-burnin]
        #         _𝚲 = convert_mat_Λ(_𝚲, theirΛ)
        #         tempR = res.R[draw_index, :]
        #         tempR[1:17] = matR["R_estimates"][1:17, draw_index-burnin]
        #         _R = Diagonal(tempR)

        #         temp𝛙 = res.𝛙[draw_index, :]
        #         temp𝛙[1:17] = matpsi["PSI_estimates"][1:17, draw_index-burnin]
        #         _𝛙 = Diagonal(temp𝛙)
        #     end

        # end

        # _θ[findfirst(x -> x == :ry1, m.parameters)] = 0.18
        # _θ[findfirst(x -> x == :σr, m.parameters)] = 0.125

        try
            for endhist ∈ rng
                forecast_by_endhist!(forecasts, target_series, draw_num, m, df, endhist, _θ, _R, _𝛙, _𝚲, n_shock_draws, hrz, speedup=speedup, data_X_start=data_X_start, core=core, zlb_level=zlb_level, excluded_series=excluded_series)
            end
            success_count = success_count + 1
        catch
            # TODO: walk model to solution
            failure_count = failure_count + 1
            # model didn't solve
            # walk_model_to_θ(m, _θ)
            if failure_count > 30000
                error("Maximum failure count reached. Successes: $(success_count)")
            end
            continue
        end
        draw_num = draw_num + 1
    end
    return forecasts
end

"""
    generate_and_plot_forecasts(
        target::String,
        rng::Union{UnitRange{<:MIT}, Vector{<:MIT}},
        m_sw::MacroModelling.ℳ,
        m_swff::MacroModelling.ℳ,
        df_sw::DataFrame,
        df_swff::DataFrame,
        mvts::MVTSeries,
        chain_sw_reg::Workspace,
        chain_swff_reg::Workspace,
        chain_sw_dfm::Workspace,
        chain_swff_dfm::Workspace;
        hrz=16,
        n_parameter_draws=50,
        n_shock_draws = 100,
        burnin=75000,
        data_X_start = 1986Q2,
        zlb_level=0.0,
        vartype=:growth,
        return_vector=false,
        subtitle_prefix="",
        legend=:bottomright,
        excluded_series=[],
        ylim=:auto, 
        layout=:auto,
        size=(1000,800))

Generate and plot forecasts for different models and specifications.

# Arguments
- `target::String`: Target variable to forecast.
- `rng::Union{UnitRange{<:MIT}, Vector{<:MIT}}`: Range of end histroies for forecasting.
- `m_sw::MacroModelling.ℳ`: SW model object.
- `m_swff::MacroModelling.ℳ`: SWFF model object.
- `df_sw::DataFrame`: DataFrame for SW model.
- `df_swff::DataFrame`: DataFrame for SWFF model.
- `mvts::MVTSeries`: MVTSeries object containing actual data.
- `chain_sw_reg::Workspace`: Workspace for SW model (regular estimation).
- `chain_swff_reg::Workspace`: Workspace for SWFF model (regular estimation).
- `chain_sw_dfm::Workspace`: Workspace for SW model (DFM estimation).
- `chain_swff_dfm::Workspace`: Workspace for SWFF model (DFM estimation).
- `hrz=16`: Forecast horizon.
- `n_parameter_draws=50`: Number of parameter draws.
- `n_shock_draws = 100`: Number of shock draws.
- `burnin=75000`: Number of burn-in samples.
- `data_X_start = 1986Q2`: Start date for data.
- `zlb_level=0.0`: Level of the zero lower bound.
- `vartype=:growth`: Type of variable (e.g., :growth, :level).
- `return_vector=false`: Whether to return a vector of plots.
- `subtitle_prefix=""`: Prefix for plot subtitles.
- `legend=:bottomright`: Legend position.
- `excluded_series=[]`: List of excluded series.
- `ylim=:auto`: Y-axis limits.
- `layout=:auto`: Plot layout.
- `size=(1000,800)`: Plot size.
"""

function generate_and_plot_forecasts(
    target::String,
    rng::Union{UnitRange{<:MIT}, Vector{<:MIT}},
    m_sw::MacroModelling.ℳ,
    m_swff::MacroModelling.ℳ,
    df_sw::DataFrame,
    df_swff::DataFrame,
    mvts::MVTSeries,
    chain_sw_reg::Workspace,
    chain_swff_reg::Workspace,
    chain_sw_dfm::Workspace,
    chain_swff_dfm::Workspace;
    hrz=16,
    n_parameter_draws=50,
    n_shock_draws = 100,
    burnin=75000,
    data_X_start = 1986Q2,
    zlb_level=0.0,
    vartype=:growth,
    return_vector=false,
    subtitle_prefix="",
    legend=:bottomright,
    excluded_series=[],
    ylim=:auto, 
    layout=:auto,
    size=(1000,800))

    forecasts_sw_reg =    gen_forecasts_for_range(m_sw, target, chain_sw_reg, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level, excluded_series=excluded_series)
    forecasts_swff_reg =  gen_forecasts_for_range(m_swff, target, chain_swff_reg, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level, excluded_series=excluded_series)
    forecasts_sw_dfm =    gen_forecasts_for_range(m_sw, target, chain_sw_dfm, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=true, data_X_start=data_X_start, core=false, zlb_level=zlb_level, excluded_series=excluded_series)
    forecasts_swff_dfm =  gen_forecasts_for_range(m_swff, target, chain_swff_dfm, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=true, data_X_start=data_X_start, core=false, zlb_level=zlb_level, excluded_series=excluded_series)

    q_plots = Vector{Any}()
    for mit ∈ keys(forecasts_sw_reg)
        default_colours = get_color_palette(:auto, plot_color(:white))
        # default_colours =  cgrad(:curl, 5, categorical = true)
        # default_colours =  cgrad(:matter, 4, categorical = true)
        _mit = eval(Meta.parse(string(mit)))
        if vartype == :growth
            med_sw_reg, std_sw_reg = get_growth_medians_stddev(forecasts_sw_reg[mit])
            med_sw_dfm, std_sw_dfm = get_growth_medians_stddev(forecasts_sw_dfm[mit])
            med_swff_reg, std_swff_reg = get_growth_medians_stddev(forecasts_swff_reg[mit])
            med_swff_dfm, std_swff_dfm = get_growth_medians_stddev(forecasts_swff_dfm[mit])
            plot_hrz = hrz
            divisor = 100.0
            times = string.(_mit:_mit+(plot_hrz))
        elseif vartype == :level
            med_sw_reg, std_sw_reg = get_medians_stddev(forecasts_sw_reg[mit])
            med_sw_dfm, std_sw_dfm = get_medians_stddev(forecasts_sw_dfm[mit])
            med_swff_reg, std_swff_reg = get_medians_stddev(forecasts_swff_reg[mit])
            med_swff_dfm, std_swff_dfm = get_medians_stddev(forecasts_swff_dfm[mit])
            plot_hrz = hrz
            divisor = 1.0
            times = string.(_mit:_mit+(plot_hrz))
        end
        
        
        p = Plots.plot(string.(_mit:_mit+(plot_hrz)), med_sw_dfm ./ divisor; ribbon=std_sw_dfm  ./ divisor, label="SW-Rich", xrotation = 60, color=default_colours[1], width=2, ylim=ylim)
        Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), med_swff_dfm  ./ divisor;ribbon=std_swff_dfm  ./ divisor, label="SWFF-Rich", xrotation = 60, color=default_colours[6], width=2, fillstyle=:x)
        Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), med_sw_reg  ./ divisor; label="SW", xrotation = 60, title="$(subtitle_prefix)$(string(mit))", legend=legend, titlefont=font(10,"Bookman Light"), color=default_colours[1], width=2, linestyle=:dash)
        Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), med_swff_reg  ./ divisor; label="SWFF", xrotation = 60, color=default_colours[6], width=2, linestyle=:dash)
        # p = Plots.plot(string.(_mit:_mit+(16)), y_growth_medians; ribbon=y_growth_stddev, label=string(mit), xrotation = 60)

        last_dp = min(_mit+(plot_hrz),lastdate(mvts[Symbol(target)]))
        actuals = copy(mvts[Symbol(target)])
        if last_dp < _mit+plot_hrz
            actuals[last_dp+1:_mit+plot_hrz+1] .= NaN
        end

        if vartype == :growth
            Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), values(diff(actuals)[_mit:_mit+(plot_hrz)]), label="actual", color="black", width=2, opacity=0.5)
        elseif vartype == :level
            Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), values(actuals[_mit:_mit+(plot_hrz)]), label="actual", color="black", width=2, opacity=0.5)
        end
            
        push!(q_plots, p)
    end
    if return_vector
        return q_plots
    end

    _layout = layout
    if layout == :auto
       _layout = length(q_plots)
    end

    forecasts_plot = Plots.plot(q_plots..., layout=_layout, size=size, suptitle="Q/Q $(target) Growth forecasts", plot_titlefontfamily="Bookman Light", bottom_margin=30Plots.px)
    return forecasts_plot
end


"""
    generate_and_plot_forecasts_mes(
        target::String,
        rng::UnitRange{<:MIT},
        m_sw::MacroModelling.ℳ,
        m_swff::MacroModelling.ℳ,
        df_sw::DataFrame,
        df_swff::DataFrame,
        mvts::MVTSeries,
        chain_sw_reg::Workspace,
        chain_swff_reg::Workspace,
        chain_sw_mes::Workspace,
        chain_swff_mes::Workspace;
        hrz=16,
        n_parameter_draws=50,
        n_shock_draws = 100,
        burnin=75000,
        data_X_start = 1986Q2,
        zlb_level=0.0,
        vartype=:growth,
        return_vector=false,
        subtitle_prefix="",
        legend=:bottomright,
        excluded_series=[],
        ylim=:auto
    )

Generate and plot forecasts for different models and specifications, including MES results.

# Arguments
- `target::String`: Target variable to forecast.
- `rng::UnitRange{<:MIT}`: Range of end histories for forecasting.
- `m_sw::MacroModelling.ℳ`: SW model object.
- `m_swff::MacroModelling.ℳ`: SWFF model object.
- `df_sw::DataFrame`: DataFrame for SW model.
- `df_swff::DataFrame`: DataFrame for SWFF model.
- `mvts::MVTSeries`: MVTSeries object containing actual data.
- `chain_sw_reg::Workspace`: Workspace for SW model (regular estimation).
- `chain_swff_reg::Workspace`: Workspace for SWFF model (regular estimation).
- `chain_sw_mes::Workspace`: Workspace for SW model (MES estimation).
- `chain_swff_mes::Workspace`: Workspace for SWFF model (MES estimation).
- `hrz=16`: Forecast horizon.
- `n_parameter_draws=50`: Number of parameter draws.
- `n_shock_draws = 100`: Number of shock draws.
- `burnin=75000`: Number of burn-in samples.
- `data_X_start = 1986Q2`: Start date for data.
- `zlb_level=0.0`: Level of the zero lower bound.
- `vartype=:growth`: Type of variable (e.g., :growth, :level).
- `return_vector=false`: Whether to return a vector of plots.
- `subtitle_prefix=""`: Prefix for plot subtitles.
- `legend=:bottomright`: Legend position.
- `excluded_series=[]`: List of excluded series.
- `ylim=:auto`: Y-axis limits.
"""
function generate_and_plot_forecasts_mes(
    target::String,
    rng::UnitRange{<:MIT},
    m_sw::MacroModelling.ℳ,
    m_swff::MacroModelling.ℳ,
    df_sw::DataFrame,
    df_swff::DataFrame,
    mvts::MVTSeries,
    chain_sw_reg::Workspace,
    chain_swff_reg::Workspace,
    chain_sw_mes::Workspace,
    chain_swff_mes::Workspace;
    hrz=16,
    n_parameter_draws=50,
    n_shock_draws = 100,
    burnin=75000,
    data_X_start = 1986Q2,
    zlb_level=0.0,
    vartype=:growth,
    return_vector=false,
    subtitle_prefix="",
    legend=:bottomright,
    excluded_series=[],
    ylim=:auto
    )

    forecasts_sw_reg =    gen_forecasts_for_range(m_sw, target, chain_sw_reg, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level, excluded_series=excluded_series)
    forecasts_swff_reg =  gen_forecasts_for_range(m_swff, target, chain_swff_reg, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level, excluded_series=excluded_series)
    forecasts_sw_mes =    gen_forecasts_for_range(m_sw, target, chain_sw_mes, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level, excluded_series=excluded_series)
    forecasts_swff_mes =  gen_forecasts_for_range(m_swff, target, chain_swff_mes, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, hrz, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level, excluded_series=excluded_series)

    q_plots = Vector{Any}()
    for mit ∈ keys(forecasts_sw_reg)
        default_colours = get_color_palette(:auto, plot_color(:white))
        _mit = eval(Meta.parse(string(mit)))
        if vartype == :growth
            med_sw_reg, std_sw_reg = get_growth_medians_stddev(forecasts_sw_reg[mit])
            med_sw_mes, std_sw_mes = get_growth_medians_stddev(forecasts_sw_mes[mit])
            med_swff_reg, std_swff_reg = get_growth_medians_stddev(forecasts_swff_reg[mit])
            med_swff_mes, std_swff_mes = get_growth_medians_stddev(forecasts_swff_mes[mit])
            plot_hrz = hrz
            divisor = 100.0
            times = string.(_mit:_mit+(plot_hrz))
        elseif vartype == :level
            med_sw_reg, std_sw_reg = get_medians_stddev(forecasts_sw_reg[mit])
            med_sw_mes, std_sw_mes = get_medians_stddev(forecasts_sw_mes[mit])
            med_swff_reg, std_swff_reg = get_medians_stddev(forecasts_swff_reg[mit])
            med_swff_mes, std_swff_mes = get_medians_stddev(forecasts_swff_mes[mit])
            plot_hrz = hrz
            divisor = 1.0
            times = string.(_mit:_mit+(plot_hrz))
        end
        
        
        p = Plots.plot(string.(_mit:_mit+(plot_hrz)), med_sw_mes ./ divisor; ribbon=std_sw_mes  ./ divisor, label="SW-Sparse", xrotation = 60, color=default_colours[1], width=2, ylim=ylim)
        Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), med_swff_mes  ./ divisor;ribbon=std_swff_mes  ./ divisor, label="SWFF-Sparse", xrotation = 60, color=default_colours[6], width=2, fillstyle=:x)
        Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), med_sw_reg  ./ divisor; label="SW", xrotation = 60, title="$(subtitle_prefix)$(string(mit))", legend=legend, titlefont=font(10,"Bookman Light"), color=default_colours[1], width=2, linestyle=:dash)
        Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), med_swff_reg  ./ divisor; label="SWFF", xrotation = 60, color=default_colours[6], width=2, linestyle=:dash)
        # p = Plots.plot(string.(_mit:_mit+(16)), y_growth_medians; ribbon=y_growth_stddev, label=string(mit), xrotation = 60)

        last_dp = min(_mit+(plot_hrz),lastdate(mvts[Symbol(target)]))
        actuals = copy(mvts[Symbol(target)])
        if last_dp < _mit+plot_hrz
            actuals[last_dp+1:_mit+plot_hrz+1] .= NaN
        end

        if vartype == :growth
            Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), values(diff(actuals)[_mit:_mit+(plot_hrz)]), label="actual", color="black", width=2, opacity=0.5)
        elseif vartype == :level
            Plots.plot!(p, string.(_mit:_mit+(plot_hrz)), values(actuals[_mit:_mit+(plot_hrz)]), label="actual", color="black", width=2, opacity=0.5)
        end

        push!(q_plots, p)
    end
    if return_vector
        return q_plots
    end

    

    forecasts_plot = Plots.plot(q_plots..., layout=length(keys(forecasts_sw_reg)), size=(1000,800), suptitle="Q/Q $(target) Growth forecasts", plot_titlefontfamily="Bookman Light")
    return forecasts_plot
end

"""
    generate_all_forecasts(
        target::String,
        rng::UnitRange{<:MIT},
        m_sw::MacroModelling.ℳ, 
        m_swff::MacroModelling.ℳ,
        chain_sw_reg::Workspace,
        chain_swff_reg::Workspace,
        chain_sw_dfm::Workspace,
        chain_swff_dfm::Workspace,
        chain_sw_mes::Workspace,
        chain_swff_mes::Workspace,
        df_sw::DataFrame,
        df_swff::DataFrame;
        horizon = 4,
        burnin=75000,
        n_parameter_draws=50,
        n_shock_draws = 100,
        data_X_start = 1986Q2,
        zlb_level=0.0,
        vartype=:growth,
        excluded_series=[]
    )

Generate all forecasts for a given target variable and range of end histories.

# Arguments
- `target::String`: Target variable to forecast.
- `rng::UnitRange{<:MIT}`: Range of end histories for forecasting.
- `m_sw::MacroModelling.ℳ`: SW model object.
- `m_swff::MacroModelling.ℳ`: SWFF model object.
- `chain_sw_reg::Workspace`: Workspace for SW model (regular estimation).
- `chain_swff_reg::Workspace`: Workspace for SWFF model (regular estimation).
- `chain_sw_dfm::Workspace`: Workspace for SW model (DFM estimation).
- `chain_swff_dfm::Workspace`: Workspace for SWFF model (DFM estimation).
- `chain_sw_mes::Workspace`: Workspace for SW model (MES estimation).
- `chain_swff_mes::Workspace`: Workspace for SWFF model (MES estimation).
- `df_sw::DataFrame`: DataFrame for SW model.
- `df_swff::DataFrame`: DataFrame for SWFF model.
- `horizon=4`: Forecast horizon.
- `burnin=75000`: Number of burn-in samples.
- `n_parameter_draws=50`: Number of parameter draws.
- `n_shock_draws=100`: Number of shock draws.
- `data_X_start=1986Q2`: Start date for data.
- `zlb_level=0.0`: Level of the zero lower bound.
- `vartype=:growth`: Type of variable (e.g., :growth, :level).
- `excluded_series=[]`: List of excluded series.
"""
function generate_all_forecasts(
        target::String,
        rng::UnitRange{<:MIT},
        m_sw::MacroModelling.ℳ, 
        m_swff::MacroModelling.ℳ,
        chain_sw_reg::Workspace,
        chain_swff_reg::Workspace,
        chain_sw_dfm::Workspace,
        chain_swff_dfm::Workspace,
        chain_sw_mes::Workspace,
        chain_swff_mes::Workspace,
        df_sw::DataFrame,
        df_swff::DataFrame;
        horizon = 4,
        burnin=75000,
        n_parameter_draws=50,
        n_shock_draws = 100,
        data_X_start = 1986Q2,
        zlb_level=0.0,
        vartype=:growth,
        excluded_series=[]
    )

    t() = Dates.format(Dates.Time(Dates.now()), "HH:MM:SS")
    println("$(t()) Generating $target forecasts for sw_reg")
    forecasts_sw_reg =    gen_forecasts_for_range(m_sw, target, chain_sw_reg, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, horizon, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level)
    println("$(t()) Generating $target forecasts for swff_reg")
    forecasts_swff_reg =  gen_forecasts_for_range(m_swff, target, chain_swff_reg, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, horizon, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level)
    println("$(t()) Generating $target forecasts for sw_dfm")
    forecasts_sw_dfm =    gen_forecasts_for_range(m_sw, target, chain_sw_dfm, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, horizon, speedup=true, data_X_start=data_X_start, core=false, zlb_level=zlb_level, excluded_series=excluded_series)
    println("$(t()) Generating $target forecasts for swff_dfm")
    forecasts_swff_dfm =  gen_forecasts_for_range(m_swff, target, chain_swff_dfm, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, horizon, speedup=true, data_X_start=data_X_start, core=false, zlb_level=zlb_level, excluded_series=excluded_series)
    println("$(t()) Generating $target forecasts for sw_mes")
    forecasts_sw_mes =    gen_forecasts_for_range(m_sw, target, chain_sw_mes, df_sw, burnin, n_parameter_draws, n_shock_draws, rng, horizon, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level)
    println("$(t()) Generating $target forecasts for swff_mes")
    forecasts_swff_mes =  gen_forecasts_for_range(m_swff, target, chain_swff_mes, df_swff, burnin, n_parameter_draws, n_shock_draws, rng, horizon, speedup=false, data_X_start=data_X_start, core=true, zlb_level=zlb_level)
    println("$(t()) Done.")

    all_forecasts = Workspace(
        :target => target,
        :vartype => vartype,
        :excluded_series => excluded_series,
        :forecasts_sw_reg => forecasts_sw_reg,
        :forecasts_swff_reg => forecasts_swff_reg,
        :forecasts_sw_dfm => forecasts_sw_dfm,
        :forecasts_swff_dfm => forecasts_swff_dfm,
        :forecasts_sw_mes => forecasts_sw_mes,
        :forecasts_swff_mes => forecasts_swff_mes,
    )

    
    return all_forecasts
end



"""
    make_comparison_chart(
        target::String,
        rng::UnitRange{<:MIT},
        m_sw::MacroModelling.ℳ,
        m_swff::MacroModelling.ℳ,
        bundle_07full::Workspace,
        bundle_07part::Workspace;
        hrz=16,
        n_parameter_draws=50,
        n_shock_draws = 100,
        burnin=75000,
        data_X_start = 1986Q2,
        vartype=:growth,
        legend=:bottomright,
        ylim=:auto)

Create a comparison chart of forecasts from different models and data bundles.

# Arguments
- `target::String`: Target variable for forecasting.
- `rng::UnitRange{<:MIT}`: Range of end histories for forecasting.
- `m_sw::MacroModelling.ℳ`: SW model object.
- `m_swff::MacroModelling.ℳ`: SWFF model object.
- `bundle_07full::Workspace`: Workspace containing full 2007 data and results.
- `bundle_07part::Workspace`: Workspace containing partial 2007 data and results.
- `hrz=16`: Forecast horizon.
- `n_parameter_draws=50`: Number of parameter draws.
- `n_shock_draws=100`: Number of shock draws.
- `burnin=75000`: Number of burn-in samples.
- `data_X_start=1986Q2`: Start date for data.
- `vartype=:growth`: Type of variable (e.g., :growth, :level).
- `legend=:bottomright`: Legend position.
- `ylim=:auto`: Y-axis limits.
"""
function make_comparison_chart(
    target::String,
    rng::UnitRange{<:MIT},
    m_sw::MacroModelling.ℳ,
    m_swff::MacroModelling.ℳ,
    bundle_07full::Workspace,
    bundle_07part::Workspace;
    hrz=16,
    n_parameter_draws=50,
    n_shock_draws = 100,
    burnin=75000,
    data_X_start = 1986Q2,
    vartype=:growth,
    legend=:bottomright,
    ylim=:auto)

    
    forecasts_07full = generate_and_plot_forecasts(target, rng, m_sw, m_swff, bundle_07full.df_sw, bundle_07full.df_swff, bundle_07full.mvts_sw, bundle_07full.chain_sw_reg, bundle_07full.chain_swff_reg, bundle_07full.chain_sw_dfm, bundle_07full.chain_swff_dfm; n_parameter_draws=n_parameter_draws, n_shock_draws=n_shock_draws, zlb_level=bundle_07full.zlb_level, burnin=burnin, hrz=hrz, data_X_start=data_X_start, vartype=vartype, subtitle_prefix="Full, ", legend=legend, ylim=ylim, return_vector=true)

    forecasts_07mes = generate_and_plot_forecasts_mes(target, rng, m_sw, m_swff, bundle_07full.df_sw, bundle_07full.df_swff, bundle_07full.mvts_sw, bundle_07full.chain_sw_reg, bundle_07full.chain_swff_reg, bundle_07full.chain_sw_mes, bundle_07full.chain_swff_mes; n_parameter_draws=n_parameter_draws, n_shock_draws=n_shock_draws, zlb_level=bundle_07full.zlb_level, burnin=burnin, hrz=hrz, data_X_start=data_X_start, vartype=vartype, subtitle_prefix="Measurement Error, ", legend=legend, ylim=ylim, return_vector=true)

    forecasts_07part = generate_and_plot_forecasts(target, rng, m_sw, m_swff, bundle_07part.df_sw, bundle_07part.df_swff, bundle_07part.mvts_sw, bundle_07part.chain_sw_reg, bundle_07part.chain_swff_reg, bundle_07part.chain_sw_dfm, bundle_07part.chain_swff_dfm; n_parameter_draws=n_parameter_draws, n_shock_draws=n_shock_draws, zlb_level=bundle_07part.zlb_level, burnin=burnin, hrz=hrz, data_X_start=data_X_start, vartype=vartype, subtitle_prefix="Fewer Series, ", legend=legend, ylim=ylim, excluded_series=bundle_07part.excluded_series, return_vector=true, )

    q_plots = [
        forecasts_07full[1],
        forecasts_07mes[1],
        forecasts_07part[1],
        forecasts_07full[2],
        forecasts_07mes[2],
        forecasts_07part[2],
    ]

    forecasts_plot = Plots.plot(q_plots..., layout=length(q_plots), size=(1000,800), suptitle="Q/Q $(target) Growth forecasts", plot_titlefontfamily="Bookman Light", dpi=1000)
    return forecasts_plot
end

"""
    walk_model_to_θ(m, θ; steps=100)

Walk the model parameters from their original values to the target `θ` in a specified number of steps.
This is used to check if the model solves at intermediate parameter values, which can help diagnose issues when a model fails to solve with a given `θ`.

# Arguments
- `m`: Model object.
- `θ`: Target parameter vector.
- `steps=100`: Number of steps to take from original parameters to `θ`.
"""
function walk_model_to_θ(m, θ; steps=100)
    # println("walking model.")
    if m.model_name == "SW"
        θ_orig = copy(SW.parameter_values)
    elseif m.model_name == "SWFF"
        θ_orig = copy(SWFF.parameter_values)
    else
        error("unknown model: $(m.model_name)")
    end
    failed = false
    for i = 1:steps
        θ_step = θ_orig + ((i/steps) .* (θ .- θ_orig))
        TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ_step, m)
        if !solved && !failed
            println("Walk failed on step $i")
            failed = true
            break
        end
    end
end

"""
    get_growth_medians_stddev(forecasts)

Calculate the median and standard deviation of growth rates from a matrix of level forecasts.

# Arguments
- `forecasts`: A matrix where each row is a forecast draw and each column is a time period.

# Returns
- A tuple containing two vectors: the median growth rates and the standard deviation of growth rates.
"""
function get_growth_medians_stddev(forecasts)
    growth_forecasts = (forecasts[:, 2:end] .- forecasts[:, 1:end-1]) .* 100
    growth_medians = median(growth_forecasts, dims=1)[1,:]
    growth_stddev = sqrt.(StatsBase.var(growth_forecasts, dims=1))[1,:]
    return (growth_medians, growth_stddev)
end

"""
    get_growth_medians(ws::Workspace)

Calculate the median growth rates for each series in a Workspace of forecasts.

# Arguments
- `ws::Workspace`: A Workspace where keys are time periods (MIT) and values are forecast matrices.

# Returns
- `ws2::Workspace`: A new Workspace where keys are time periods and values are vectors of median growth rates.
"""
function get_growth_medians(ws::Workspace)
    ws2 = Workspace()
    for mit ∈ keys(ws)
        med, std = get_growth_medians_stddev(ws[mit])
        ws2[mit] = med
    end
    return ws2
end

"""
    get_medians_stddev(forecasts)

Calculate the median and standard deviation of level forecasts.

# Arguments
- `forecasts`: A matrix where each row is a forecast draw and each column is a time period.

# Returns
- A tuple containing two vectors: the median levels and the standard deviation of levels.
"""
function get_medians_stddev(forecasts)
    medians = median(forecasts[:, 2:end], dims=1)[1,:]
    stddev = sqrt.(StatsBase.var(forecasts, dims=1))[1,:]
    return (medians, stddev)
end

"""
    get_medians(ws::Workspace)

Calculate the median levels for each series in a Workspace of forecasts.

# Arguments
- `ws::Workspace`: A Workspace where keys are time periods (MIT) and values are forecast matrices.

# Returns
- `ws2::Workspace`: A new Workspace where keys are time periods and values are vectors of median levels.
"""
function get_medians(ws::Workspace)
    ws2 = Workspace()
    for mit ∈ keys(ws)
        med, std = get_medians_stddev(ws[mit])
        ws2[mit] = med
    end
    return ws2
end


## Diebold Mariano stuff
"""
    autocovariance(x::Vector{Float64}, lag::Int)

Compute the autocovariance of a vector `x` at a given `lag`.

# Arguments
- `x::Vector{Float64}`: Input time series.
- `lag::Int`: The lag at which to compute the autocovariance.

# Returns
- The autocovariance value.
"""
function autocovariance(x::Vector{Float64}, lag::Int)
    n = length(x)
    x_mean = mean(x)
    sum((x[1:n-lag] .- x_mean) .* (x[lag+1:n] .- x_mean)) / n
end

"""
    dm_test(d::Vector{Float64}, h::Int=1)

Perform the Diebold-Mariano test for predictive accuracy.

# Arguments
- `d::Vector{Float64}`: Vector of loss differentials between two forecasts.
- `h::Int=1`: Forecast horizon. Used to determine the number of autocovariances to include in the variance calculation.

# Returns
- The Diebold-Mariano test statistic.
"""
function dm_test(d::Vector{Float64}, h::Int=1)
    # d here is a loss differential
    T = length(d)
    d_bar = mean(d)
    
    # Variance estimation with autocovariances up to lag h-1
    var_d = autocovariance(d, 0)
    for k in 1:h-1
        var_d += 2 * autocovariance(d, k)
    end
    var_d /= T
    
    dm_stat = d_bar / sqrt(abs.(var_d))
    p_value = 2 * (1 - Distributions.cdf(Distributions.Normal(0, 1), abs(dm_stat)))
    return (dm_stat, p_value)
end

# function get_loss_growth(forecasts, hrz, actual)
#     first_mit = eval(Meta.parse(string(first(keys(forecasts)))))
#     loss_ts = TSeries(first_mit:first_mit+length(keys(forecasts)))
#     for endhist_sym ∈ keys(forecasts)
#         endhist = eval(Meta.parse(string(endhist_sym)))
#         loss_ts[endhist] = (forecasts[endhist_sym][1+hrz] - actual[endhist+hrz]) .^ 2
#     end
#     return loss_ts
# end

"""
    get_loss(forecasts, hrz, actual)

Calculate the squared error loss for a given forecast horizon.

# Arguments
- `forecasts`: Workspace of forecasts, where keys are end_hist dates (MIT) and values are forecast vectors.
- `hrz`: Forecast horizon (integer).
- `actual`: TSeries of actual values.

# Returns
- `loss_ts`: TSeries containing the squared error loss for each forecast period.
"""
function get_loss(forecasts, hrz, actual)
    first_mit = eval(Meta.parse(string(first(keys(forecasts)))))
    loss_ts = TSeries(first_mit:first_mit+length(keys(forecasts))-1)
    for endhist_sym ∈ keys(forecasts)
        endhist = eval(Meta.parse(string(endhist_sym)))
        # println("$(forecasts[endhist_sym][1+hrz]) vs )
        loss_ts[endhist] = (forecasts[endhist_sym][1+hrz] - actual[endhist+hrz]) .^ 2
    end
    return loss_ts
end

"""
    get_cumulative_loss(forecasts, hrz, actual)

Calculate the squared error of the cumulative forecast (sum of forecasts up to horizon `hrz`) against the cumulative actual values.

# Arguments
- `forecasts`: Workspace of forecasts, where keys are end_hist dates (MIT) and values are forecast vectors.
- `hrz`: Forecast horizon (integer).
- `actual`: TSeries of actual values.

# Returns
- `loss_ts`: TSeries containing the squared error of the cumulative forecast for each forecast period.
"""
function get_cumulative_loss(forecasts, hrz, actual)
    first_mit = eval(Meta.parse(string(first(keys(forecasts)))))
    loss_ts = TSeries(first_mit:first_mit+length(keys(forecasts))-1)
    for endhist_sym ∈ keys(forecasts)
        endhist = eval(Meta.parse(string(endhist_sym)))
        # println("$(forecasts[endhist_sym][1+hrz]) vs )
        loss_ts[endhist] = (sum(forecasts[endhist_sym][1:1+hrz]) - sum(actual[endhist:endhist+hrz])) .^ 2
    end
    return loss_ts
end

"""
    print_dm_table(medians::Workspace, actual::TSeries, ranges)

Print a LaTeX-formatted table of Diebold-Mariano test statistics comparing different models and forecast horizons.

# Arguments
- `medians::Workspace`: Workspace containing median forecasts from different models (e.g., `medians.sw_reg`, `medians.sw_dfm`).
- `actual::TSeries`: TSeries of actual values.
- `ranges`: A vector of three date ranges (MIT) over which to calculate the DM statistics.
"""
function print_dm_table(medians::Workspace, actual::TSeries, ranges)
    if length(ranges) !== 3
        error("length of ranges must be 3")
    end
    
    loss_sw_reg_1 = get_loss(medians.sw_reg, 1, actual)
    loss_sw_dfm_1 = get_loss(medians.sw_dfm, 1, actual)
    loss_sw_mes_1 = get_loss(medians.sw_mes, 1, actual)
    loss_swff_reg_1 = get_loss(medians.swff_reg, 1, actual)
    loss_swff_dfm_1 = get_loss(medians.swff_dfm, 1, actual)
    loss_swff_mes_1 = get_loss(medians.swff_mes, 1, actual)

    loss_sw_reg_2 = get_loss(medians.sw_reg, 2, actual)
    loss_sw_dfm_2 = get_loss(medians.sw_dfm, 2, actual)
    loss_sw_mes_2 = get_loss(medians.sw_mes, 2, actual)
    loss_swff_reg_2 = get_loss(medians.swff_reg, 2, actual)
    loss_swff_dfm_2 = get_loss(medians.swff_dfm, 2, actual)
    loss_swff_mes_2 = get_loss(medians.swff_mes, 2, actual)

    loss_sw_reg_4 = get_loss(medians.sw_reg, 4, actual)
    loss_sw_dfm_4 = get_loss(medians.sw_dfm, 4, actual)
    loss_sw_mes_4 = get_loss(medians.sw_mes, 4, actual)
    loss_swff_reg_4 = get_loss(medians.swff_reg, 4, actual)
    loss_swff_dfm_4 = get_loss(medians.swff_dfm, 4, actual)
    loss_swff_mes_4 = get_loss(medians.swff_mes, 4, actual)

        
    # SW-Reg vs SWFF-Reg
    s1 = "SW-Reg vs SWFF-Reg & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_swff_reg_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_swff_reg_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_swff_reg_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_swff_reg_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_swff_reg_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_swff_reg_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_swff_reg_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_swff_reg_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_swff_reg_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_swff_reg_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_swff_reg_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_swff_reg_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_swff_reg_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_swff_reg_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_swff_reg_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_swff_reg_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_swff_reg_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_swff_reg_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SW-Rich vs SWFF-Rich
    s2="SW-Rich vs SWFF-Rich & " *
    "$(round(dm_test(values(loss_sw_dfm_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SW-Reg vs SW-Rich
    s3="SW-Reg vs SW-Rich & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SW-Reg vs SW-Sparse
    s4="SW-Reg vs SW-Sparse & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_mes_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_mes_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_mes_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_mes_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_mes_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_mes_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_mes_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_mes_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_mes_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_mes_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_mes_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_mes_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_mes_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_mes_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_mes_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_mes_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_mes_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_mes_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"


    # SW-Sparse vs SW-Rich
    s5="SW-Sparse vs SW-Rich & " *
    "$(round(dm_test(values(loss_sw_mes_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_mes_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_mes_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_mes_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_mes_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_mes_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_mes_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_mes_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_mes_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_mes_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"


    # SWFF-Reg vs SWFF-Rich
    s6="SWFF-Reg vs SWFF-Rich & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"



    # SWFF-Reg vs SWFF-Sparse
    s7="SWFF-Reg vs SWFF-Sparse & "*
    "$(round(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_mes_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_mes_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_mes_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_mes_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_mes_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_mes_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_mes_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_mes_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_mes_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_mes_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_mes_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_mes_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_mes_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_mes_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_mes_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_mes_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_mes_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_mes_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SWFF-Sparse vs SWFF-Rich
    s8="SWFF-Sparse vs SWFF-Rich & "*
    "$(round(dm_test(values(loss_swff_mes_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_mes_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_mes_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_mes_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_mes_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_mes_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_mes_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_mes_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_mes_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_mes_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    println("""
    $s1  \\\\
    $s2 \\\\
    \\cmidrule(lr){1-10}
    $s3 \\\\
    $s4 \\\\
    $s5 \\\\
    \\cmidrule(lr){1-10}
    $s6 \\\\
    $s7 \\\\
    $s8 \\\\
    """)
end

"""
    print_cumulative_dm_table(medians::Workspace, actual::TSeries, ranges)

Print a LaTeX-formatted table of Diebold-Mariano test statistics comparing different models and forecast horizons.

# Arguments
- `medians::Workspace`: Workspace containing median forecasts from different models (e.g., `medians.sw_reg`, `medians.sw_dfm`).
- `actual::TSeries`: TSeries of actual values.
- `ranges`: A vector of three date ranges (MIT) over which to calculate the DM statistics.
"""
function print_cumulative_dm_table(medians::Workspace, actual::TSeries, ranges)
    if length(ranges) !== 3
        error("length of ranges must be 3")
    end
    loss_sw_reg_1 = get_cumulative_loss(medians.sw_reg, 1, actual)
    loss_sw_dfm_1 = get_cumulative_loss(medians.sw_dfm, 1, actual)
    loss_sw_mes_1 = get_cumulative_loss(medians.sw_mes, 1, actual)
    loss_swff_reg_1 = get_cumulative_loss(medians.swff_reg, 1, actual)
    loss_swff_dfm_1 = get_cumulative_loss(medians.swff_dfm, 1, actual)
    loss_swff_mes_1 = get_cumulative_loss(medians.swff_mes, 1, actual)

    loss_sw_reg_2 = get_cumulative_loss(medians.sw_reg, 2, actual)
    loss_sw_dfm_2 = get_cumulative_loss(medians.sw_dfm, 2, actual)
    loss_sw_mes_2 = get_cumulative_loss(medians.sw_mes, 2, actual)
    loss_swff_reg_2 = get_cumulative_loss(medians.swff_reg, 2, actual)
    loss_swff_dfm_2 = get_cumulative_loss(medians.swff_dfm, 2, actual)
    loss_swff_mes_2 = get_cumulative_loss(medians.swff_mes, 2, actual)

    loss_sw_reg_4 = get_cumulative_loss(medians.sw_reg, 4, actual)
    loss_sw_dfm_4 = get_cumulative_loss(medians.sw_dfm, 4, actual)
    loss_sw_mes_4 = get_cumulative_loss(medians.sw_mes, 4, actual)
    loss_swff_reg_4 = get_cumulative_loss(medians.swff_reg, 4, actual)
    loss_swff_dfm_4 = get_cumulative_loss(medians.swff_dfm, 4, actual)
    loss_swff_mes_4 = get_cumulative_loss(medians.swff_mes, 4, actual)

        
    # SW-Reg vs SWFF-Reg
    s1 = "SW-Reg vs SWFF-Reg & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_swff_reg_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_swff_reg_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_swff_reg_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_swff_reg_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_swff_reg_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_swff_reg_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_swff_reg_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_swff_reg_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_swff_reg_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_swff_reg_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_swff_reg_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_swff_reg_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_swff_reg_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_swff_reg_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_swff_reg_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_swff_reg_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_swff_reg_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_swff_reg_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SW-Rich vs SWFF-Rich
    s2="SW-Rich vs SWFF-Rich & " *
    "$(round(dm_test(values(loss_sw_dfm_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SW-Reg vs SW-Rich
    s3="SW-Reg vs SW-Rich & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SW-Reg vs SW-Sparse
    s4="SW-Reg vs SW-Sparse & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_mes_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[1]]-loss_sw_mes_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_mes_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[1]]-loss_sw_mes_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_mes_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[1]]-loss_sw_mes_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_mes_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[2]]-loss_sw_mes_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_mes_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[2]]-loss_sw_mes_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_mes_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[2]]-loss_sw_mes_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_mes_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_reg_1[ranges[3]]-loss_sw_mes_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_mes_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_reg_2[ranges[3]]-loss_sw_mes_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_mes_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_reg_4[ranges[3]]-loss_sw_mes_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"


    # SW-Sparse vs SW-Rich
    s5="SW-Sparse vs SW-Rich & " *
    "$(round(dm_test(values(loss_sw_mes_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_mes_1[ranges[1]]-loss_sw_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_mes_2[ranges[1]]-loss_sw_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_mes_4[ranges[1]]-loss_sw_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_mes_1[ranges[2]]-loss_sw_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_mes_2[ranges[2]]-loss_sw_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_mes_4[ranges[2]]-loss_sw_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_mes_1[ranges[3]]-loss_sw_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_mes_2[ranges[3]]-loss_sw_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_mes_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_mes_4[ranges[3]]-loss_sw_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"


    # SWFF-Reg vs SWFF-Rich
    s6="SWFF-Reg vs SWFF-Rich & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"



    # SWFF-Reg vs SWFF-Sparse
    s7="SWFF-Reg vs SWFF-Sparse & "*
    "$(round(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_mes_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[1]]-loss_swff_mes_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_mes_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[1]]-loss_swff_mes_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_mes_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[1]]-loss_swff_mes_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_mes_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[2]]-loss_swff_mes_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_mes_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[2]]-loss_swff_mes_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_mes_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[2]]-loss_swff_mes_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_mes_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_reg_1[ranges[3]]-loss_swff_mes_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_mes_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_reg_2[ranges[3]]-loss_swff_mes_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_mes_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_reg_4[ranges[3]]-loss_swff_mes_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    # SWFF-Sparse vs SWFF-Rich
    s8="SWFF-Sparse vs SWFF-Rich & "*
    "$(round(dm_test(values(loss_swff_mes_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_mes_1[ranges[1]]-loss_swff_dfm_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_mes_2[ranges[1]]-loss_swff_dfm_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_mes_4[ranges[1]]-loss_swff_dfm_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_mes_1[ranges[2]]-loss_swff_dfm_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_mes_2[ranges[2]]-loss_swff_dfm_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_mes_4[ranges[2]]-loss_swff_dfm_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_mes_1[ranges[3]]-loss_swff_dfm_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_mes_2[ranges[3]]-loss_swff_dfm_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_mes_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_mes_4[ranges[3]]-loss_swff_dfm_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    println("""
    $s1  \\\\
    $s2 \\\\
    \\cmidrule(lr){1-10}
    $s3 \\\\
    $s4 \\\\
    $s5 \\\\
    \\cmidrule(lr){1-10}
    $s6 \\\\
    $s7 \\\\
    $s8 \\\\
    """)
end

function print_cumulative_dm_table_full_vs_part(medians_full::Workspace, medians_part::Workspace, actual::TSeries, ranges)
    if length(ranges) !== 3
        error("length of ranges must be 3")
    end
    loss_sw_reg_full_1 = get_cumulative_loss(medians_full.sw_reg, 1, actual)
    loss_sw_dfm_full_1 = get_cumulative_loss(medians_full.sw_dfm, 1, actual)
    loss_sw_mes_full_1 = get_cumulative_loss(medians_full.sw_mes, 1, actual)
    loss_swff_reg_full_1 = get_cumulative_loss(medians_full.swff_reg, 1, actual)
    loss_swff_dfm_full_1 = get_cumulative_loss(medians_full.swff_dfm, 1, actual)
    loss_swff_mes_full_1 = get_cumulative_loss(medians_full.swff_mes, 1, actual)

    loss_sw_reg_full_2 = get_cumulative_loss(medians_full.sw_reg, 2, actual)
    loss_sw_dfm_full_2 = get_cumulative_loss(medians_full.sw_dfm, 2, actual)
    loss_sw_mes_full_2 = get_cumulative_loss(medians_full.sw_mes, 2, actual)
    loss_swff_reg_full_2 = get_cumulative_loss(medians_full.swff_reg, 2, actual)
    loss_swff_dfm_full_2 = get_cumulative_loss(medians_full.swff_dfm, 2, actual)
    loss_swff_mes_full_2 = get_cumulative_loss(medians_full.swff_mes, 2, actual)

    loss_sw_reg_full_4 = get_cumulative_loss(medians_full.sw_reg, 4, actual)
    loss_sw_dfm_full_4 = get_cumulative_loss(medians_full.sw_dfm, 4, actual)
    loss_sw_mes_full_4 = get_cumulative_loss(medians_full.sw_mes, 4, actual)
    loss_swff_reg_full_4 = get_cumulative_loss(medians_full.swff_reg, 4, actual)
    loss_swff_dfm_full_4 = get_cumulative_loss(medians_full.swff_dfm, 4, actual)
    loss_swff_mes_full_4 = get_cumulative_loss(medians_full.swff_mes, 4, actual)


    loss_sw_reg_part_1 = get_cumulative_loss(medians_part.sw_reg, 1, actual)
    loss_sw_dfm_part_1 = get_cumulative_loss(medians_part.sw_dfm, 1, actual)
    loss_sw_mes_part_1 = get_cumulative_loss(medians_part.sw_mes, 1, actual)
    loss_swff_reg_part_1 = get_cumulative_loss(medians_part.swff_reg, 1, actual)
    loss_swff_dfm_part_1 = get_cumulative_loss(medians_part.swff_dfm, 1, actual)
    loss_swff_mes_part_1 = get_cumulative_loss(medians_part.swff_mes, 1, actual)

    loss_sw_reg_part_2 = get_cumulative_loss(medians_part.sw_reg, 2, actual)
    loss_sw_dfm_part_2 = get_cumulative_loss(medians_part.sw_dfm, 2, actual)
    loss_sw_mes_part_2 = get_cumulative_loss(medians_part.sw_mes, 2, actual)
    loss_swff_reg_part_2 = get_cumulative_loss(medians_part.swff_reg, 2, actual)
    loss_swff_dfm_part_2 = get_cumulative_loss(medians_part.swff_dfm, 2, actual)
    loss_swff_mes_part_2 = get_cumulative_loss(medians_part.swff_mes, 2, actual)

    loss_sw_reg_part_4 = get_cumulative_loss(medians_part.sw_reg, 4, actual)
    loss_sw_dfm_part_4 = get_cumulative_loss(medians_part.sw_dfm, 4, actual)
    loss_sw_mes_part_4 = get_cumulative_loss(medians_part.sw_mes, 4, actual)
    loss_swff_reg_part_4 = get_cumulative_loss(medians_part.swff_reg, 4, actual)
    loss_swff_dfm_part_4 = get_cumulative_loss(medians_part.swff_dfm, 4, actual)
    loss_swff_mes_part_4 = get_cumulative_loss(medians_part.swff_mes, 4, actual)

        
    # SW-Reg vs SWFF-Reg
    s1 = "SW-Rich (full vs fewer series) & " *
    "$(round(dm_test(values(loss_sw_dfm_full_1[ranges[1]]-loss_sw_dfm_part_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_1[ranges[1]]-loss_sw_dfm_part_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_2[ranges[1]]-loss_sw_dfm_part_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_2[ranges[1]]-loss_sw_dfm_part_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_4[ranges[1]]-loss_sw_dfm_part_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_4[ranges[1]]-loss_sw_dfm_part_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_1[ranges[2]]-loss_sw_dfm_part_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_1[ranges[2]]-loss_sw_dfm_part_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_2[ranges[2]]-loss_sw_dfm_part_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_2[ranges[2]]-loss_sw_dfm_part_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_4[ranges[2]]-loss_sw_dfm_part_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_4[ranges[2]]-loss_sw_dfm_part_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_1[ranges[3]]-loss_sw_dfm_part_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_1[ranges[3]]-loss_sw_dfm_part_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_2[ranges[3]]-loss_sw_dfm_part_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_2[ranges[3]]-loss_sw_dfm_part_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_sw_dfm_full_4[ranges[3]]-loss_sw_dfm_part_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_sw_dfm_full_4[ranges[3]]-loss_sw_dfm_part_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    s2 = "SWFF-Rich (full vs fewer series) & " *
    "$(round(dm_test(values(loss_swff_dfm_full_1[ranges[1]]-loss_swff_dfm_part_1[ranges[1]]), 1)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_1[ranges[1]]-loss_swff_dfm_part_1[ranges[1]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_2[ranges[1]]-loss_swff_dfm_part_2[ranges[1]]), 2)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_2[ranges[1]]-loss_swff_dfm_part_2[ranges[1]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_4[ranges[1]]-loss_swff_dfm_part_4[ranges[1]]), 4)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_4[ranges[1]]-loss_swff_dfm_part_4[ranges[1]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_1[ranges[2]]-loss_swff_dfm_part_1[ranges[2]]), 1)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_1[ranges[2]]-loss_swff_dfm_part_1[ranges[2]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_2[ranges[2]]-loss_swff_dfm_part_2[ranges[2]]), 2)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_2[ranges[2]]-loss_swff_dfm_part_2[ranges[2]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_4[ranges[2]]-loss_swff_dfm_part_4[ranges[2]]), 4)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_4[ranges[2]]-loss_swff_dfm_part_4[ranges[2]]), 4)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_1[ranges[3]]-loss_swff_dfm_part_1[ranges[3]]), 1)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_1[ranges[3]]-loss_swff_dfm_part_1[ranges[3]]), 1)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_2[ranges[3]]-loss_swff_dfm_part_2[ranges[3]]), 2)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_2[ranges[3]]-loss_swff_dfm_part_2[ranges[3]]), 2)[2] < 0.05 ? "*" : "")" * " & " *
    "$(round(dm_test(values(loss_swff_dfm_full_4[ranges[3]]-loss_swff_dfm_part_4[ranges[3]]), 4)[1], digits=1))$(dm_test(values(loss_swff_dfm_full_4[ranges[3]]-loss_swff_dfm_part_4[ranges[3]]), 4)[2] < 0.05 ? "*" : "")"

    println("""
    $s1  \\\\
    $s2 \\\\
    \\cmidrule(lr){1-10}
    """)
end

function print_cumulative_dm_table_multiple_targets(medians::Vector{Workspace}, actuals::Vector{<:TSeries}, range)
    
    losses = Vector{Workspace}()
    for i = 1:length(medians)
        l = Workspace()
        medians_full = medians[i]
        actual = actuals[i]

        l.loss_sw_reg_1 = get_cumulative_loss(medians_full.sw_reg, 1, actual)
        l.loss_sw_dfm_1 = get_cumulative_loss(medians_full.sw_dfm, 1, actual)
        l.loss_sw_mes_1 = get_cumulative_loss(medians_full.sw_mes, 1, actual)
        l.loss_swff_reg_1 = get_cumulative_loss(medians_full.swff_reg, 1, actual)
        l.loss_swff_dfm_1 = get_cumulative_loss(medians_full.swff_dfm, 1, actual)
        l.loss_swff_mes_1 = get_cumulative_loss(medians_full.swff_mes, 1, actual)

        l.loss_sw_reg_2 = get_cumulative_loss(medians_full.sw_reg, 2, actual)
        l.loss_sw_dfm_2 = get_cumulative_loss(medians_full.sw_dfm, 2, actual)
        l.loss_sw_mes_2 = get_cumulative_loss(medians_full.sw_mes, 2, actual)
        l.loss_swff_reg_2 = get_cumulative_loss(medians_full.swff_reg, 2, actual)
        l.loss_swff_dfm_2 = get_cumulative_loss(medians_full.swff_dfm, 2, actual)
        l.loss_swff_mes_2 = get_cumulative_loss(medians_full.swff_mes, 2, actual)

        l.loss_sw_reg_4 = get_cumulative_loss(medians_full.sw_reg, 4, actual)
        l.loss_sw_dfm_4 = get_cumulative_loss(medians_full.sw_dfm, 4, actual)
        l.loss_sw_mes_4 = get_cumulative_loss(medians_full.sw_mes, 4, actual)
        l.loss_swff_reg_4 = get_cumulative_loss(medians_full.swff_reg, 4, actual)
        l.loss_swff_dfm_4 = get_cumulative_loss(medians_full.swff_dfm, 4, actual)
        l.loss_swff_mes_4 = get_cumulative_loss(medians_full.swff_mes, 4, actual)

        push!(losses, l)
    end
    
    # SW-Reg vs SWFF-Reg
    s1 = "SW-Reg vs SWFF-Reg" 
    for i = 1:length(medians)
        s1 = s1 * " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_1[range]-losses[i].loss_swff_reg_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_1[range]-losses[i].loss_swff_reg_1[range]), 1)[2] < 0.05 ? "*" : "")"  *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_2[range]-losses[i].loss_swff_reg_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_2[range]-losses[i].loss_swff_reg_2[range]), 2)[2] < 0.05 ? "*" : "")"  *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_4[range]-losses[i].loss_swff_reg_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_4[range]-losses[i].loss_swff_reg_4[range]), 4)[2] < 0.05 ? "*" : "")" 
    end

    # SW-Rich vs SWFF-Rich
    s2 = "SW-Rich vs SWFF-Rich" 
     for i = 1:length(medians)
        s2 = s2 * " & " *  "$(round(dm_test(values(losses[i].loss_sw_dfm_1[range]-losses[i].loss_swff_dfm_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_sw_dfm_1[range]-losses[i].loss_swff_dfm_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_dfm_2[range]-losses[i].loss_swff_dfm_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_sw_dfm_2[range]-losses[i].loss_swff_dfm_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_dfm_4[range]-losses[i].loss_swff_dfm_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_sw_dfm_4[range]-losses[i].loss_swff_dfm_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end

   # SW-Reg vs SW-Rich
    s3="SW-Reg vs SW-Rich" 
     for i = 1:length(medians)
        s3 = s3 * " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_1[range]-losses[i].loss_sw_dfm_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_1[range]-losses[i].loss_sw_dfm_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_2[range]-losses[i].loss_sw_dfm_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_2[range]-losses[i].loss_sw_dfm_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_4[range]-losses[i].loss_sw_dfm_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_4[range]-losses[i].loss_sw_dfm_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end

    # SW-Reg vs SW-Sparse
    s4="SW-Reg vs SW-Sparse" 
    for i = 1:length(medians)
        s4 = s4 * " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_1[range]-losses[i].loss_sw_mes_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_1[range]-losses[i].loss_sw_mes_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_2[range]-losses[i].loss_sw_mes_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_2[range]-losses[i].loss_sw_mes_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_reg_4[range]-losses[i].loss_sw_mes_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_sw_reg_4[range]-losses[i].loss_sw_mes_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end

    # SW-Sparse vs SW-Rich
    s5="SW-Sparse vs SW-Rich"
    for i = 1:length(medians)
        s5 = s5 * " & " *  "$(round(dm_test(values(losses[i].loss_sw_mes_1[range]-losses[i].loss_sw_dfm_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_sw_mes_1[range]-losses[i].loss_sw_dfm_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_mes_2[range]-losses[i].loss_sw_dfm_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_sw_mes_2[range]-losses[i].loss_sw_dfm_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_sw_mes_4[range]-losses[i].loss_sw_dfm_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_sw_mes_4[range]-losses[i].loss_sw_dfm_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end
    

    # SWFF-Reg vs SWFF-Rich
    s6="SWFF-Reg vs SWFF-Rich"
    for i = 1:length(medians)
        s6 = s6 * " & " *  "$(round(dm_test(values(losses[i].loss_swff_reg_1[range]-losses[i].loss_swff_dfm_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_swff_reg_1[range]-losses[i].loss_swff_dfm_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_swff_reg_2[range]-losses[i].loss_swff_dfm_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_swff_reg_2[range]-losses[i].loss_swff_dfm_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_swff_reg_4[range]-losses[i].loss_swff_dfm_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_swff_reg_4[range]-losses[i].loss_swff_dfm_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end

    # SWFF-Reg vs SWFF-Sparse
    s7="SWFF-Reg vs SWFF-Sparse"
    for i = 1:length(medians)
        s7 = s7 * " & " *  "$(round(dm_test(values(losses[i].loss_swff_reg_1[range]-losses[i].loss_swff_mes_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_swff_reg_1[range]-losses[i].loss_swff_mes_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_swff_reg_2[range]-losses[i].loss_swff_mes_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_swff_reg_2[range]-losses[i].loss_swff_mes_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_swff_reg_4[range]-losses[i].loss_swff_mes_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_swff_reg_4[range]-losses[i].loss_swff_mes_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end

    # SWFF-Sparse vs SWFF-Rich
    s8="SWFF-Sparse vs SWFF-Rich"
    for i = 1:length(medians)
        s8 = s8 * " & " *  "$(round(dm_test(values(losses[i].loss_swff_mes_1[range]-losses[i].loss_swff_dfm_1[range]), 1)[1], digits=1))$(dm_test(values(losses[i].loss_swff_mes_1[range]-losses[i].loss_swff_dfm_1[range]), 1)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_swff_mes_2[range]-losses[i].loss_swff_dfm_2[range]), 2)[1], digits=1))$(dm_test(values(losses[i].loss_swff_mes_2[range]-losses[i].loss_swff_dfm_2[range]), 2)[2] < 0.05 ? "*" : "")" *
                  " & " *  "$(round(dm_test(values(losses[i].loss_swff_mes_4[range]-losses[i].loss_swff_dfm_4[range]), 4)[1], digits=1))$(dm_test(values(losses[i].loss_swff_mes_4[range]-losses[i].loss_swff_dfm_4[range]), 4)[2] < 0.05 ? "*" : "")"
    end


    println("""
    $s1  \\\\
    $s2 \\\\
    \\cmidrule(lr){1-7}
    $s3 \\\\
    $s4 \\\\
    $s5 \\\\
    \\cmidrule(lr){1-7}
    $s6 \\\\
    $s7 \\\\
    $s8 \\\\
    """)

end
 
"""
    get_most_median_draw(chain; burnin=75000)

Find the parameter draw that is closest to the median of all parameters after a burn-in period.

This function calculates the median value for each parameter across all post-burn-in draws.
Then, for each individual draw, it computes the sum of squared differences between its
parameter values and the calculated medians. The draw with the smallest sum of squared
differences is considered the "most median" draw.

# Arguments
- `chain`: A Workspace object containing the MCMC chain results, specifically `chain.𝛉` which is a matrix of parameter draws (draws x parameters).
- `burnin=75000`: The number of initial draws to discard as burn-in.

# Returns
- `lowest_index::Int`: The index (row number) in `chain.𝛉` corresponding to the draw most similar to the median parameter values.
"""
function get_most_median_draw(chain; burnin=75000)
    medians = StatsBase.median(chain.𝛉[burnin+1:end, :], dims=1)
    lowest_err = Inf
    lowest_index = burnin+1
    for i = burnin+1:size(chain.𝛉, 1)
        err = sum((chain.𝛉[i, :] .- medians) .^ 2)
        if err < lowest_err
            lowest_err = err
            lowest_index = i
        end
    end
    return lowest_index
end

"""
    get_state_errors(m, θ, S, time_index)

Calculate the structural shock values (state errors) for a given time period.

This function backs out the structural shocks `ϵ` from the state transition equation:
`Sₜ = G Sₜ₋₁ + H ϵₜ`.
It assumes that shocks can have AR(1) processes, i.e., `εˢʰᵏ_t = ρˢʰᵏ * εˢʰᵏ_{t-1} + σˢʰᵏ * ϵˢʰᵏ_t`,
where `εˢʰᵏ` is the shock level (a state variable) and `ϵˢʰᵏ` is the i.i.d. structural shock.
If no `ρ` parameter is found for a shock, it's assumed `ρ = 0`.

# Arguments
- `m`: The model object (`MacroModelling.ℳ`), containing parameter and shock definitions.
- `θ`: A vector of parameter values.
- `S`: A KeyedArray representing the time series of states (variables x time).
- `time_index`: The specific time index for which to calculate the shocks.

# Returns
- `shock_vals::Vector{Float64}`: A vector of calculated structural shock values for the given `time_index`, corresponding to `m.exo`.
"""
function get_state_errors(m, θ, S, time_index)
    current_state = S(time=time_index)
    prev_state = S(time=time_index-1)
    shock_vals = zeros(length(m.exo))
    for (i, shk) ∈ enumerate(m.exo)
        shk_rho_sym = Symbol("ρ$(replace(string(shk), "ϵ"=>""))")
        shk_rho_idx = findfirst(s -> s == shk_rho_sym, m.parameters)
        shk_sigma_sym = Symbol("σ$(replace(string(shk), "ϵ"=>""))")
        shk_sigma_idx = findfirst(s -> s == shk_sigma_sym, m.parameters)
        shk_lvl_sym =  Symbol(replace(string(shk), "ϵ"=>"ε"))
        if shk_rho_idx !== nothing
            shock_vals[i] = (current_state(Variable=shk_lvl_sym) - θ[shk_rho_idx]*prev_state(Variable=shk_lvl_sym)) / θ[shk_sigma_idx]
        else
            shock_vals[i] = current_state(Variable=shk_lvl_sym) / θ[shk_sigma_idx]
        end
    end
    return shock_vals
end

"""
    generate_shock_decomp(m::MacroModelling.ℳ, df::DataFrame, chain::Workspace, startdate::MIT, enddate::MIT; speedup=true, data_X_start=1986Q2, core=false, excluded_series=[], burnin=75000)

Generate a historical shock decomposition using the median parameter draw from an MCMC chain.

This function first identifies the parameter draw closest to the median of the posterior distribution
(after burn-in) from the `chain` Workspace. It then calls the main `generate_shock_decomp`
method with these median parameters to compute the decomposition.

# Arguments
- `m::MacroModelling.ℳ`: The model object.
- `df::DataFrame`: DataFrame containing the historical data.
- `chain::Workspace`: Workspace containing the MCMC estimation results (parameters `𝛉`, `R`, `𝛙`, `𝚲`).
- `startdate::MIT`: The start date for the shock decomposition analysis.
- `enddate::MIT`: The end date for the shock decomposition analysis.
- `speedup=true`: Boolean indicating whether to use speedup techniques (e.g., for DFM).
- `data_X_start=1986Q2`: The start date of the data used for estimation.
- `core=false`: Boolean indicating if only core series were used for estimation.
- `excluded_series=[]`: A list of series to exclude from the data.
- `burnin=75000`: Number of initial MCMC draws to discard as burn-in when finding the median draw.

# Returns
- `decomp::Workspace`: A Workspace containing the shock decomposition results. For each shock in `m.exo`, it stores a sub-Workspace with `:states` (contribution to state variables) and `:data` (contribution to observed data).
"""
function generate_shock_decomp(m::MacroModelling.ℳ, df::DataFrame, chain::Workspace, startdate::MIT, enddate::MIT; 
    speedup=true, data_X_start=1986Q2, core=false, excluded_series=[], burnin=75000)

    median_draw = get_most_median_draw(chain, burnin=burnin)
    θ = chain.𝛉[median_draw,:]
    R = Diagonal(chain.R[median_draw,:])
    𝛙 = Diagonal(chain.𝛙[median_draw,:])
    𝚲 = chain.𝚲[median_draw]

    decomp = generate_shock_decomp(m, df, startdate, enddate, θ, R, 𝛙, 𝚲; 
        speedup=speedup, data_X_start=data_X_start, core=core, excluded_series=excluded_series)
    return decomp

end

"""
    generate_shock_decomp(m::MacroModelling.ℳ, df::DataFrame, startdate::MIT, enddate::MIT, θ::Vector, R::Diagonal, 𝛙::Diagonal, 𝚲::Matrix; speedup=false, data_X_start=1986Q2, core=false, excluded_series=[])

Generate a historical shock decomposition for a given set of parameters and data.

This function computes the contribution of each structural shock to the historical evolution
of state variables and observed data. It does this by:
1. Filtering states using the Kalman filter (`generate_states`).
2. For each shock:
    a. Simulate the model forward from `startdate` to `enddate`, setting that specific shock to zero for all periods.
    b. The difference between the filtered historical series and this counterfactual simulation represents the contribution of that shock.

# Arguments
- `m::MacroModelling.ℳ`: The model object.
- `df::DataFrame`: DataFrame containing the historical data.
- `startdate::MIT`: The start date for the shock decomposition analysis.
- `enddate::MIT`: The end date for the shock decomposition analysis.
- `θ::Vector`: Vector of model parameter values.
- `R::Diagonal`: Measurement error covariance matrix.
- `𝛙::Diagonal`: Autoregressive coefficient matrix for measurement errors.
- `𝚲::Matrix`: Observation matrix linking states to observables.
- `speedup=false`: Boolean indicating whether to use speedup techniques (e.g., for DFM).
- `data_X_start=1986Q2`: The start date of the data used for estimation/filtering.
- `core=false`: Boolean indicating if only core series are used in `samp_data`.
- `excluded_series=[]`: A list of series to exclude from the data.

# Returns
- `forecasts_without_shock::Workspace`: A Workspace where keys are shock symbols (from `m.exo`). Each key maps to another Workspace containing:
    - `:states`: A KeyedArray (variables x time) of the contribution of that shock to each state variable.
    - `:data`: A DataFrame (time x variables) of the contribution of that shock to each observed data series.
"""
function generate_shock_decomp(m::MacroModelling.ℳ, df::DataFrame, startdate::MIT, enddate::MIT, θ::Vector, R::Diagonal, 𝛙::Diagonal, 𝚲::Matrix; speedup=false, data_X_start=1986Q2, core=false, excluded_series=[])
    TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    if !solved
        error("model doesn't solve")
    end
    G_keyed, H_keyed, SS = get_sparse_solution(m, 𝐒)
    G = Matrix(G_keyed)
    H = Matrix(H_keyed)
    
    # constrain data
    df_constrained = df[1:length(data_X_start:enddate),:]
    
    samp_data = get_data_obj(m, df_constrained, core, excluded_series)

    S, _ = generate_states(m, θ, 𝐒, samp_data, 𝚲, 𝛙, R; backwards_pass=false, speedup=speedup)
    # forecast
    

    measurement_error_levels = Matrix(samp_data.X[2:end, :])' - 𝚲* S
    measurement_errors = measurement_error_levels[:,2:end] - 𝛙*Matrix(measurement_error_levels[:,1:end-1])
    

    forecasts_without_shock = Workspace()
    forecast_start = length(startdate-data_X_start)+3

    for shk ∈ m.exo
        shk_idx = findfirst(s -> s == shk, m.exo)
        states_forecast = copy(S(:, 1:length(data_X_start:enddate)))
        states_forecast .= 0
        states_forecast[:,1:forecast_start-1] = copy(S[:,1:forecast_start-1])

        data_forecast = copy(samp_data.X[2:end,:])
        data_forecast .= 0
        data_forecast[1:forecast_start-1,:] = copy(samp_data.X[2:forecast_start, :])

        measure_lag = Vector(samp_data.X[3,:]) - 𝚲* S[:,2] 

        for i = forecast_start:size(S, 2)
            pred_states = rekey(Matrix(G) * states_forecast(time=i-1), 1 => G_keyed.Variables)
            v = get_state_errors(m, θ, S, i)
            # v = get_state_errors3(S, G, H, i)
            # v = get_state_errors2(pred_states_perfect, S, H, i)
            v[shk_idx] = 0
            states_forecast[:,i] = pred_states + H * v
            # println(states_forecast[:,i])
            # println(S[:,i])
            # println(v)
            # println(v2)
            # @assert states_forecast[:,i] ≈ S[:,i]

            measure = 𝛙*measure_lag + measurement_errors(time=i);
            data_forecast[i,:] = 𝚲 * states_forecast[:,i] + measure;

            # next period's lagged measurement error is this period's measurment error
            measure_lag=measure;

        end

        forecasts_without_shock[shk] = Workspace(
            :states => copy(S .- states_forecast),
            :data => copy(samp_data.X[2:end, :] .- data_forecast),
        )

    end

    return forecasts_without_shock

end

"""
    get_shock_decomp_mvts(components, startdate, enddate, actuals, shocks=[:ϵG, :ϵI, :ϵa, :ϵb, :ϵp, :ϵq, :ϵF, :ϵr, :ϵw])

Organize shock decomposition results into MVTSeries objects for plotting and analysis.

This function takes the raw shock decomposition contributions (differences from a no-shock counterfactual)
and transforms them into time series representing the cumulative impact of each shock, normalized
to start at zero at the specified `startdate`. It also calculates the contribution of
measurement error as the residual between the sum of identified shock contributions and the actual data.

# Arguments
- `components`: A Workspace generated by `generate_shock_decomp`, where keys are shock symbols and values are Workspaces containing `:data` (DataFrame of shock contributions to observables).
- `startdate`: The desired start date for the output MVTSeries. Contributions will be normalized to be zero at this date.
- `enddate`: The desired end date for the output MVTSeries.
- `actuals`: A Workspace or similar structure containing the TSeries of actual historical data for each observable. Keys should be variable symbols (e.g., `:RGDP`).
- `shocks=[:ϵG, ..., :ϵw]`: A vector of shock symbols to include in the decomposition. Defaults to a standard set.

# Returns
- `decomps::Workspace`: A Workspace where keys are observable variable symbols (e.g., `:RGDP`). Each key maps to an `MVTSeries`.
    The `MVTSeries` has time chainning from `startdate` to `enddate`. Its columns correspond to the specified `shocks` plus an additional `:m.e.` (measurement error) column.
    Values represent the cumulative contribution of each shock (or m.e.) to the deviation of the variable from its value at `startdate`.
"""
function get_shock_decomp_mvts(components, startdate, enddate, actuals, shocks=[:ϵG, :ϵI, :ϵa, :ϵb, :ϵp, :ϵq, :ϵF, :ϵr, :ϵw])
    first_shock = first(keys(components))
    vars = names(components[first_shock].data)
    rng = enddate-size(components[first_shock].data,1)+1:enddate
    # rng = enddate-size(components[first_shock].data,1)+2:enddate
    @assert length(rng) == size(components[first_shock].data,1)
    decomps = Workspace()
    mvts_cols = filter(x -> x ∈ keys(components), shocks)
    push!(mvts_cols, Symbol("m.e."))
    for var ∈ vars
        mvts = MVTSeries(rng, mvts_cols, zeros)
        for shk ∈ shocks
            if shk ∉ keys(components)
                continue
            end
            mvts[rng, shk] = components[shk].data[!, var]
            mvts[shk] = mvts[shk] .- mvts[shk][startdate]
            # mvts[shk] = mvts[shk] ./ mvts[shk][startdate]
            # mvts[rng, shk] = components[shk].states(Variable=:Y)
        end
        total = TSeries(rng, sum(mvts, dims=2))
        actual = actuals[Symbol(var)] .- actuals[Symbol(var)][startdate]
        # actual = actuals[Symbol(var)] ./ actuals[Symbol(var)][startdate]
        mvts[Symbol("m.e.")] = actual - total

        mvts = mvts[startdate:enddate]
        # println(mvts[startdate])
        # mvst.values = mvts.values .- mvts[startdate]
        decomps[Symbol(var)] = deepcopy(mvts)
    end
    return decomps
end

"""
    plot_mvts(mvts::MVTSeries; bar_width=0.7, title="", kwargs...)

Create a stacked bar plot from an MVTSeries, typically used for visualizing shock decompositions.

This function generates a grouped bar chart where each group corresponds to a time period
in the `mvts` and bars within each group are stacked, representing the values of different
series (e.g., shock contributions) at that time. A line plot overlaying the bars shows the
sum of the series at each time point.

# Arguments
- `mvts::MVTSeries`: The multivariate time series data to plot. Rows are time periods, columns are different series.
- `bar_width=0.7`: The width of the bars in the plot.
- `title=""`: The title for the plot.
- `kwargs...`: Additional keyword arguments to be passed to `StatsPlots.groupedbar` and `Plots.plot`.

# Returns
- `p::Plots.Plot`: A plot object representing the stacked bar chart.

# Example
```julia
# Assuming `decomp_mvts` is an MVTSeries from a shock decomposition
# plotting_variable = :Rgdp
# p = plot_mvts(decomps[plotting_variable], title="Shock Decomposition for \$plotting_variable")
# display(p)
```
"""
function plot_mvts(mvts; bar_width=0.7, title="", legend=:best, kwargs...)
    num_entries = length(keys(mvts))
    rng = rangeof(mvts)
    times = string.(collect(rng))
    color=reshape(collect(cgrad(:roma, num_entries, categorical = true)), (1,num_entries))
    println(kwargs)
    p = StatsPlots.groupedbar(times, mvts.values, 
        bar_position = :stack, 
        labels=reshape(string.(collect(keys(mvts))), (1, length(keys(mvts)))),
        bar_width=bar_width, 
        color=color, xrotation=60, title=title, titlefont=font(10,"Bookman Light"), legend=legend)
    StatsPlots.plot!(p, times, sum(mvts.values, dims=2), color="black", label=nothing, width=3)
    return p

end
#  palette(get_color_palette(:auto, plot_color(:white)))
#  palette(get_color_palette(:matter, plot_color(:white)))

#  palette(cgrad(:matter, 4, categorical = true))
#  palette(cgrad(:curl, 4, categorical = true))


#  cgrad(:matter, 5, categorical = true)

"""
    prior_posterior_charts(
        models::Vector{MacroModelling.ℳ},
        chains::Vector{Workspace};
        burnin=75000,
        param_indices=1:length(models[1].parameters),
        chain_titles=repeat([""], length(chains)),
        layout=:auto,
        size=(1000, 800)
    )
    prior_posterior_charts(model::MacroModelling.ℳ, chain::Workspace; kwargs...)
    prior_posterior_charts(model::MacroModelling.ℳ, chains::Vector{Workspace}; kwargs...)

Generate and plot density comparisons of prior and posterior distributions for model parameters.

This function creates a grid of plots, where each plot corresponds to a specific model parameter.
In each plot, it shows the density of the prior distribution (derived from the first model in `models`)
and the posterior densities obtained from one or more MCMC `chains`.

# Arguments
- `models::Vector{MacroModelling.ℳ}`: A vector of model objects. The prior distribution is taken from `models[1]`. If multiple chains are from different (but comparable) models, provide them accordingly.
- `chains::Vector{Workspace}`: A vector of Workspace objects, each containing MCMC chain results (specifically `chain.𝛉` for parameter draws).
- `burnin=75000`: The number of initial MCMC draws to discard from each chain for posterior density estimation.
- `param_indices=1:length(models[1].parameters)`: A collection of indices specifying which parameters (by their order in `models[1].parameters`) to plot. Defaults to all parameters.
- `chain_titles=repeat([""], length(chains))`: A vector of strings to append to the legend entries for posteriors from different chains, helping to distinguish them.
- `layout=:auto`: The layout of the subplots in the final plot. Defaults to an automatic arrangement. Can be an integer or a tuple (rows, cols).
- `size=(1000, 800)`: The overall size of the generated plot.

# Returns
- `priors_plot::Plots.Plot`: A plot object containing the grid of prior-posterior density comparisons.

# Notes
- The prior for each parameter `p` is sampled using an internal call like `rand(prior_θ(Val(Symbol(models[1].model_name)), Val(models[1].parameters[p])), 75000)`. Ensure `prior_θ` is defined and accessible.
- Convenience methods are provided for single model/single chain and single model/multiple chains scenarios.
"""
prior_posterior_charts(model::MacroModelling.ℳ, chain::Workspace; kwargs...) = prior_posterior_charts([model], [chain]; kwargs...)
prior_posterior_charts(model::MacroModelling.ℳ, chains::Vector{Workspace}; kwargs...) = prior_posterior_charts(repeat([model], length(chains)), chains; kwargs...)
function prior_posterior_charts(models::Vector{MacroModelling.ℳ}, chains::Vector{Workspace}; burnin=75000, param_indices=1:length(models[1].parameters), 
    chain_titles=repeat([""], length(chains)), 
    layout=:auto,
    size=(1000, 800))
    q_plots = Vector{Any}()
    for param_index ∈ param_indices
        for (i, chain) ∈ enumerate(chains)
            m = models[i]
            
            if i == 1
                p = StatsPlots.density(chain.𝛉[burnin+1:end, param_index], label="$(m.parameters[param_index])_posterior$(chain_titles[i])", title="$(m.parameters[param_index])", titlefont=font(10,"Bookman Light"))
            else
                p = StatsPlots.density!(p, chain.𝛉[burnin+1:end, param_index], label="$(m.parameters[param_index])_posterior$(chain_titles[i])", title="$(m.parameters[param_index])", titlefont=font(10,"Bookman Light"))
            end 
        end
        # add prior
        p = StatsPlots.density!(p, rand(prior_θ(Val(Symbol(models[1].model_name)), Val(models[1].parameters[param_index])), 75000), label="$(models[1].parameters[param_index])_prior")
        push!(q_plots, p)
    end

    _layout = layout
    if layout == :auto
       _layout = length(q_plots)
    end

    priors_plot = Plots.plot(q_plots..., layout=_layout, size=size, suptitle="Priors and Posteriors", plot_titlefontfamily="Bookman Light", bottom_margin=30Plots.px)
    return priors_plot
end

# General idea:
# 1. Backout all the shocks.
# 2. Run the simulation forward with all the shock, but setting 1 to zero.
# 3. Shock contribution in each quarter is Y - Y*
# 4. Do this for each shock...
function make_shock_decomp_charts(
    target::Symbol,
    m_sw::MacroModelling.ℳ,
    m_swff::MacroModelling.ℳ,
    bundle::Workspace,
    mit_start::MIT,
    mit_end::MIT;
    data_X_start=1986Q2,
    target_name = "_DEFAULT",
    suptitle="",
    legend=:best
    )

    if target_name == "_DEFAULT"
        target_name = string(target)
    end

    shock_decomp_components_SWReg = generate_shock_decomp(m_sw,  bundle.df_sw, bundle.chain_sw_reg, mit_start, mit_end; speedup=false, data_X_start=data_X_start, core=true, excluded_series=[])
    shock_decomp_components_SWRich = generate_shock_decomp(m_sw,  bundle.df_sw, bundle.chain_sw_dfm, mit_start, mit_end; speedup=true, data_X_start=data_X_start, core=false, excluded_series=[])
    shock_decomp_components_SWFFReg = generate_shock_decomp(m_swff,  bundle.df_swff, bundle.chain_swff_reg, mit_start, mit_end; speedup=false, data_X_start=data_X_start, core=true, excluded_series=[])
    shock_decomp_components_SWFFRich = generate_shock_decomp(m_swff,  bundle.df_swff, bundle.chain_swff_dfm, mit_start, mit_end; speedup=true, data_X_start=data_X_start, core=false, excluded_series=[])

    decomps_SWReg = get_shock_decomp_mvts(shock_decomp_components_SWReg, mit_start, mit_end, bundle.mvts_sw)
    decomps_SWRich = get_shock_decomp_mvts(shock_decomp_components_SWRich, mit_start, mit_end, bundle.mvts_sw)
    decomps_SWFFReg = get_shock_decomp_mvts(shock_decomp_components_SWFFReg, mit_start, mit_end, bundle.mvts_sw)
    decomps_SWFFRich = get_shock_decomp_mvts(shock_decomp_components_SWFFRich, mit_start, mit_end, bundle.mvts_sw)

    p_SWReg = plot_mvts(decomps_SWReg[target], title="SW-Reg", legend=legend)
    p_SWRich = plot_mvts(decomps_SWRich[target], title="SW-Rich", legend=legend)
    p_SWFFReg = plot_mvts(decomps_SWFFReg[target], title="SWFF-Reg", legend=legend)
    p_SWFFRich = plot_mvts(decomps_SWFFRich[target], title="SWFF-Rich", legend=legend)

    q_plots = [p_SWReg, p_SWRich, p_SWFFReg, p_SWFFRich]
    Plots.plot(q_plots..., layout=length(q_plots), size=(1000,800), suptitle=suptitle, plot_titlefontfamily="Bookman Light")
end



# Functions for loading Gelfer's results
gelfer_map_sw = Workspace(
    :Ψ => 1,
    :ιp => 2,
    :ιw => 3,
    :ξp => 4,
    :ξw => 5,
    :νl => 6,
    :σC => 7,
    :h => 8,
    :φ => 9,
    :Φ => 10,
    :rπ1 => 11,
    :ry1 => 12,
    :ρ => 13,
    :ρa => 14,
    :ρb => 15,
    :ρG => 16,
    :ρI => 17,
    :σa => 18,
    :σb => 19,
    :σG => 20,
    :σr => 21,
    :σI => 22,
    :σp => 23,
    :σw => 24,
    :σq => 25,
    :rπ2 => 26,
    :ry2 => 27,
    :ρp => 28,
    :ρw => 29,
)

gelfer_map_swff = Workspace(
    :Ψ => 1,
    :ιp => 2,
    :ιw => 3,
    :ξp => 4,
    :ξw => 5,
    :νl => 6,
    :σC => 7,
    :h => 8,
    :φ => 9,
    :Φ => 10,
    :rπ1 => 11,
    :ry1 => 12,
    :ρ => 13,
    :ρa => 14,
    :ρb => 15,
    :ρG => 16,
    :ρI => 17,
    :σa => 18,
    :σb => 19,
    :σG => 20,
    :σr => 21,
    :σI => 22,
    :σp => 23,
    :σw => 24,
    :σF => 25,
    :rπ2 => 26,
    :ry2 => 27,
    :ρF => 28,
    :χ => 29,
    :ρp => 30,
    :ρw => 31,
)

"""
    process_gelfer_mat(m, gelfer_mat, params_map, pad = 150000)

Structures the Gelfer output like our output.

# Arguments
- `m::MacroModelling.ℳ`: The model which produced the results.
- `gelfer_mat::Matrix{Float64}`: The parmeter posteriors from the loaded MAT file.
- `params_map::Workspace`: A Workspace mapping the model parameters to their indices in the MAT results.
- `pad::Integer`: The minimum length of the resultant Workspace.

# Returns
- `ret::Workspace`: A Workspace with a 𝛉 key containing the remapped parameter posteriors.
"""
function process_gelfer_mat(m, gelfer_mat, params_map, pad = 150000)
    # structure the Gelfer output like our output
    new_mat = zeros(size(gelfer_mat,1), length(m.parameter_values));
    for i ∈ 1:size(gelfer_mat, 2)
        param = m.parameters[i]
        gelfer_idx = params_map[param]
        new_mat[:, i] = gelfer_mat[:, gelfer_idx]
        if param == :χ
            # χ_aux = 0.0225+(0.0825*χ)
            new_mat[:, i] = 0.0225 .+ 0.0825 .* (gelfer_mat[:, gelfer_idx])
        end
    end
    z_mat = zeros(pad-size(new_mat, 1), size(new_mat,2))
    ret = Workspace(
        :𝛉 => vcat(z_mat, new_mat)
    )
    return ret
end

"""
    param_box_plots(m, chains, chain_names, gelfer_means; suptitle="Title", plot_size=(1000,9000), param_range=1:29,  burnin=75000)

Function for making figure 1.
"""
function param_box_plots(m, chains, chain_names, gelfer_means; suptitle="Title", plot_size=(1000,9000), param_range=1:29,  burnin=75000)
    df_combined = nothing
    for i ∈ 1:length(chains)
        df = DataFrame(chains[i].𝛉[burnin+1:end, :], m.parameters)
        df_long = stack(df)
        # df_long.chain .= "$(i)."
        df_long.chain .= "$(chain_names[i])."
        if i == 1
            df_combined = copy(df_long)
        else
            df_combined = vcat(df_combined, df_long)
        end
    end
    plots_vector = Vector{Any}()
    for param ∈ m.parameters[param_range]
        df_param = subset(df_combined, :variable => ByRow(x -> x == string(param))) 
        p_param = StatsPlots.@df df_param StatsPlots.violin(String.(:chain), :value, title=string(param), color="grey", legend=false, guidefontsize=4, tickfontsize=8, titlefont=font(10,"Bookman Light"))
        StatsPlots.@df df_param StatsPlots.boxplot!(p_param, String.(:chain), :value, fillalpha=0.60, color="white", outliers=false, guidefontsize=4, tickfontsize=8)
        for i ∈ 1:length(gelfer_means)
            if gelfer_means[i] isa Dict
                Plots.plot!(p_param, [i-1;i], [gelfer_means[i][param][1]; gelfer_means[i][param][1]], linewidth = 3, linestyle=:dash, color="#3399CC") # was #3399CC
            end
        end
        push!(plots_vector, p_param)
    end

    # make plot
    p_param_boxes = Plots.plot(plots_vector..., 
        layout=length(plots_vector), size=plot_size, legend=false,
        suptitle=suptitle,  plot_titlefontfamily="Bookman Light", guidefontsize=4, tickfontsize=8
        )

    return p_param_boxes

end