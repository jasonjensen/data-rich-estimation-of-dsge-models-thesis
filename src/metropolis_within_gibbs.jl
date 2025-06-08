"""
    metropolis_within_gibbs(m::MacroModelling.ℳ, data=Workspace(), Σ⁻¹::Union{Symmetric,Nothing}=nothing; niter=100, burnin=10, c=0.5, w=0.02, 
        target_acceptance_rate=0.25, 
        Λ_constraints = Dict{String, Pair{Symbol,Symbol}}(),
        c_min = w,
        verbose = false,
        fixed_parameters = Vector{Symbol}(),
        restore_from = nothing,
        restore_from_current_iter=0,
        speedup=false,
        measurement_error = true,
        progress_plots=false,
    )

Run the Metropolis-within-Gibbs sampler for model estimation.

# Arguments
- `m::MacroModelling.ℳ`: Model object.
- `data=Workspace()`: Workspace containing the data.
- `Σ⁻¹::Union{Symmetric,Nothing}=nothing`: Inverse of the covariance matrix for the proposal distribution.
- `niter=100`: Number of iterations.
- `burnin=10`: Number of burn-in iterations.
- `c=0.5`: Initial scaling factor for the proposal distribution.
- `w=0.02`: Step size for adjusting `c`.
- `target_acceptance_rate=0.25`: Target acceptance rate for the Metropolis-Hastings step.
- `Λ_constraints = Dict{String, Pair{Symbol,Symbol}}()`: Dictionary of constraints on the observation matrix.
- `c_min = w`: Minimum value for `c`.
- `verbose = false`: Whether to print verbose output.
- `fixed_parameters = Vector{Symbol}()`: Vector of fixed parameters.
- `restore_from = nothing`: Path to a file to restore the sampler state from.
- `restore_from_current_iter=0`: Iteration to restore from if `restore_from` is specified.
- `speedup=false`: Whether to use speedup techniques.
- `measurement_error = true`: Whether to estimate measurement error.
- `progress_plots=false`: Whether to generate progress plots.

# Returns
- `ret::Workspace`: Workspace containing the MCMC chain and other results.
"""
function metropolis_within_gibbs(m::MacroModelling.ℳ, data=Workspace(), Σ⁻¹::Union{Symmetric,Nothing}=nothing; niter=100, burnin=10, c=0.5, w=0.02, 
    target_acceptance_rate=0.25, 
    Λ_constraints = Dict{String, Pair{Symbol,Symbol}}(),
    c_min = w,
    verbose = false,
    fixed_parameters = Vector{Symbol}(),
    restore_from = nothing,
    restore_from_current_iter=0,
    speedup=false,
    measurement_error = true,
    progress_plots=false,
)

    

    local S
    local R_Θ
    local Λ_μ
    local 𝛙_μ
    local 𝛙_σ

    # initialization
    θ = copy(m.parameter_values)
    TT, SS_and_pars, 𝐒mat, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)

    P₀=get_P₀(m, 𝐒mat)
    M = get_measurement_matrix(m, 𝐒mat, data, core=true)
    A = get_measurement_const(m, θ, data)
    # SR = Diagonal([(0.2*0.579923)^2, (0.2*1.470832)^2, (0.2*2.237937)^2]) # here is this from?!?!?!
    SR = I(length(data.core_series));
    Q = get_Q(m, θ)
    Q = I(length(m.exo))
    # should return the full set of states
    S, _ = kalman(m, 𝐒mat, data, M, SR, A, Q; core=true, backward=true)
    
    R_Θ, Λ_μ, 𝛙_μ, 𝛙_σ = get_initial_Γ(S, data, Λ_constraints)
    
    # iterate a few times:
    for i = 1:10
        S, loglikX = generate_states(m, θ, 𝐒mat, data, Λ_μ, 𝛙_μ, R_Θ)
        R_Θ, Λ_μ, 𝛙_μ, 𝛙_σ = get_initial_Γ(S, data, Λ_constraints)
    end

    nAccepted = 0
    𝛉 = zeros(niter, length(m.parameters))
    𝛉_proposed = zeros(niter, length(m.parameters))
    𝐜 = zeros(niter)
    𝛚 = zeros(niter)
    𝐋 = zeros(niter)
    𝐋𝛉= zeros(niter)
    𝐋𝚪= zeros(niter)
    𝚲draws = repeat([zeros(size(Λ_μ))], niter) 
    Rdraws = repeat([zeros(size(R_Θ))], niter) 
    𝛙draws = repeat([zeros(size(𝛙_μ))], niter) 
    accepted_record = zeros(niter)
    acceptance_rate = 1.0
    rerolled_loglikelihood = falses(niter)
    c_move = "stable"

    _w = w
    
    # iter = 1:niter
    loglikθ = get_log_likelihood_θ(m, fixed_parameters)
    loglikΓ = 0
    loglik = get_log_likelihood_X(m, Matrix(data.X[2:end,:]), S[:,1:end-1], Λ_μ, 𝛙_μ, R_Θ, 2) + loglikθ
    highest_loglik = -Inf
    num_errors = 0
    Λdraw = Λ_μ
    Rdraw = R_Θ
    𝛙draw = 𝛙_μ
    ret = Workspace()
    start_time = Dates.format(Dates.now(),"yyyy-mm-ddTHHMM")
    ret.max_iterations = niter
    ret.burnin = burnin
    current_iter = 0
    sol_guess_init = copy(m.solution.perturbation.qme_solution)
    local_acceptance_rate = 0.0;
    short_accepted_record = falses(100)

    

    if restore_from !== nothing
        current_iter, 𝛉, 𝐜, 𝛚, 𝐋, 𝐋𝛉, 𝐋𝚪, 𝛉_proposed, 𝚲draws, Rdraws, 𝛙draws, burnin, niter, Λdraw, Rdraw, 𝛙draw, loglikθ, loglikΓ, loglik, start_time, c, Λ_constraints, SR, fixed_parameters, S, accepted_record, nAccepted, θ, num_errors, speedup, short_accepted_record, measurement_error = restore_metropolis_within_gibbs(restore_from, restore_from_current_iter)
        m.parameter_values = 𝛉[current_iter,:]
        TT, SS_and_pars, 𝐒mat, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    end

    iter = ProgressBar(current_iter+1:niter)
    last_error = ""
    last_accepted = 1
    next_to_last_accepted = 1
    last_Smat = 𝐒mat
    next_to_last_Smat = 𝐒mat
    Σ⁻¹actual = modifyΣ(Σ⁻¹, m, θ)

    for i ∈ iter
        f = metropolis_within_gibbs_step
        accepted = false
        ω = 0
        
        try
            # i == 1 && @warn "NO TRY/CATCH"
            accepted, θ, 𝐒mat, S, ω, loglik,loglikΓ, loglikθ, θdraw, Λdraw, Rdraw, 𝛙draw = f(
                    m, θ, 𝐒mat, S, data, Σ⁻¹actual, 
                    𝛙draw, 
                    𝛙_σ, 
                    Λ_constraints, c, loglikθ, loglik,
                    Λdraw, 
                    Rdraw, 
                    𝛙draw;
                    auto_accept=(i==1),
                    verbose = verbose,
                    fixed_parameters=fixed_parameters,
                    speedup=speedup,
                    measurement_error=measurement_error,
                    P₀=P₀,
                )
                
        catch err
            last_error = typeof(err)
            # rethrow(err)
            # m.solution.perturbation.qme_solution = copy(sol_guess_init)
            num_errors += 1
            θ = 𝛉[next_to_last_accepted, :]
            Rdraw = Rdraws[next_to_last_accepted]
            𝛙draw = 𝛙draws[next_to_last_accepted] 
            𝛙draw = 𝛙draws[next_to_last_accepted] 
            𝐒mat = next_to_last_Smat
            accepted = false
        end
        
        if accepted
            nAccepted += 1
            accepted_record[i] = 1.0
            next_to_last_accepted = last_accepted
            last_accepted = i
            next_to_last_Smat = last_Smat
            last_Smat = 𝐒mat
        end
        𝐋[i] = loglik
        𝐋𝛉[i] = loglikθ
        𝐋𝚪[i] = loglikΓ
        highest_loglik = max(highest_loglik, loglik)
        popfirst!(short_accepted_record)
        push!(short_accepted_record, accepted)
       
        acceptance_rate = round(nAccepted/i; digits=2)
        
        
        if i % 10 == 0 && i > burnin 
            GC.gc()
            if i > 100
                local_acceptance_rate = round(mean(short_accepted_record); digits=2)
            else
                rng = 1:i
                local_acceptance_rate = round(sum(accepted_record[rng])/length(rng); digits=2)
            end

            # only adjust c after burnin, and every once-in-a-while
            if local_acceptance_rate > target_acceptance_rate #&& c < 100.0 - _w
                if c < _w
                    c = c*1.5
                    c_move = "growing (x1.5)"
                else
                    c_move = "growing ($_w)"
                    c += _w
                end
            elseif local_acceptance_rate < target_acceptance_rate
                if c - _w < _w
                    c_move = "shrinking (x0.7)"
                    c = c / 1.3
                else
                    c_move = "shrinking ($_w)"
                    c = max(c - _w, c_min)
                end
                if c < c_min#
                    # expand a little
                    c_move = "growing (x3)"
                    c = c_min*3
                end
            else
                c_move = "stable"
            end
            
            # decrease _w
            _w = w*((niter-i)/niter)

            # update Σ⁻¹
            Σ⁻¹actual = modifyΣ(Σ⁻¹, m, θ)
        end
        
        set_multiline_postfix(iter, "c: $(c)\nc move: $(c_move)\nAcceptance rate: $(local_acceptance_rate)\nLogLik: $(loglik)\nTotal accepted: $(nAccepted)\nHighest loglik: $(highest_loglik)\nnum_errors: $(num_errors)\nlast_error: $(last_error)\n ")
        
        𝐜[i] = c
        𝛚[i] = ω
        𝛉[i, :] = θ
        𝛉_proposed[i, :] = θ
        𝚲draws[i] = Λdraw
        Rdraws[i] = Rdraw 
        𝛙draws[i] = 𝛙draw
        
        # store progress, make charts
        if i % 5000 == 0
            ret.𝛉 = 𝛉
            ret.𝐜 = 𝐜
            ret.𝛚 = 𝛚
            ret.𝐋 = 𝐋
            ret.𝐋𝛉 = 𝐋𝛉
            ret.𝐋𝚪 = 𝐋𝚪
            ret.𝛉_proposed = 𝛉_proposed
            ret.𝚲draws = 𝚲draws
            ret.Rdraws = Rdraws
            ret.𝛙draws = 𝛙draws
            ret.burnin = burnin
            ret.current_iteration = i
            ret.start_time = start_time
            ret.Λ_constraints = Λ_constraints
            ret.SR=SR
            ret.fixed_parameters=fixed_parameters
            ret.S = S
            ret.accepted_record = accepted_record
            ret.num_errors = num_errors
            ret.speedup = speedup
            ret.measurement_error = measurement_error
            jldsave("""data/results/partial_gelfer_$(start_time).jld2"""; ret)

            if progress_plots
                make_progress_plots(m, 𝛉, start_time, i)
            end

        end
    end
    println("Overall acceptance rate: $(round(nAccepted/niter; digits=2))")
    println("Errors: $(num_errors)")
   
    
    ret.𝛉 = 𝛉
    ret.𝐜 = 𝐜
    ret.𝛚 = 𝛚
    ret.𝐋 = 𝐋
    ret.𝐋𝛉 = 𝐋𝛉
    ret.𝐋𝚪 = 𝐋𝚪
    ret.𝛉_proposed = 𝛉_proposed
    ret.𝚲draws = 𝚲draws
    ret.Rdraws = Rdraws
    ret.𝛙draws = 𝛙draws
    ret.burnin = burnin
    ret.w = w
    ret._w = _w
    ret.start_time = start_time
    ret.current_iteration = niter
    ret.speedup = speedup
    ret.measurement_error = measurement_error

    m.parameter_values = θ
    
    return m, ret
end

"""
    metropolis_within_gibbs_step(m, θ, 𝐒mat, S, data,  Σ⁻¹, 𝛙_μ  , 𝛙_σ, Λ_constraints, c, loglikθ, prev_loglik,
        Λ_prev=nothing, R_prev=nothing, 𝛙_prev=nothing; 
        auto_accept=false,
        verbose = false,
        fixed_parameters = Vector{Symbol}(),
        speedup=false,
        measurement_error = true,
        P₀=P₀,
    )

Perform a single step of the Metropolis-within-Gibbs sampler.

# Arguments
- `m`: Model object.
- `θ`: Current parameter vector.
- `𝐒mat`: Current state-space matrices.
- `S`: Current state vector.
- `data`: Workspace containing the data.
- `Σ⁻¹`: Inverse of the covariance matrix for the proposal distribution.
- `𝛙_μ`: Mean of the measurement error autocorrelation.
- `𝛙_σ`: Standard deviation of the measurement error autocorrelation.
- `Λ_constraints`: Dictionary of constraints on the observation matrix.
- `c`: Scaling factor for the proposal distribution.
- `loglikθ`: Log-likelihood of the parameters.
- `prev_loglik`: Previous log-likelihood.
- `Λ_prev=nothing`: Previous observation matrix.
- `R_prev=nothing`: Previous measurement error covariance matrix.
- `𝛙_prev=nothing`: Previous measurement error autocorrelation matrix.
- `auto_accept=false`: Whether to automatically accept the proposal.
- `verbose = false`: Whether to print verbose output.
- `fixed_parameters = Vector{Symbol}()`: Vector of fixed parameters.
- `speedup=false`: Whether to use speedup techniques.
- `measurement_error = true`: Whether to estimate measurement error.
- `P₀=P₀`: Initial state covariance matrix.

# Returns
- A tuple containing: `accepted`, `θ`, `𝐒mat`, `S`, `ω`, `loglik`, `loglikΓ`, `loglikθ`, `θdraw`, `Λdraw`, `Rdraw`, `𝛙draw`.
"""
function metropolis_within_gibbs_step(m, θ, 𝐒mat, S, data,  Σ⁻¹, 𝛙_μ  , 𝛙_σ, Λ_constraints, c, loglikθ, prev_loglik,
    Λ_prev=nothing, R_prev=nothing, 𝛙_prev=nothing; 
        auto_accept=false,
        verbose = false,
        fixed_parameters = Vector{Symbol}(),
        speedup=false,
        measurement_error = true,
        P₀=P₀,
    )
    # draw new states every iteration
    
    @timeit to "generate_states" Sᵍ, loglikX, _sumSpred = generate_states(m, θ, 𝐒mat, data, Λ_prev, 𝛙_prev, R_prev; speedup=speedup,P₀=P₀)
    if abs(_sumSpred) < 1e-9
        return false, θ, 𝐒mat, S, 0.0, prev_loglik, -Inf, loglikθ, θ, Λ_prev, R_prev, 𝛙_prev
    end
    # hold on to the solution guess
    sol_guess = copy(m.solution.perturbation.qme_solution)

    @timeit to "draw_Γ" Λ, 𝛙, R, loglikΓ = draw_Γ(Sᵍ, data, Diagonal(𝛙_μ), 𝛙_σ, Λ_constraints, measurement_error, verbose)
    if isinf(loglikΓ) 
        return false, θ, 𝐒mat, Sᵍ, 0.0, prev_loglik, loglikΓ, loglikθ, θ, Λ_prev, R_prev, 𝛙_prev
    end
    
    @timeit to "perturb_θ" _θ, _loglikθ, _𝐒mat, solved = perturb_θ(m, θ, Σ⁻¹, c, fixed_parameters=fixed_parameters)
    if !solved > 0 
        m.solution.perturbation.qme_solution = sol_guess
        return false, θ, 𝐒mat, Sᵍ, 0.0, prev_loglik, loglikΓ, loglikθ, θ, Λ, R, 𝛙
    end

    @timeit to "generate_states" _S, _loglikX = generate_states(m, _θ, _𝐒mat, data, Λ, 𝛙, R; speedup=speedup, backwards_pass=false, P₀=P₀)
    @timeit to "generate_states" Sᵍ,  loglikX = generate_states(m, θ, 𝐒mat, data, Λ, 𝛙, R; speedup=speedup, backwards_pass=false,P₀=P₀)
    
    
    _loglik = _loglikX  + _loglikθ
    loglik = loglikX + loglikθ
    
    log_ratio = _loglik - loglik
    ω = max(min(exp(log_ratio), 1), 0)
    if verbose
        println("ω: $ω")
        println("r: $log_ratio")
        println("========================")
        println("θ: $_loglikθ vs $loglikθ")
        println("Γ: $_loglikΓ vs $loglikΓ")
        println("X: $_loglikX vs $loglikX")
        println("T: $(_loglik) vs $(loglik)")
        println("========================")
        println("========================")
        println("========================")
        println("========================")
        println("========================")
        println("========================")
    end
    accept_prob = rand(Uniform(0,1))
    if ω > accept_prob || auto_accept == true
        return true, _θ, _𝐒mat, _S, ω, _loglik, 0, _loglikθ,  _θ, Λ, R, 𝛙
    end
    
    return false, θ, 𝐒mat, Sᵍ, ω, loglik, loglikΓ, loglikθ, _θ, Λ, R, 𝛙
        
   
end

"""
    restore_metropolis_within_gibbs(restore_from, restore_from_current_iter)

Restore the state of a Metropolis-within-Gibbs sampler from a saved JLD2 file.

This function loads a previously saved `Workspace` (assumed to be stored under the key "ret"
in the JLD2 file) containing the state of an MCMC estimation. It then reconstructs
the necessary variables to resume the `metropolis_within_gibbs` sampler from a specified
iteration.

# Arguments
- `restore_from::String`: The file path to the JLD2 file containing the saved sampler state.
- `restore_from_current_iter::Int`: The iteration number from which to restore the sampler.
                                   If `0`, it restores from the last saved iteration in the file.

# Returns
- `current_iter::Int`: The iteration number from which the sampler is restored.
- `𝛉::Matrix{Float64}`: Matrix of all parameter draws up to the last saved iteration.
- `𝐜::Vector{Float64}`: Vector of scaling factors `c` for each iteration.
- `𝛚::Vector{Float64}`: Vector of acceptance probabilities `ω` for each iteration.
- `𝐋::Vector{Float64}`: Vector of total log-likelihoods for each iteration.
- `𝐋𝛉::Vector{Float64}`: Vector of parameter log-likelihoods (priors) for each iteration.
- `𝐋𝚪::Vector{Float64}`: Vector of measurement equation parameter log-likelihoods for each iteration.
- `𝛉_proposed::Matrix{Float64}`: Matrix of proposed parameter draws for each iteration.
- `𝚲draws::Vector{Matrix{Float64}}`: Vector of drawn observation matrices `Λ` for each iteration.
- `Rdraws::Vector{Diagonal{Float64}}`: Vector of drawn measurement error covariance matrices `R` for each iteration.
- `𝛙draws::Vector{Diagonal{Float64}}`: Vector of drawn measurement error autocorrelation matrices `ψ` for each iteration.
- `burnin::Int`: The burn-in period used in the saved sampler.
- `niter::Int`: The total number of iterations planned for the saved sampler.
- `Λdraw::Matrix{Float64}`: The observation matrix `Λ` at the `current_iter`.
- `Rdraw::Diagonal{Float64}`: The measurement error covariance `R` at the `current_iter`.
- `𝛙draw::Diagonal{Float64}`: The measurement error autocorrelation `ψ` at the `current_iter`.
- `loglikθ::Float64`: The parameter log-likelihood at `current_iter`.
- `loglikΓ::Float64`: The measurement equation parameter log-likelihood at `current_iter`.
- `loglik::Float64`: The total log-likelihood at `current_iter`.
- `start_time::String`: The start timestamp of the original saved sampler.
- `c::Float64`: The scaling factor `c` at `current_iter`.
- `Λ_constraints`: Constraints on the observation matrix from the saved sampler.
- `SR`: Measurement error covariance for core series from the saved sampler (likely fixed).
- `fixed_parameters::Vector{Symbol}`: List of fixed parameters from the saved sampler.
- `S`: State vector draws from the saved sampler, corresponding to `current_iter`.
- `accepted_record::Vector{Float64}`: Record of accepted proposals (1.0 for accepted, 0.0 otherwise).
- `nAccepted::Float64`: Total number of accepted proposals up to `current_iter`.
- `θ::Vector{Float64}`: Parameter vector at `current_iter`.
- `num_errors::Int`: Number of errors encountered in the saved sampler.
- `speedup::Bool`: Speedup flag from the saved sampler.
- `short_accepted_record::Vector{Bool}`: Record of acceptance for the last 100 iterations.
- `measurement_error::Bool`: Flag indicating if measurement error was estimated in the saved sampler.

This function also updates the global model `m` with the parameters from the restored iteration
and re-solves the model to get the state-space representation `𝐒mat`.
"""
function restore_metropolis_within_gibbs(restore_from, restore_from_current_iter)
    ret = load(restore_from)["ret"]

    current_iter = ret.current_iteration
    if restore_from_current_iter !== 0
        current_iter = restore_from_current_iter
    end

    𝛉 = ret.𝛉
    𝐜 = ret.𝐜
    𝛚 = ret.𝛚
    𝐋 = ret.𝐋
    𝐋𝛉 = ret.𝐋𝛉
    𝐋𝚪 = ret.𝐋𝚪
    𝛉_proposed = ret.𝛉_proposed
    𝚲draws = ret.𝚲draws
    Rdraws = ret.Rdraws
    𝛙draws = ret.𝛙draws
    burnin = ret.burnin
    niter = ret.max_iterations
    
    
    Λdraw = 𝚲draws[current_iter]
    Rdraw = Rdraws[current_iter]
    𝛙draw = 𝛙draws[current_iter]
    loglikθ = 𝐋𝛉[current_iter] 
    loglikΓ = 𝐋𝚪[current_iter]
    loglik = 𝐋[current_iter]
    start_time = ret.start_time
    c = ret.𝐜[current_iter]
    Λ_constraints = ret.Λ_constraints
    SR=ret.SR
    fixed_parameters=ret.fixed_parameters
    S = ret.S
    m.parameter_values = 𝛉[current_iter,:]
    accepted_record = ret.accepted_record
    nAccepted = sum(accepted_record)
    θ = ret.𝛉[current_iter, :]
    num_errors = ret.num_errors
    speedup = ret.speedup
    TT, SS_and_pars, 𝐒mat, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    short_accepted_record = Bool.(accepted_record[end-99:end])
    measurement_error = ret.measurement_error

    return (current_iter, 𝛉, 𝐜, 𝛚, 𝐋, 𝐋𝛉, 𝐋𝚪, 𝛉_proposed, 𝚲draws, Rdraws, 𝛙draws, burnin, niter, Λdraw, Rdraw, 𝛙draw, loglikθ, loglikΓ, loglik, start_time, c, Λ_constraints, SR, 
    fixed_parameters, S, accepted_record, nAccepted, θ, num_errors, speedup, short_accepted_record, measurement_error)  
end

"""
    make_progress_plots(m::MacroModelling.ℳ, 𝛉::Matrix, start_time::String)

Generate and save progress plots during the MCMC estimation.

This function creates two types of plots at intervals during the `metropolis_within_gibbs` run:
1.  **Parameter Trace Plots**: For each estimated parameter in the model (`m.parameters`),
    it plots the history of its sampled values from the MCMC chain `𝛉` up to the
    current iteration `i`. These plots help visualize the mixing and convergence of
    each parameter. The plots are saved to `graphs/params_progress_{start_time}.png`.
2.  **Impulse Response Function (IRF) Plots**: For a specific shock (e.g., `:ϵr`), it plots
    IRFs using:
    a. The model's default parameters (e.g., `SW.parameter_values`).
    b. Pre-defined "working parameters" (e.g., `working_params_sw_dfm`).
    c. The current parameter draw `𝛉[i,:]` from the MCMC chain.
    These are saved to `graphs/irfs_progress_{start_time}.png`.
    If the number of iterations `i` is greater than 50,000, it also plots IRFs using the
    mean of the parameter draws from iteration 50,000 to `i`, saved to
    `graphs/irfs_progress_mean_{start_time}.png`.

The specific parameters plotted and the layout depend on `m.model_name` (either "SW" or "SWFF").
Plotting uses a `plot_step` to avoid plotting every single draw, for efficiency.

# Arguments
- `m::MacroModelling.ℳ`: The model object.
- `𝛉::Matrix`: A matrix containing the MCMC draws for all parameters up to the current
                 iteration. Dimensions are `(current_iteration, num_parameters)`.
                 The function internally accesses the current iteration `i` (implicitly, as this function is called within a loop over `i`).
- `start_time::String`: A timestamp string used in naming the output plot files,
                        typically indicating when the MCMC estimation started.
"""
function make_progress_plots(m::MacroModelling.ℳ, 𝛉, start_time, i)
    # plot stuff
    if m.model_name == "SW"
        working_params_sw_dfm = [[gelfer_means_sw_dfm_estimation[v][1] for v ∈ m.parameters[1:29]]..., m.parameter_values[30:end]...]
        working_params_sw_dfm[findfirst(x -> x == :φ, m.parameters)] = SW.parameter_values[findfirst(x -> x == :φ, m.parameters)]
    elseif m.model_name == "SWFF"
        working_params_swff_dfm = [[gelfer_means_swff_dfm_estimation[v][1] for v ∈ m.parameters[1:31]]..., m.parameter_values[32:end]...]
    end

    plot_step=10
    if m.model_name == "SW"
        p_param_progress = Plots.plot(
            Plots.plot(𝛉[1:plot_step:i,1], title=m.parameters[1]), 
            Plots.plot(𝛉[1:plot_step:i,2], title=m.parameters[2]), 
            Plots.plot(𝛉[1:plot_step:i,3], title=m.parameters[3]), 
            Plots.plot(𝛉[1:plot_step:i,4], title=m.parameters[4]), 
            Plots.plot(𝛉[1:plot_step:i,5], title=m.parameters[5]), 
            Plots.plot(𝛉[1:plot_step:i,6], title=m.parameters[6]), 
            Plots.plot(𝛉[1:plot_step:i,7], title=m.parameters[7]), 
            Plots.plot(𝛉[1:plot_step:i,8], title=m.parameters[8]), 
            Plots.plot(𝛉[1:plot_step:i,9], title=m.parameters[9]), 
            Plots.plot(𝛉[1:plot_step:i,10], title=m.parameters[10]), 
            Plots.plot(𝛉[1:plot_step:i,11], title=m.parameters[11]), 
            Plots.plot(𝛉[1:plot_step:i,12], title=m.parameters[12]), 
            Plots.plot(𝛉[1:plot_step:i,13], title=m.parameters[13]), 
            Plots.plot(𝛉[1:plot_step:i,14], title=m.parameters[14]), 
            Plots.plot(𝛉[1:plot_step:i,15], title=m.parameters[15]), 
            Plots.plot(𝛉[1:plot_step:i,16], title=m.parameters[16]), 
            Plots.plot(𝛉[1:plot_step:i,17], title=m.parameters[17]), 
            Plots.plot(𝛉[1:plot_step:i,18], title=m.parameters[18]), 
            Plots.plot(𝛉[1:plot_step:i,19], title=m.parameters[19]),
            Plots.plot(𝛉[1:plot_step:i,20], title=m.parameters[20]), 
            Plots.plot(𝛉[1:plot_step:i,21], title=m.parameters[21]), 
            Plots.plot(𝛉[1:plot_step:i,22], title=m.parameters[22]), 
            Plots.plot(𝛉[1:plot_step:i,23], title=m.parameters[23]), 
            Plots.plot(𝛉[1:plot_step:i,24], title=m.parameters[24]), 
            Plots.plot(𝛉[1:plot_step:i,25], title=m.parameters[25]), 
            Plots.plot(𝛉[1:plot_step:i,26], title=m.parameters[26]), 
            Plots.plot(𝛉[1:plot_step:i,27], title=m.parameters[27]), 
            Plots.plot(𝛉[1:plot_step:i,28], title=m.parameters[28]), 
            Plots.plot(𝛉[1:plot_step:i,29], title=m.parameters[29]),
        layout=(6, 5), legend=false, size=(1000,1000));
        Plots.savefig(p_param_progress, "graphs/params_progress_$(start_time).png");

        try
            pϵr1 = plot_irfs(m, 
                [SW.parameter_values, working_params_sw_dfm, 𝛉[i,:]],
                ["Reg", "DFM", "Est"],
                :ϵr)
            # pϵr1 = plot_irfs(m, SW.parameter_values, 𝛉[i,:], :ϵr);
            Plots.savefig(pϵr1, "graphs/irfs_progress_$(start_time).png");
            if i > 50000
                # pϵr2 = plot_irfs(m, SW.parameter_values, mean(𝛉[50000:i,:], dims=1)[:], :ϵr)
                pϵr2 = plot_irfs(m, 
                    [SW.parameter_values, working_params_sw_dfm, mean(𝛉[50000:end,:], dims=1)[:]],
                    ["Reg", "DFM", "Est"],
                    :ϵr)
                Plots.savefig(pϵr2, "graphs/irfs_progress_mean_$(start_time).png");
            end
        catch
        end
    end
    if m.model_name == "SWFF"
        p_param_progress = Plots.plot(
            Plots.plot(𝛉[1:plot_step:i,1], title=m.parameters[1]), 
            Plots.plot(𝛉[1:plot_step:i,2], title=m.parameters[2]), 
            Plots.plot(𝛉[1:plot_step:i,3], title=m.parameters[3]), 
            Plots.plot(𝛉[1:plot_step:i,4], title=m.parameters[4]), 
            Plots.plot(𝛉[1:plot_step:i,5], title=m.parameters[5]), 
            Plots.plot(𝛉[1:plot_step:i,6], title=m.parameters[6]), 
            Plots.plot(𝛉[1:plot_step:i,7], title=m.parameters[7]), 
            Plots.plot(𝛉[1:plot_step:i,8], title=m.parameters[8]), 
            Plots.plot(𝛉[1:plot_step:i,9], title=m.parameters[9]), 
            Plots.plot(𝛉[1:plot_step:i,10], title=m.parameters[10]), 
            Plots.plot(𝛉[1:plot_step:i,11], title=m.parameters[11]), 
            Plots.plot(𝛉[1:plot_step:i,12], title=m.parameters[12]), 
            Plots.plot(𝛉[1:plot_step:i,13], title=m.parameters[13]), 
            Plots.plot(𝛉[1:plot_step:i,14], title=m.parameters[14]), 
            Plots.plot(𝛉[1:plot_step:i,15], title=m.parameters[15]), 
            Plots.plot(𝛉[1:plot_step:i,16], title=m.parameters[16]), 
            Plots.plot(𝛉[1:plot_step:i,17], title=m.parameters[17]), 
            Plots.plot(𝛉[1:plot_step:i,18], title=m.parameters[18]), 
            Plots.plot(𝛉[1:plot_step:i,19], title=m.parameters[19]),
            Plots.plot(𝛉[1:plot_step:i,20], title=m.parameters[20]), 
            Plots.plot(𝛉[1:plot_step:i,21], title=m.parameters[21]), 
            Plots.plot(𝛉[1:plot_step:i,22], title=m.parameters[22]), 
            Plots.plot(𝛉[1:plot_step:i,23], title=m.parameters[23]), 
            Plots.plot(𝛉[1:plot_step:i,24], title=m.parameters[24]), 
            Plots.plot(𝛉[1:plot_step:i,25], title=m.parameters[25]), 
            Plots.plot(𝛉[1:plot_step:i,26], title=m.parameters[26]), 
            Plots.plot(𝛉[1:plot_step:i,27], title=m.parameters[27]), 
            Plots.plot(𝛉[1:plot_step:i,28], title=m.parameters[28]), 
            Plots.plot(𝛉[1:plot_step:i,29], title=m.parameters[29]),
            Plots.plot(𝛉[1:plot_step:i,30], title=m.parameters[30]),
        layout=(6, 5), legend=false, size=(1000,1000));
        Plots.savefig(p_param_progress, "graphs/params_progress_$(start_time).png");

        try
            pϵr1 = plot_irfs(m, 
                [SWFF.parameter_values, working_params_swff_dfm, 𝛉[i,:]],
                ["Reg", "DFM", "Est"],
                :ϵr)
            # pϵr1 = plot_irfs(m, SW.parameter_values, 𝛉[i,:], :ϵr);
            Plots.savefig(pϵr1, "graphs/irfs_progress_$(start_time).png");
            if i > 50000
                # pϵr2 = plot_irfs(m, SW.parameter_values, mean(𝛉[50000:i,:], dims=1)[:], :ϵr)
                pϵr2 = plot_irfs(m, 
                    [SWFF.parameter_values, working_params_swff_dfm, mean(𝛉[50000:end,:], dims=1)[:]],
                    ["Reg", "DFM", "Est"],
                    :ϵr)
                Plots.savefig(pϵr2, "graphs/irfs_progress_mean_$(start_time).png");
            end
        catch
        end
    end
end

"""
    get_jacobian(m, θ, Λ, Ψ, R, data, fixed_parameters, δ=0.005, speedup=false)

Calculate the Jacobian of the log-likelihood function.

# Arguments
- `m`: Model object.
- `θ`: Parameter vector.
- `Λ`: Observation matrix.
- `Ψ`: Autocorrelation matrix for measurement errors.
- `R`: Covariance matrix for measurement errors.
- `data`: Workspace containing the data.
- `fixed_parameters`: Vector of fixed parameters.
- `δ=0.005`: Step size for finite differences.
- `speedup=false`: Whether to use speedup techniques.

# Returns
- `J`: Jacobian vector.
"""
function get_jacobian(m, θ, Λ, Ψ, R, data, fixed_parameters, δ=0.005, speedup=false)
    J = zeros(length(m.parameters))
    loglike_orig = 0
    _, _, 𝐒, _, _ = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    # M = get_measurement_matrix(m, 𝐒, data, core=true)
    # A = get_measurement_const(m, θ, data)
    # SR = I(length(data.core_series)) # not sure about this!
    # Q = I(length(m.exo))
    # should return the full set of states
    
    # Threads.@threads for i ∈ 1:n_draws
    # _S, loglike_orig = kalman(m, 𝐒, data, M, SR, A, Q; core=true, backward=false)
    _S, loglike_orig = generate_states(m, θ, 𝐒, data, Λ, Ψ, R, backwards_pass=false, speedup=speedup)
    @show loglike_orig
    # @show δ

    for i ∈ 1:length(J)
        if m.parameters[i] ∈ fixed_parameters
            continue
        end
        _δ = copy(δ)
        _θ = copy(θ)
        _θ[i] = θ[i] + _δ
        while clamp_θ(Val(Symbol(m.model_name)), Val(m.parameters[i]), _θ[i]) != _θ[i]
            _δ = _δ/2
            _θ[i] = θ[i] + _δ
        end
        _, _, _𝐒, _, _ = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), _θ, m)
        # _M = get_measurement_matrix(m, _𝐒, data, core=true)
        # _A = get_measurement_const(m, _θ, data)
        # _S, loglike_param = kalman(m, _𝐒, data, _M, SR, _A, Q; core=true, backward=false)
        _S, loglike_param = generate_states(m, _θ, _𝐒, data, Λ, Ψ, R, backwards_pass=false, speedup=speedup )

        dLLdθ = (loglike_param - loglike_orig)/_δ
        J[i] = dLLdθ
    end
    
    return J
end

"""
    get_hessian(m, θ, Λ, Ψ, R, data, fixed_parameters, δ = 0.025, speedup=false)

Calculate the Hessian of the log-likelihood function.

# Arguments
- `m`: Model object.
- `θ`: Parameter vector.
- `Λ`: Observation matrix.
- `Ψ`: Autocorrelation matrix for measurement errors.
- `R`: Covariance matrix for measurement errors.
- `data`: Workspace containing the data.
- `fixed_parameters`: Vector of fixed parameters.
- `δ = 0.025`: Step size for finite differences.
- `speedup=false`: Whether to use speedup techniques.

# Returns
- `H`: Hessian matrix.
"""
function get_hessian(m, θ, Λ, Ψ, R, data, fixed_parameters, δ = 0.025, speedup=false)
    init_δ = copy(δ)
    H = zeros(length(m.parameters),length(m.parameters))
    J = get_jacobian(m, θ, Λ, Ψ, R, data, fixed_parameters, init_δ, speedup)
    # the k'th entry in J is how much the LokLik goes up when the k'th parameter is increased
    # J[k] is dL/dθₖ
    for i ∈ 1:length(m.parameters)
        if m.parameters[i] ∈ fixed_parameters
            continue
        end
        _δ = copy(δ)
        _θ = deepcopy(θ)
        _θ[i] = θ[i] + _δ
        while clamp_θ(Val(Symbol(m.model_name)), Val(m.parameters[i]), _θ[i]) != _θ[i]
            _δ = _δ/2
            _θ[i] = θ[i] + _δ
        end
        _J = get_jacobian(m, _θ, Λ, Ψ, R, data, fixed_parameters, init_δ, speedup)
        # _J[k] is dL\dθₖ
        # the k'th entry in _J is how much the LogLik goes up when the k-th parmeter is increased, conditional on having increased the i'th parameter
        for j ∈ 1:length(m.parameters)
            @show i, j
            if m.parameters[j] ∈ fixed_parameters
                continue
            end
            # the order of i,j doesn't matter here because it's symmetric.
            # when I change parameter i, how much does the derivative (of the loglik) with respects to parameter j change
            # H[i,j] = (J[j] - _J[j]) / _δ
            # NOTE: the correct order of the subtraction has been reversed here because we want the negative of the Hessian values
            H[i,j] = (J[j] - _J[j]) / _δ
        end
        display(H)
        if sum(H) == 0
            error("Loglikelihood not changing")
        end
    end
    return H
    

end

"""
    modifyΣ(v::Val{:SW}, m, Σ, θ)
    modifyΣ(v::Val{:SWFF}, m, Σ, θ)
    modifyΣ(Σ⁻¹, m, θ)

Modify the proposal covariance matrix `Σ` (or its inverse `Σ⁻¹`) based on current parameter values `θ`.

This function adjusts the proposal covariance matrix used in the Metropolis-Hastings algorithm.
The primary goal is to improve mixing and convergence by reducing or eliminating
cross-correlations in the proposal distribution when parameters are near their bounds.

For specific parameters (e.g., autoregressive coefficients like `:ρ`, `:ρa`, or parameters
like `:h`, `:νl`, `:φ`, `:rπ1`), if their current value in `θ` is close to a predefined
boundary (e.g., > 0.95 or < 0.05 for AR coefficients), the corresponding off-diagonal
elements (covariances) in `Σ` are set to zero. The diagonal element (variance)
is preserved.

After modifications, the function ensures that the resulting matrix `_Σ` is
positive semi-definite (PSD) by calling `make_psd(_Σ)` if necessary.

The function is dispatched based on the model type (`Val{:SW}` or `Val{:SWFF}`) or can
accept `Σ⁻¹` directly. If `Σ⁻¹` is passed, it's assumed to be the inverse of the
covariance matrix and similar logic is applied.

# Arguments
- `v::Val`: A `Val` type specifying the model (e.g., `Val(:SW)` or `Val(:SWFF)`). Used for dispatch.
- `m`: The model object (`MacroModelling.ℳ`).
- `Σ`: The proposal covariance matrix (or `Σ⁻¹` if passing the inverse).
- `θ`: The current vector of parameter values.

# Returns
- `_Σ::Symmetric`: The modified symmetric proposal covariance matrix (or its inverse if `Σ⁻¹` was input).
                  The matrix is ensured to be positive semi-definite.
"""
modifyΣ(Σ, m, θ) = modifyΣ(Val(Symbol(m.model_name)), m, Σ, θ)
function modifyΣ(v::Val{:SW}, m, Σ, θ)
    _Σ = Matrix(Σ)
    # kill cross-correlations
    for (i, param) ∈ enumerate(m.parameters)
        if param ∈ (:ιp, :ιw, :ρ, :ρa, :ρb, :ρG, :ρI, :ρp, :ρw)
            if θ[i] > 0.95 || θ[i] < 0.05
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        elseif param ∈ (:h,:ρw, :ξw, :ξp)
            if θ[i] > 0.9 || θ[i] < 0.1
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        elseif param ∈ (:νl,)
            if θ[i] < 0.5
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        elseif param ∈ (:φ,)
            if θ[i] < 0.3
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        elseif param ∈ (:rπ1,)
            if θ[i] < 1.1
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        end
    end
    _Σ = Symmetric(_Σ)
    if !isposdef(_Σ)
        _Σ = make_psd(_Σ)
    end
    return _Σ
end
function modifyΣ(v::Val{:SWFF}, m, Σ, θ)
    _Σ = Matrix(Σ)
    # kill cross-correlations
    for (i, param) ∈ enumerate(m.parameters)
        if param ∈ (:ιp, :ιw, :ρ, :ρa, :ρb, :ρG, :ρI, :ρp, :ρw)
            if θ[i] > 0.95 || θ[i] < 0.05
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        elseif param ∈ (:h,:ρw, :ξw, :ξp)
            if θ[i] > 0.9 || θ[i] < 0.1
                _Σ[i,:] .= 0.0
                _Σ[:,i] .= 0.0
                _Σ[i,i] = Σ[i,i] 
            end
        end
    end
    _Σ = Symmetric(_Σ)
    if !isposdef(_Σ)
        _Σ = make_psd(_Σ)
    end
    return _Σ
end

"""
    get_measurement_matrix(m::MacroModelling.ℳ, 𝐒mat::Matrix, data::Workspace; core=true)

Construct the measurement matrix `M` (often denoted as `Λ` or `H` in state-space models)
that links the model's state variables to the observed data series.

The measurement equation is typically `Xₜ = M Sₜ + constant + errorₜ`. This function
determines the `M` matrix based on the model definition and the specified data.

It identifies which model variables correspond to the observed series listed in
`data.series` (or `data.core_series` if `core=true`). For each observed series,
it finds the corresponding state variable in the model's solution.
The `M` matrix will have dimensions `(number_of_observed_series, number_of_states)`.
An element `M[i, j]` will be 1 if observed series `i` corresponds to state variable `j`,
and 0 otherwise. This assumes a direct observation of states without scaling.

# Arguments
- `m::MacroModelling.ℳ`: The model object, containing solved state variables.
- `𝐒mat::Matrix`: The steady-state matrix obtained from the model solution (used to get state variable names via `𝐒mat.Variables`).
- `data::Workspace`: A Workspace containing data information, specifically `data.series`
                     or `data.core_series` which lists the names of the observed variables.
- `core=true`: If `true`, uses `data.core_series` to determine observed variables.
               If `false`, uses `data.series`.

# Returns
- `M::Matrix{Float64}`: The `(number_of_observed_series x number_of_states)` measurement matrix.
"""
function get_measurement_matrix(m, 𝐒, data; core = false)

    G, H, SS = get_sparse_solution(m, 𝐒)
    states = G.Variables
    nVars = data.nVars

    vars = names(data.X)
    if core == true 
        vars = names(data.X_core)
        nVars = length(vars)
    end

    Λ = zeros(nVars, size(G,2))
    
    for (i,var) in enumerate(vars)
        var_sym = Symbol(var)
        if var ∈ data.core_series
            j = findfirst(x -> x == Symbol(var), data.core_series_modelsyms)
            var_sym = data.core_series_modelsyms[j]
        end
        if var_sym ∈ G.Variables
            state_index = findfirst(x -> x == var_sym, G.Variables)
            Λ[i, state_index] = 1
            if var_sym == :R
                Λ[i, state_index] = 4
            end
        end
    end

    return Λ
end


"""
    get_measurement_const(m::MacroModelling.ℳ, θ::Vector, data::Workspace)

Calculate the constant term `A` in the measurement equation.

This function is intended to compute the constant vector `A` in a measurement equation
of the form `Xₜ = M Sₜ + A + errorₜ`. However, the current implementation always
returns a zero vector of the same length as `data.series`. This implies that the
measurement equation used in the context of this function assumes no constant term,
or that the constant is handled elsewhere (e.g., by de-meaning the data or incorporating
it into the steady state of the states `Sₜ`).

# Arguments
- `m::MacroModelling.ℳ`: The model object (currently unused by the function).
- `θ::Vector`: The vector of model parameter values (currently unused by the function).
- `data::Workspace`: A Workspace containing data information, specifically `data.series`
                     to determine the dimension of the constant vector.

# Returns
- `A::Vector{Float64}`: A zero vector of length equal to the number of series in `data.series`.
"""
get_measurement_const(m::MacroModelling.ℳ, θ::Vector{Float64}, data::Workspace) = get_measurement_const(Val(Symbol(m.model_name)), m, θ, data.core_series_modelsyms)
function get_measurement_const(v::Val{:SW}, m::MacroModelling.ℳ, θ::Vector{Float64},  syms::Vector{Symbol}) 
    TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    return SS_and_pars
end
function get_measurement_const(v::Val{:SWFF}, m::MacroModelling.ℳ, θ::Vector{Float64}, syms::Vector{Symbol}) 
    TT, SS_and_pars, 𝐒, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    return SS_and_pars
end

"""
    get_Q(m::MacroModelling.ℳ, θ::Vector)

Construct the covariance matrix `Q` of the structural shocks (state innovations) `ϵₜ`.

The state transition equation is typically `Sₜ = G Sₜ₋₁ + H ϵₜ`, where `ϵₜ ~ N(0, Q)`.
This function creates `Q` based on the model's shock definitions (`m.exo`) and their
standard deviations specified in the parameter vector `θ`.

It iterates through the exogenous shocks defined in `m.exo`. For each shock `shk`,
it looks for a corresponding parameter named `σ{shk_name_without_ϵ}` (e.g., if `shk` is `:ϵr`,
it looks for `:σr`). The square of this standard deviation parameter becomes the
diagonal element of `Q` for that shock. Off-diagonal elements are assumed to be zero,
meaning shocks are uncorrelated.

If a standard deviation parameter is not found for a shock, the corresponding diagonal
element of `Q` defaults to 1.0.

# Arguments
- `m::MacroModelling.ℳ`: The model object, containing exogenous shock definitions (`m.exo`)
                         and parameter names (`m.parameters`).
- `θ::Vector`: The vector of model parameter values.

# Returns
- `Q::Diagonal{Float64, Vector{Float64}}`: A diagonal matrix representing the covariance
                                          matrix of the structural shocks.
"""
function get_Q(v::Val{:SW}, m::MacroModelling.ℳ,  θ)
    return Diagonal(θ[22:29] .^ 2)
end
function get_Q(v::Val{:SWFF}, m::MacroModelling.ℳ,  θ)
    return Diagonal(θ[24:31] .^ 2)
end
get_Q(m::MacroModelling.ℳ, θ) = get_Q(Val(Symbol(m.model_name)), m, θ)