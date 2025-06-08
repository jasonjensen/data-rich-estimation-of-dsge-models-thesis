"""
    get_dist_R_Λ(X, ρ, S::KeyedArray, s₀=0.001, v₀=3, 𝚲₀=zeros(size(S, 1)), 𝐌₀=I(size(S, 1)), c=1.0; constraint=:none, state_index=0)

Calculate the conditional posterior distributions for R (measurement error variance) and Λ (observation matrix coefficients).

# Arguments
- `X`: Time series of observed data.
- `ρ`: Autocorrelation coefficient for detrending.
- `S::KeyedArray`: Time series of states.
- `s₀=0.001`: Prior scale for R.
- `v₀=3`: Prior degrees of freedom for R.
- `𝚲₀=zeros(size(S, 1))`: Prior mean for Λ.
- `𝐌₀=I(size(S, 1))`: Prior precision for Λ.
- `constraint=:none`: Type of constraint on Λ (:none, :one, :four, :any, :anyfour).
- `state_index=0`: Index of the state variable if a constraint is applied.

# Returns
- `R_dist`: InverseGamma distribution for R.
- `Λ_dist`: Normal or MvNormal distribution for Λ.
"""
function get_dist_R_Λ(X, ρ, S::KeyedArray, s₀=0.001, v₀=3, 𝚲₀=zeros(size(S, 1)), 𝐌₀=I(size(S, 1)); constraint=:none, state_index=0)
   
    Xᵡ =  X[2:end] - ρ*X[1:end-1] # discard first obs
    Sᵡ     = Matrix(S')[2:end, :] .- ρ .* Matrix(S')[1:end-1, :]
    
    if constraint ∈ (:one, :four, :any, :anyfour)
        Sᵡ = reshape(Sᵡ[:,state_index], (size(Sᵡ,1),1))
        𝐌₀ = reshape([1.0], (1,1))
    end
    
    nObs = size(Sᵡ,1)
    
    𝚲fit = nothing
    𝚲hat = zeros(size(S, 1))

    if constraint == :none
        𝚲fit = lm(Sᵡ, Xᵡ)
        𝚲hat = coef(𝚲fit)
    elseif constraint == :one
        𝚲hat = [1.0]
        𝚲₀ = [1]
    elseif constraint == :four
        𝚲hat = [4.0]
        𝚲₀ = [4]
    elseif constraint == :any
        𝚲fit = lm(Sᵡ, Xᵡ)
        𝚲hat = [coef(𝚲fit)[1]]
        𝚲₀ = [1]
    elseif constraint == :anyfour
        𝚲fit = lm(Sᵡ, Xᵡ)
        𝚲hat = [coef(𝚲fit)[1]]
        𝚲₀ = [4]
    end


    𝐌bar = 𝐌₀ + (Sᵡ' * Sᵡ)
    𝚲bar = inv(𝐌bar) * ((Sᵡ' * Sᵡ)*𝚲hat + 𝐌₀* 𝚲₀)
    XᵀX = Xᵡ' * Xᵡ
    
    if XᵀX isa Matrix{Float64}
        s₀ =  reshape([s₀], (1,1))
    end
    
    if constraint == :none
        sbar = s₀ + 0.5*(XᵀX - 𝚲bar' * (𝐌₀ + (Sᵡ' * Sᵡ)) * 𝚲bar + 𝚲₀' * 𝐌₀ * 𝚲₀)
        sbar = first(sbar)
    else
        sbar = first(s₀) + 0.5*((Xᵡ - Sᵡ * 𝚲hat)' * (Xᵡ - Sᵡ * 𝚲hat))
    end
    vbar = v₀ + (nObs - 1)/2

    # InverseGamma(α::shape, θ::Scale) ; shape controls tail heaviness, scale controls spread
    # The shape parameter is tied to degrees of freedom
    R_dist = Distributions.InverseGamma(vbar, sbar)
    
    Λ_dist = nothing
    if constraint ∈ (:one, :four, :any, :anyfour)
        Λ_dist = Distributions.Normal(𝚲bar[1], sqrt.(inv(𝐌bar[1])))
    else
        𝐌bar⁻¹ = isposdef(𝐌bar) ? inv(𝐌bar) : pinv(𝐌bar)
        if !isposdef(𝐌bar⁻¹)
            𝐌bar⁻¹ = make_psd(𝐌bar⁻¹)
        end
        Λ_dist = Distributions.MvNormal(𝚲bar, 𝐌bar⁻¹)
    end
    
    return R_dist, Λ_dist
end

"""
    get_dist_ψ(X, S, Λₖ, Rₖₖ, ρ₀, σ², c=1.0)

Calculate the conditional posterior distribution for ψ (autocorrelation of measurement error).

# Arguments
- `X`: Time series of observed data.
- `S`: Time series of states.
- `Λₖ`: Observation matrix coefficient.
- `Rₖₖ`: Measurement error variance.
- `ρ₀`: Prior mean for ψ.
- `σ²`: Prior variance for ψ.

# Returns
- `ψ_dist`: Truncated Normal distribution for ψ.
"""
function get_dist_ψ(X, S, Λₖ, Rₖₖ, ρ₀, σ²)
    
    eₖ = X - Matrix(S'[1:end-1,:]) * Λₖ
    
    𝛙hatₖₖ = inv(eₖ[1:end-1]' * eₖ[1:end-1]) * eₖ[1:end-1]' * eₖ[2:end]

    Vbar = inv(inv(Rₖₖ * inv(eₖ[1:end-1]' * eₖ[1:end-1])) +  inv(σ²))
    𝛙bar = Vbar * (inv(Rₖₖ * inv(eₖ[1:end-1]' * eₖ[1:end-1])) * 𝛙hatₖₖ + inv(σ²)*ρ₀)
    
    ψ_dist = truncated(Distributions.Normal(𝛙bar, sqrt(Vbar)), 0.0 + eps(Float64), 1.0-eps(Float64))

    return ψ_dist
end

"""
     draw_ψ(X, S, Λₖ, Rₖₖ, ρ₀, σ², c=1.0)

Sample ψ (autocorrelation of measurement error) from its conditional posterior distribution.

# Arguments
- `X`: Time series of observed data.
- `S`: Time series of states.
- `Λₖ`: Observation matrix coefficient.
- `Rₖₖ`: Measurement error variance.
- `ρ₀`: Prior mean for ψ.
- `σ²`: Prior variance for ψ.

# Returns
- `ψ_draw`: A draw for ψ.
- `loglikψ`: Log-likelihood of the draw.
"""
function draw_ψ(X, S, Λₖ, Rₖₖ, ρ₀, σ²)
    ψ_dist = get_dist_ψ(X, S, Λₖ, Rₖₖ, ρ₀, σ²)
    ψ_draw = rand(ψ_dist)
    ψ_prior = Distributions.Normal(0,1)
    loglikψ = logpdf(ψ_prior, ψ_draw)
    return ψ_draw, loglikψ
end

"""
    get_initial_Γ(S::KeyedArray, data::Workspace, Λ_constraints=Dict{String, Pair{Symbol,Symbol}}())

Helper function to call `get_initial_Γ` with data from a Workspace.

# Arguments
- `S::KeyedArray`: KeyedArray of states.
- `data::Workspace`: Workspace object containing data.
- `Λ_constraints=Dict{String, Pair{Symbol,Symbol}}()`: Dictionary of constraints for Λ.

# Returns
- Tuple: Contains initial estimates for R, Λ, 𝛙, and 𝛙_σ.
"""
function get_initial_Γ(S::KeyedArray, data::Workspace, Λ_constraints=Dict{String, Pair{Symbol,Symbol}}())
    return get_initial_Γ(data.X, S, data.nVars, data.nStates, data.core_series, data.core_series_modelsyms, Λ_constraints)
end

# TODO: clean up get_initial_Γ vs get_initial_Γ_old
"""
    get_initial_Γ(X::DataFrame, S::KeyedArray, nVars::Int64=size(X,2), nStates::Int64=size(S,1), core_series=Vector{String}(), core_series_modelsyms=Vector{Symbol}(), Λ_constraints=Dict{String, Pair{Symbol,Symbol}}())

Calculate initial estimates for R, Λ, 𝛙, and 𝛙_σ.

This function iterates through the variables in `X`, calculates initial estimates for the measurement error variance (R),
observation matrix (Λ), and autocorrelation of measurement error (𝛙 and 𝛙_σ) for each variable.
It uses `get_initial_Γ_old` for the main calculations and then refines R and Λ using `get_dist_R_Λ`.

# Arguments
- `X::DataFrame`: DataFrame of observed data.
- `S::KeyedArray`: KeyedArray of states.
- `nVars::Int64`: Number of variables.
- `nStates::Int64`: Number of states.
- `core_series=Vector{String}()`: Vector of core series names.
- `core_series_modelsyms=Vector{Symbol}()`: Vector of core series model symbols.
- `Λ_constraints=Dict{String, Pair{Symbol,Symbol}}()`: Dictionary of constraints for Λ.

# Returns
- `R`: Initial estimate for the measurement error variance matrix.
- `Λ`: Initial estimate for the observation matrix.
- `𝛙`: Initial estimate for the autocorrelation of measurement error matrix.
- `𝛙_σ`: Initial estimate for the standard deviation of 𝛙.
"""
function get_initial_Γ(X::DataFrame, S::KeyedArray, nVars::Int64=size(X,2), nStates::Int64=size(S,1), core_series=Vector{String}(), core_series_modelsyms=Vector{Symbol}(), Λ_constraints=Dict{String, Pair{Symbol,Symbol}}())
    R, Λ, 𝛙, 𝛙_σ = get_initial_Γ_old(X,S,nVars,nStates, core_series, core_series_modelsyms, Λ_constraints)

    # Set Psi to very low for core series
    for (i, key) in enumerate(names(X))
        # @show key
        state_var, constraint = first(S.Variable) => :none
        state_index = 1
        if haskey(Λ_constraints, key)
            state_var, constraint = Λ_constraints[key]
            state_index = findfirst(x -> x == state_var, S.Variable)
        end
        if constraint == :one || constraint == :four
            𝛙[i,i] = 𝛙[i,i] * 1e-6
        end

        
        R_dist, Λ_dist = get_dist_R_Λ(X[!, key], 𝛙[i,i], S; constraint=constraint, state_index=state_index)
        R[i,i] = Distributions.mean(R_dist)
        if constraint == :none
            Λ_dist_alt = Distributions.MvNormal(Λ_dist.μ, sqrt(R[i,i] .* Λ_dist.Σ)) 
            Λ[i,:] = Distributions.mean(Λ_dist_alt)
        else
            Λ_dist_alt = Distributions.Normal(Λ_dist.μ, sqrt(R[i,i]) .* Λ_dist.σ) 
            Λ[i,:] .= 0
            Λ[i,state_index] = Distributions.mean(Λ_dist_alt)
        end
    end
    return (R, Λ, 𝛙, 𝛙_σ)
end

"""
    get_initial_Γ_old(X::DataFrame, S::KeyedArray, nVars::Int64=size(X,2), nStates::Int64=size(S,1), core_series=Vector{String}(), core_series_modelsyms=Vector{Symbol}(), Λ_constraints=Dict{String, Pair{Symbol,Symbol}}())

Calculate initial estimates for R, Λ, 𝛙, and 𝛙_σ using ordinary least squares (OLS).

This function iterates through each variable in `X` and performs the following steps:
1. Extracts the relevant state variables based on constraints.
2. Fits an OLS model to estimate the relationship between the variable and states.
3. Calculates the autocorrelation (ρ) of the residuals from the initial fit.
4. Refits the OLS model using GLS-transformed data (correcting for autocorrelation).
5. Stores the estimated coefficients (Λ_μ), standard errors (Λ_σ), variance of residuals (R_Θ), autocorrelation (𝛙_μ), and standard error of autocorrelation (𝛙_σ).

# Arguments
- `X::DataFrame`: DataFrame of observed data.
- `S::KeyedArray`: KeyedArray of states.
- `nVars::Int64`: Number of variables.
- `nStates::Int64`: Number of states.
- `core_series=Vector{String}()`: Vector of core series names (unused in current implementation, but kept for compatibility).
- `core_series_modelsyms=Vector{Symbol}()`: Vector of core series model symbols (unused in current implementation, but kept for compatibility).
- `Λ_constraints=Dict{String, Pair{Symbol,Symbol}}()`: Dictionary of constraints for Λ. 

# Returns
- `R_Θ`: Initial estimate for the measurement error variance matrix.
- `Λ_μ`: Initial estimate for the observation matrix (mean coefficients).
- `𝛙_μ`: Initial estimate for the autocorrelation of measurement error matrix.
- `𝛙_σ`: Initial estimate for the standard deviation of 𝛙_μ.
"""
function get_initial_Γ_old(X::DataFrame, S::KeyedArray, nVars::Int64=size(X,2), nStates::Int64=size(S,1), core_series=Vector{String}(), core_series_modelsyms=Vector{Symbol}(), Λ_constraints=Dict{String, Pair{Symbol,Symbol}}())
    Λ_μ = zeros(nVars, nStates)
    Λ_σ = zeros(nVars, nStates)
    Λ_vcov = zeros(nVars, nStates, nStates)
    R_Θ = Diagonal(zeros(nVars,nVars))
    𝛙_μ= Diagonal(zeros(nVars, nVars))
    𝛙_σ= Diagonal(zeros(nVars, nVars))

    for (i, key) in enumerate(names(X))
        
        # Retrieve constraint
        state_var, constraint = first(S.Variable) => :none
        state_index = 1
        if haskey(Λ_constraints, key)
            state_var, constraint = Λ_constraints[key]
            state_index = findfirst(x -> x == state_var, S.Variable)
        end
        

        # extract regression matrices
        # obs! no intercept
        _S = Matrix(S')
        
        Xₖ = X[2:end, key] # discard first obs
        
       
        # adjust regression matrices
        if constraint ∈ (:any, :one, :four, :anyfour)
            _S = _S[:,state_index:state_index]
        end
        if size(_S,1) !== size(Xₖ,1)
            _S = _S[end-size(Xₖ,1)+1:end, :]
        end
        
        initial_fit = lm(_S, Xₖ)
        
        # autocorrelation of errors
        resid = residuals(initial_fit)
        if constraint == :one
            resid = Xₖ - _S
        elseif constraint == :four
            resid = Xₖ - (4 * _S)
        end
        ρ_fit = lm(reshape(resid[1:end-1], (length(resid)-1, 1)), resid[2:end])
        ρ = ρ_fit.pp.beta0[1]
        𝛙_μ[i,i] = ρ
        𝛙_σ[i,i] = stderror(ρ_fit)[1]
        
        
        # revised fit
        y_revised =  initial_fit.rr.y[2:end] - ρ*initial_fit.rr.y[1:end-1]
        # would need to add a const*(1-ρ) component if there is a constant.
        X_revised = initial_fit.pp.X[2:end, :] .- ρ .* initial_fit.pp.X[1:end-1, :]
        if constraint == :one || constraint == :four
            y_revised = _S[2:end] - ρ*_S[1:end-1]
            X_revised = Xₖ[2:end] - ρ*Xₖ[1:end-1]
            X_revised = reshape(X_revised, (length(X_revised), 1))
        end
        
        λ_fit = lm(X_revised, y_revised)
        R_Θ[i,i] = StatsBase.var(residuals(λ_fit))

        if constraint == :none
            Λ_μ[i,:] = coef(λ_fit)
            Λ_σ[i,:] = stderror(λ_fit)
            Λ_vcov[i, :, :] = vcov(λ_fit)
        elseif constraint == :one
            Λ_μ[i,:] = zeros(size(Λ_μ, 2))
            Λ_σ[i,:] = zeros(size(Λ_μ, 2))
            Λ_μ[i,state_index] = 1.0
        elseif constraint == :four
            Λ_μ[i,:] = zeros(size(Λ_μ, 2))
            Λ_σ[i,:] = zeros(size(Λ_μ, 2))
            Λ_μ[i,state_index] = 4.0
        elseif constraint == :any
            Λ_μ[i,:] = zeros(size(Λ_μ, 2))
            Λ_σ[i,:] = zeros(size(Λ_μ, 2))
            Λ_μ[i,state_index] = ρ_fit.pp.beta0[1]
            Λ_σ[i,state_index] = stderror(λ_fit)[1]
        elseif constraint == :anyfour
            Λ_μ[i,:] = zeros(size(Λ_μ, 2))
            Λ_σ[i,:] = zeros(size(Λ_μ, 2))
            Λ_μ[i,state_index] = ρ_fit.pp.beta0[1]
            Λ_σ[i,state_index] = stderror(λ_fit)[1]
        end
    end

    return R_Θ, Λ_μ, 𝛙_μ, 𝛙_σ
end



"""
    draw_R_Λ(X, 𝛙_μ, S, c;  constraint=:none, state_index=0, var="")

Draw R and Λ from their conditional posterior distributions.

# Arguments
- `X`: Time series of observed data.
- `𝛙_μ`: Mean of the autocorrelation of measurement error.
- `S`: Time series of states.
- `constraint=:none`: Type of constraint on Λ (:none, :one, :four, :any, :anyfour).
- `state_index=0`: Index of the state variable if a constraint is applied.
- `var=""`: Name of the variable (used for special rules for FedFunds and RINV).

# Returns
- `R_draw`: Draw of the measurement error variance.
- `Λ_draw`: Draw of the observation matrix coefficients.
- `0.0`: Placeholder for log-likelihood of 𝛙 (not calculated in this function).
- `0.0`: Placeholder for log-likelihood of 𝛙_σ (not calculated in this function).
"""
function draw_R_Λ(X, 𝛙_μ, S;  constraint=:none, state_index=0, var="")
   
    R_dist, Λ_dist = get_dist_R_Λ(X, 𝛙_μ, S, constraint=constraint, state_index=state_index)
    R_draw = 1.0 .* rand(R_dist)
    
    # special rules for FedFunds and RINV
    if var ∈ ("FedFunds",)
        while R_draw > 0.05
            vbar = R_dist.invd.α * 2
            sbar = R_dist.invd.θ * 2
            R_dist = Distributions.InverseGamma((vbar+1)/2, 2/sbar)
            R_draw = rand(R_dist)
        end
    end
    if var == "RINV"
        while R_draw > 2.0
            vbar = R_dist.invd.α * 2
            sbar = R_dist.invd.θ * 2
            R_dist = Distributions.InverseGamma((vbar+1)/2, 2/sbar)
            R_draw = rand(R_dist)
        end
    end
    
    # loglikR = logpdf(Distributions.InverseGamma(0.001,3), R_draw)
    
    Λ_draw = zeros(size(S,1))
    if constraint == :none
        Λ_dist_corrected = Distributions.MvNormal(Λ_dist.μ, R_draw .* Λ_dist.Σ)
        Λ_draw = rand(Λ_dist_corrected)
    else
        Λ_dist_corrected = Distributions.Normal(Λ_dist.μ, sqrt(R_draw) .* Λ_dist.σ) 
        Λ_draw[state_index] = rand(Λ_dist_corrected)
    end
    
    return R_draw, Λ_draw, 0.0, 0.0
end


"""
    draw_Γ(X::DataFrame, S, nVars, nStates, 𝛙_μ, 𝛙_σ, core_series=Vector{String}(), core_series_modelsyms=Vector{Symbol}(), Λ_constraints=Dict{String, Pair{Symbol,Symbol}}(), c=1.0, measurement_error=true, verbose=false)

Draw Γ (Λ, 𝛙, R) from their conditional posterior distributions.

This function iterates through each variable in `X` and performs the following:
1. Draws R and Λ using `draw_R_Λ`.
2. Applies constraints to Λ if specified.
3. If `measurement_error` is true, samples 𝛙 using `draw_ψ`.
4. Calculates the total log-likelihood for Γ.

# Arguments
- `X::DataFrame`: DataFrame of observed data.
- `S`: Time series of states.
- `nVars`: Number of variables.
- `nStates`: Number of states.
- `𝛙_μ`: Mean of the prior for 𝛙.
- `𝛙_σ`: Standard deviation of the prior for 𝛙.
- `core_series=Vector{String}()`: Vector of core series names.
- `core_series_modelsyms=Vector{Symbol}()`: Vector of core series model symbols.
- `Λ_constraints=Dict{String, Pair{Symbol,Symbol}}()`: Dictionary of constraints for Λ.
- `measurement_error=true`: Boolean indicating whether to model measurement error.
- `verbose=false`: Boolean indicating whether to print verbose output.
"""
function draw_Γ(X::DataFrame, S, nVars, nStates, 𝛙_μ, 𝛙_σ, core_series=Vector{String}(), core_series_modelsyms=Vector{Symbol}(), Λ_constraints=Dict{String, Pair{Symbol,Symbol}}(), measurement_error=true, verbose=false)
    local Λ = zeros(nVars, nStates)
    local R = Diagonal(zeros(nVars, nVars))
    local 𝛙= Diagonal(zeros(nVars, nVars))
    local loglikΛ = zeros(nVars)
    local loglikR = zeros(nVars)
    local loglikψ = zeros(nVars)

    for (i, key) in enumerate(names(X))
        constraint = :none
        state_index=0
        if haskey(Λ_constraints, key)
            state_var, constraint = Λ_constraints[key]
            state_index = findfirst(x -> x == state_var, S.Variable)
        end

        R[i,i], Λ[i,:], loglikR[i], loglikΛ[i]  = draw_R_Λ(X[!, key], 𝛙_μ[i,i], S; constraint=constraint, state_index=state_index, var=key)
        
        if haskey(Λ_constraints, key)
            state_var, constraint = Λ_constraints[key]
            state_index = findfirst(x -> x == state_var, S.Variable)
            old_val =  Λ[i,state_index]
            Λ[i,:] = zeros(size(Λ, 2))
            if constraint == :one
                Λ[i,state_index] = 1.0
                loglikR[i] = 0
                loglikΛ[i] = 0
            elseif constraint == :four
                Λ[i,state_index] = 4.0
                loglikR[i] = 0
                loglikΛ[i] = 0
            elseif constraint ∈ (:any, :anyfour)
                Λ[i,state_index] = old_val
            end
        end

        if measurement_error
            𝛙[i,i], loglikψ[i] = draw_ψ(X[2:end, key], S, Λ[i,:], R[i,i], 𝛙_μ[i,i], 𝛙_σ[i,i])
        end
    end
    
    loglikΓ = sum(loglikΛ) + sum(loglikR) + sum(loglikψ)
   
    if verbose
        println("==================")
        println("Λ: $(sum(loglikΛ))")
        display(Λ)
        println("R: $(sum(loglikR))")
        display(R)
        println("ψ: $(sum(loglikψ))")
        display(𝛙)
        println("==================")
    end
    
    if !measurement_error
        R = Diagonal(zeros(nVars, nVars))
        loglikΓ = 0 # simplifying assumption that Λ is fixed when measurement_error==false
    end


    return Λ, 𝛙, R, loglikΓ
end

"""
    draw_Γ(S::KeyedArray, data::Workspace, 𝛙_μ::Diagonal, 𝛙_σ::Diagonal, Λ_constraints=Dict{String, Pair{Symbol,Symbol}}(), measurement_error=true, verbose=false)

Draw Γ (Λ, 𝛙, R) from their conditional posterior distributions using data from a Workspace.

# Arguments
- `S`: Time series of states.
- `data`: Workspace containing X and other information.
- `𝛙_μ`: Mean of the prior for 𝛙.
- `𝛙_σ`: Standard deviation of the prior for 𝛙.
- `Λ_constraints=Dict{String, Pair{Symbol,Symbol}}()`: Dictionary of constraints for Λ.
- `measurement_error=true`: Boolean indicating whether to model measurement error.
- `verbose=false`: Boolean indicating whether to print verbose output.

# Returns
- `R`: Draw of the measurement error variance matrix.
- `Λ`: Draw of the observation matrix.
- `𝛙`: Draw of the autocorrelation of measurement error matrix.
- `𝛙_σ`: Draw of the standard deviation of 𝛙.
- `loglik`: Total log-likelihood for Γ.
"""
function draw_Γ(S::KeyedArray, data::Workspace, 𝛙_μ::Diagonal, 𝛙_σ::Diagonal, Λ_constraints=Dict{String, Pair{Symbol,Symbol}}(), measurement_error=true, verbose=false)
    return draw_Γ(data.X, S, data.nVars, data.nStates, 𝛙_μ, 𝛙_σ, data.core_series, data.core_series_modelsyms, Λ_constraints, measurement_error, verbose)
end


"""
    perturb_θ(m, θprev, Σ⁻¹, c; fixed_parameters=Vector{Symbol}())

Perturb the previous parameter vector `θprev` and evaluate the prior likelihood.

# Arguments
- `m`: The model object.
- `θprev`: The previous parameter vector.
- `Σ⁻¹`: Inverse of the covariance matrix for the proposal distribution.
- `c`: Scaling factor for the proposal distribution.
- `fixed_parameters=Vector{Symbol}()`: A vector of symbols indicating parameters that should not be perturbed.

# Returns
- `θ`: The new (potentially perturbed) parameter vector.
- `total_loglik_new`: The log-likelihood of the new parameters under their priors.
- `𝐒mat`: The new steady-state solution if the model solves with the new parameters.
- `true` if the model solves with the new parameters, `false` otherwise.
"""
function perturb_θ(m, θprev, Σ⁻¹, c; fixed_parameters=Vector{Symbol}())
    dist = MvNormal(θprev, c * Σ⁻¹)
    perturbations = rand(dist)

    θ = copy(θprev)
    
    total_loglik_new = 0.0
    for (i, param, new_value, prev_value) in zip(1:length(m.parameters), m.parameters, perturbations, θprev)
        if param ∉ fixed_parameters
            θ[i] = clamp_θ(Val(Symbol(m.model_name)), Val(param), new_value)
            if θ[i] != new_value 
                return (θprev, 0, nothing, false)
            end
        end
        
        # compute loglik
        if param ∉ fixed_parameters
            loglik = logpdf(
                prior_θ(Val(Symbol(m.model_name)), Val(param)), 
                θ[i]
            )
            total_loglik_new += loglik
        end
    end

    # make sure it works
    sol_guess = copy(m.solution.perturbation.qme_solution)
    TT, SS_and_pars, 𝐒mat, state, solved = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ, m)
    if !solved
        m.solution.perturbation.qme_solution = copy(sol_guess)
        return θprev, 0, nothing, false
    end

    # total_loglik_new = sum(logpdf(dist, perturbations))
    return θ, total_loglik_new, 𝐒mat, true
end



