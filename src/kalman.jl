"""
    generate_states(m, θ, 𝐒, data, Λ, 𝛙, R; burnin=0, backwards_pass=true, speedup = false, P₀=nothing)

Generate states using a Kalman filter and smoother, with an option for a speedup.

This function implements a Kalman filter and an optional backward smoothing pass (Carter-Kohn style)
to generate states for a given model `m`, parameters `θ`, steady-state `𝐒`, data,
and measurement equation matrices `Λ`, `𝛙`, `R`.

It includes an optional speedup based on Jungbacker and Koopman (2015).

# Arguments
- `m`: The model object.
- `θ`: The parameter vector.
- `𝐒`: The steady-state solution.
- `data`: A Workspace containing the observed data (`X`, `nObs`, `nStates`, `nVars`).
- `Λ`: The observation matrix (loading factors).
- `𝛙`: The autocorrelation matrix for measurement errors.
- `R`: The covariance matrix for measurement errors.
- `burnin=0`: Number of initial periods to discard from the likelihood calculation.
- `backwards_pass=true`: Whether to perform the backward smoothing pass. If `false`, only filtered states are returned.
- `speedup=false`: Whether to use the Jungbacker and Koopman (2015) speedup for the Kalman filter.
- `P₀=nothing`: Initial covariance matrix for the states. If `nothing`, it's initialized to zeros or computed from the model.

# Returns
- If `backwards_pass == true`:
    - `Sₜ`: A KeyedArray of smoothed states.
    - `log_likelihood`: The log-likelihood of the data given the parameters.
    - `sum_predicted_states`: Sum of predicted states (used for debugging/diagnostics).
- If `backwards_pass == false`:
    - `Sₜ`: A KeyedArray of filtered states.
    - `log_likelihood`: The log-likelihood of the data given the parameters.
"""
function generate_states(m, θ, 𝐒, data, Λ, 𝛙, R; burnin=0, backwards_pass=true, speedup = false, P₀=nothing)
   
   @timeit to "preamble" begin
        # countables
        nObs = data.nObs -1
        nStates = data.nStates
        nVars = data.nVars
        
        # nudges to help with invertibility
        nudge   = Symmetric(1e-12*I(nStates*2))
        nudgeᵡ  = 1e-7*I(nStates)
        nudgeᶠ  = 1e-9*I(nStates*2)
        nudgeᴿ  = 1e-9*I(nVars)
        if speedup
            nudgeᴿ = 1e-9*I(nStates*2)
        end

        # measurement matrix
        M = zeros(nStates*2, nStates*2)
        M[1:nStates,1:nStates] = Matrix(I(nStates))

        # get the transition matrix
        zeros_nStates = zeros(nStates,nStates)
        G_keyed, H_keyed, SS = get_sparse_solution(m, 𝐒)

        # raw transition matrices
        G = Matrix(G_keyed)
        H = Matrix(H_keyed)
        
        # padded transition matrices
        G̃ᵡ = G
        G̃ = [G̃ᵡ          zeros_nStates
            I(nStates)  zeros_nStates]
        G̃ᵀ = G̃'  
        H̃ = vcat(Matrix(H), zeros(size(H)))

        # shock matrix sandwich
        Q = I(size(H̃,2)) # Q is just the Identity matrix. H depends on θ
        H̃QH̃ᵀ = Symmetric(H̃*Q*H̃')
        Σᵤ = H̃QH̃ᵀ

        # corrected Λ
        Λ̃  = [Λ -𝛙*Λ]
        Λᵀ = (Λ̃ )'

        # corrected data
        X = Matrix(data.X)'
        X̃ = X[:,2:end] - 𝛙*X[:,1:end-1]

        # Jungbacker and Koopman (2008) speedup
        C = I(nStates)
        Aᴸ = I(nVars)
        Rᴸ = copy(R)
        
        if speedup
            R⁻¹ = pinv(R)
            # Jungbacker and Koopman (2008) speedup
            @timeit to "C" C = inv(Λᵀ * R⁻¹ * Λ̃  + nudge*1000)
            @timeit to "Aᴸ" Aᴸ = C * Λᵀ * R⁻¹
            # @assert C' ≈ C
        end
        
        X̃ᴸ = Aᴸ * X̃
        
        AᴸΛ = Aᴸ * Λ̃ 
        AᴸΛᵀ = AᴸΛ'
        Rᴸ = Symmetric(Aᴸ * R * Aᴸ') + nudgeᴿ
        Qᵡ = (M*Σᵤ*M')[1:nStates, 1:nStates]

        
        nudgeᶠ = 1e-9*I(size(Rᴸ,1))
        

        # effective number of variables
        nVarsᴸ = size(X̃ᴸ,2)

        # loglike baseline
        loglike_baseline = 0
        if speedup
            vSum = 0
            # Zˣ = Λ̃  * C
            # Λᵥ = Zˣ * inv(Zˣ' * R⁻¹ * Zˣ) * Zˣ' *R⁻¹
            Λᵥ = Λ̃  * inv(Λᵀ * R⁻¹ * Λ̃ ) * Λᵀ *R⁻¹
            for i = 1:nObs
                v = view(X̃,:,i) - Λᵥ*view(X̃,:,i)
                vAdd = v' * R⁻¹ * v
                vSum = vSum + vAdd
            end
            loglike_baseline = -0.5*(nVars-(nStates*2)*nObs)*log(2*pi) - 0.5*nObs*log(det(R)/det(Rᴸ)) -0.5*vSum
        end
        
        ## data structures
        @timeit to "prealloc" begin
            # prediction step
            S̃ₜ₊₁₎ₜ = repeat([zeros(nStates*2)], nObs) # Predicted states
            P̃ₜ₊₁₎ₜ = repeat([Symmetric(zeros( nStates*2, nStates*2))], nObs) 
            ηₜ₊₁₎ₜ =  repeat([zeros(nStates*2)], nObs) # prediction errors
            # fₜ₊₁₎ₜ =  repeat([Symmetric(zeros( nStates*2, nStates*2))], nObs) # prediction error variance
            
            # update step
            S̃ₜ₎ₜ = repeat([zeros(nStates*2)], nObs) # Smoothed states (filtered state meands)
            S̃ₜ₎ₜ = deepcopy(S̃ₜ₎ₜ)
            P̃ₜ₎ₜ = repeat([Symmetric(zeros( nStates*2, nStates*2))], nObs) # Smoothed state covariances

            # Kalman gain
            Kₜ = repeat([zeros( size(AᴸΛᵀ))], nObs)
          
            # log likelihood
            LX̃ᴸ = zeros(nObs)
            
            # Initial conditions
            S̃₀ = repeat(SS, 2)
            if P₀ === nothing
                P̃₀ = zeros(size(G̃))
            else
                P̃₀ = P₀
            end
          
            # workspaces
            K_workspace = Matrix{Float64}(undef, size(AᴸΛᵀ))
        end
    end

    # Forward pass
    @timeit to "forward pass"  @inbounds @simd for t in 1:nObs-1

        # Make predictions
        @timeit to "predict" begin
            if t == 1
                S̃ₜ₊₁₎ₜ[t] = G̃ * S̃₀
                P̃ₜ₊₁₎ₜ[t] = Symmetric(G̃*P̃₀*G̃ᵀ) + H̃QH̃ᵀ
            else
                S̃ₜ₊₁₎ₜ[t] = G̃ * S̃ₜ₎ₜ[t]
                P̃ₜ₊₁₎ₜ[t] = Symmetric(G̃ * P̃ₜ₎ₜ[t] * G̃ᵀ) + H̃QH̃ᵀ #+ nudge    
            end
            ηₜ₊₁₎ₜ[t] = view(X̃ᴸ,:,t+1) -  AᴸΛ * S̃ₜ₊₁₎ₜ[t] # prediction error
            F =  lu(AᴸΛ * P̃ₜ₊₁₎ₜ[t] * AᴸΛᵀ + Rᴸ + nudgeᶠ) # variance of prediction error
        end
        
        # Calculate Kalman gain
        @timeit to "inversion" begin
            invf = inv(F)
            mul!(K_workspace, P̃ₜ₊₁₎ₜ[t], AᴸΛᵀ)
            Kₜ[t] = K_workspace * invf
            # update_kalman_gain!(Kₜ[t], P̃ₜ₊₁₎ₜ[t], AᴸΛᵀ, fₜ₊₁₎ₜ[t], K_workspace)
        end
        
        # Update
        @timeit to "update" begin
            if t < nObs
               S̃ₜ₎ₜ[t+1] = S̃ₜ₊₁₎ₜ[t] + Kₜ[t] * ηₜ₊₁₎ₜ[t]
               P̃ₜ₎ₜ[t+1] = Symmetric(P̃ₜ₊₁₎ₜ[t] - Kₜ[t] * AᴸΛ * P̃ₜ₊₁₎ₜ[t] + nudge)
            end
            
        end

        #calculate loglikelihood
        @timeit to "loglik" begin 
            detF = det(F)
            if detF > 0.0
                LX̃ᴸ[t] = -(nVarsᴸ/2)*log(2*pi) -0.5*log(detF) - 0.5*((ηₜ₊₁₎ₜ[t])'*(invf*ηₜ₊₁₎ₜ[t]));
            else
                LX̃ᴸ[t] = -Inf;
            end
        end
    end

    if backwards_pass == false
        # return states and baseline loglikelihood
        @timeit to "hcat" Sₜ = KeyedArray(hcat(S̃ₜ₎ₜ...)[1:nStates,:], Variable = G_keyed.Variables, time=1:nObs)
        return (Sₜ, sum(LX̃ᴸ[burnin+1:end]) + loglike_baseline)
    end
    
    # Backward pass
    S̃ₜ = repeat([zeros(nStates)], nObs)
    @timeit to "backward pass" @inbounds @simd for t in nObs:-1:1
        if t == nObs
            tempP̃ₜ₎ₜ = Matrix(P̃ₜ₎ₜ[t][1:nStates,1:nStates])
            if isposdef(tempP̃ₜ₎ₜ) == false
                tempP̃ₜ₎ₜ = make_psd(tempP̃ₜ₎ₜ)
            end
            S̃ₜ[t] = rand(MvNormal(S̃ₜ₎ₜ[t][1:nStates], tempP̃ₜ₎ₜ))
        elseif t>0
            S̃ₜ₎ₜᵡ = @view S̃ₜ₎ₜ[t][1:nStates]
            P̃ₜ₎ₜᵡ = @view P̃ₜ₎ₜ[t][1:nStates,1:nStates]
        
            _partial = P̃ₜ₎ₜᵡ * G̃ᵡ' * inv(Symmetric(G̃ᵡ * P̃ₜ₎ₜᵡ * G̃ᵡ') + Qᵡ + nudgeᵡ) 
            
            S̃ₜ₎ₜₛ = S̃ₜ₎ₜᵡ +  _partial * (S̃ₜ[t+1] - G̃ᵡ*S̃ₜ₎ₜᵡ)
            P̃ₜ₎ₜₛ = Symmetric(P̃ₜ₎ₜᵡ - _partial * G̃ᵡ * P̃ₜ₎ₜᵡ + nudgeᵡ)
            if isposdef(P̃ₜ₎ₜₛ) == false
                P̃ₜ₎ₜₛ = make_psd(P̃ₜ₎ₜₛ)
            end

            # These states are being drawn conditional on error variances informed by Γ
            @timeit to "MvNormal" S̃ₜ[t] = rand(MvNormal(S̃ₜ₎ₜₛ, P̃ₜ₎ₜₛ))
        end
    end

    # return states
    Sₜ = copy(S̃ₜ)
    insert!(Sₜ, 1, zeros(nStates))
    Sₜ[1] = Sₜ[2]

    @timeit to "hcat" Sₜ = KeyedArray(hcat(Sₜ...), Variable = G_keyed.Variables, time=1:nObs+1)
    
    return (Sₜ, sum(LX̃ᴸ) + loglike_baseline, sum(sum(S̃ₜ₊₁₎ₜ)))

end

"""
    update_kalman_gain!(K::SubArray{Float64}, P::SubArray{Float64}, Λ::Adjoint{Float64}, f::Matrix{Float64}, _workspace::Matrix{Float64})
    update_kalman_gain!(K::Matrix{Float64}, P::Symmetric{Float64}, Λ::Adjoint{Float64}, f::Symmetric{Float64}, _workspace::Matrix{Float64})

Update the Kalman gain `K` in-place.

Calculates `K = P * Λ' * inv(f)` using a workspace to avoid allocations.

# Arguments
- `K`: Kalman gain matrix (output).
- `P`: Covariance matrix of the predicted state.
- `Λ`: Transpose of the observation matrix.
- `f`: Covariance matrix of the prediction error.
- `_workspace`: Pre-allocated workspace matrix.
"""
@views function update_kalman_gain!(K::SubArray{Float64}, P::SubArray{Float64}, Λ::Adjoint{Float64}, f::Matrix{Float64}, _workspace::Matrix{Float64})
    mul!(_workspace, P, Λ)
    K[:,:] = _workspace / f
end
@views function update_kalman_gain!(K::Matrix{Float64}, P::Symmetric{Float64}, Λ::Adjoint{Float64}, f::Symmetric{Float64}, _workspace::Matrix{Float64})
    mul!(_workspace, P, Λ)
    K[:,:] = _workspace / f
end


"""
    kalman(m, 𝐒, data, Λ, R=I(data.nVars), A=zeros(data.nVars), Q = I(data.nVars); debug=false, core=false, backward=false, burnin=0)

Perform Kalman filtering and smoothing.

Implements the standard Kalman filter and an optional backward smoothing pass.
The state-space model is defined as:
State transition: `Sₜ = G*Sₜ₋₁ + H*εₜ` where `εₜ ~ N(0,Q)`
Measurement: `X = Λ*Sₜ + A + νₜ` where `νₜ ~ N(0,R)` (Note: `A` is an intercept/offset term in the measurement equation).

# Arguments
- `m`: The model object.
- `𝐒`: The steady-state solution.
- `data`: A Workspace containing the observed data (`X`, `nObs`, `nStates`, `nVars`).
- `Λ`: The observation matrix.
- `R=I(data.nVars)`: Covariance matrix for measurement errors. Defaults to identity.
- `A=zeros(data.nVars)`: Intercept term in the measurement equation. Defaults to zeros.
- `Q=I(data.nVars)`: Covariance matrix for state shocks. Defaults to identity (Note: the docstring for `H` in the model implies `Q` might be implicitly handled by `H`).
- `debug=false`: If true, prints debugging information.
- `core=false`: If true, uses `data.X_core` instead of `data.X`.
- `backward=false`: If true, performs a backward smoothing pass.
- `burnin=0`: Number of initial periods to discard from the likelihood calculation.

# Returns
- `Sₜ`: A KeyedArray of filtered (if `backward=false`) or smoothed (if `backward=true`) states.
- `log_likelihood`: The log-likelihood of the data given the parameters and model.
"""
function kalman(m, 𝐒, data, Λ, R=I(data.nVars), A=zeros(data.nVars), Q = I(data.nVars); debug=false, core=false, backward=false, burnin=0)
    """
    state transition equations:
    Sₜ = G*Sₜ₋₁ + H*εₜ      where εₜ ~ N(0,I)
    measurement error equation:
    X = Λ*Sₜ + νₜ           where νₜ ~ N(0,R)
    """
    G_keyed, H_keyed, SS = get_sparse_solution(m, 𝐒)


    nObs = data.nObs
    nStates = data.nStates
    nVars = data.nVars
    
    G = Matrix(G_keyed)
    Gᵀ = G'
    H = Matrix(H_keyed)
    X = Matrix(data.X)'
    if core == true
        X = Matrix(data.X_core)'
        nVars = size(X,1)
    end

    H̃QH̃ᵀ = H*Q*H'  
    Λᵀ = Λ'

    S̃ₜ₊₁₎ₜ = zeros(nObs, nStates) # Predicted states
    P̃ₜ₊₁₎ₜ = zeros(nObs, nStates, nStates) # Predicted state covariances
    ηₜ₊₁₎ₜ = zeros(nObs, nVars) # prediction errors
    fₜ₊₁₎ₜ =  zeros(nObs, nVars, nVars) #something?
    
    S̃ₜ₎ₜ = zeros(nObs, nStates) # Smoothed states (filtered state meands)
    P̃ₜ₎ₜ = zeros(nObs, nStates, nStates) # Smoothed state covariances  
    P̃₀ = zeros(size(G))
    Kₜ = zeros(nObs, nStates, nVars)
    LX = zeros(nObs)
    S̃₀ = SS
    
    nudgeᶠ = 1e-12*I(nVars)
    nudgeᴾ = 1e-12*I(nStates)

    # workspaces:
    K_workspace = Matrix{Float64}(undef, nStates, nVars)

    # Forward pass
    @timeit to "forward pass" for t in 1:nObs-1
        
        # Prepare views
        ᵥS̃ₜ₎ₜ = @view S̃ₜ₎ₜ[t,:]
        ᵥP̃ₜ₎ₜ = @view P̃ₜ₎ₜ[t,:,:]
        ᵥS̃ₜ₊₁₎ₜ = @view S̃ₜ₊₁₎ₜ[t,:]
        ᵥP̃ₜ₊₁₎ₜ = @view P̃ₜ₊₁₎ₜ[t,:,:]
        ᵥηₜ₊₁₎ₜ = @view ηₜ₊₁₎ₜ[t,:]
        ᵥfₜ₊₁₎ₜ = @view fₜ₊₁₎ₜ[t,:,:]
        ᵥKₜ = @view Kₜ[t,:,:]

        # Make predictions
        if t == 1
            # update
            ᵥS̃ₜ₊₁₎ₜ .= G * S̃₀
            ᵥP̃ₜ₊₁₎ₜ .= G*P̃₀*Gᵀ + H̃QH̃ᵀ
        else
            # update
            ᵥS̃ₜ₊₁₎ₜ .= G * ᵥS̃ₜ₎ₜ
            ᵥP̃ₜ₊₁₎ₜ .= G * ᵥP̃ₜ₎ₜ * Gᵀ + H̃QH̃ᵀ    
            
        end

        #errors
        ᵥηₜ₊₁₎ₜ .= view(X,:,t+1) - (Λ*A + Λ * ᵥS̃ₜ₊₁₎ₜ)
        ᵥfₜ₊₁₎ₜ .= Λ * ᵥP̃ₜ₊₁₎ₜ * Λᵀ + R + nudgeᶠ      

        # Calculate Kalman gain
        if !isposdef(ᵥfₜ₊₁₎ₜ)
            ᵥfₜ₊₁₎ₜ .= make_psd(ᵥfₜ₊₁₎ₜ)
        end
        @timeit to "inversion" update_kalman_gain!(ᵥKₜ, ᵥP̃ₜ₊₁₎ₜ, Λᵀ, fₜ₊₁₎ₜ[t,:,:], K_workspace)
        
        LX[t] = -0.5*nVars*log(2*pi) - 0.5*log(det(ᵥfₜ₊₁₎ₜ)) - 0.5 * ᵥηₜ₊₁₎ₜ' * inv(copy(ᵥfₜ₊₁₎ₜ)) * ᵥηₜ₊₁₎ₜ

        # Update
        if t < nObs 
            ₙS̃ₜ₎ₜ = @view S̃ₜ₎ₜ[t+1,:]
            ₙP̃ₜ₎ₜ = @view P̃ₜ₎ₜ[t+1,:,:]

            ₙS̃ₜ₎ₜ .= ᵥS̃ₜ₊₁₎ₜ +  ᵥKₜ * ᵥηₜ₊₁₎ₜ
            ₙP̃ₜ₎ₜ .=  Symmetric(ᵥP̃ₜ₊₁₎ₜ - ᵥKₜ * Λ * ᵥP̃ₜ₊₁₎ₜ + nudgeᴾ) 
        end
    end

    Sₜ = KeyedArray(S̃ₜ₎ₜ', Variable = G_keyed.Variables, time=1:nObs)

    if backward
        Sₜ = zeros(nObs, nStates)
        LX = zeros(nObs)
        loglikS = 0
        @timeit to "backward pass" for t in nObs:-1:1
            ᵥS̃ₜ = @view Sₜ[t,:]
            
            if t == nObs
                if isposdef(P̃ₜ₎ₜ[t,1:nStates,1:nStates]) == false
                    P̃ₜ₎ₜ[t,1:nStates,1:nStates] = make_psd(P̃ₜ₎ₜ[t,1:nStates,1:nStates])
                end
                ᵥS̃ₜ .= rand(MvNormal(S̃ₜ₎ₜ[t,1:nStates], P̃ₜ₎ₜ[t,1:nStates,1:nStates]))
            elseif t>2
                ₙS̃ₜ =  @view Sₜ[t+1, 1:nStates]
                # these means are huge near the start
                S̃ₜ₎ₜᵡ = @view S̃ₜ₎ₜ[t, 1:nStates]
                P̃ₜ₎ₜᵡ = @view P̃ₜ₎ₜ[t,1:nStates,1:nStates]
                if isposdef(P̃ₜ₎ₜᵡ) == false
                    P̃ₜ₎ₜᵡ = make_psd(P̃ₜ₎ₜᵡ) + nudgeᴾ
                end
                # TODO: might need to be HQH (non-transpose)
                _partial = P̃ₜ₎ₜᵡ * G' * inv(G * P̃ₜ₎ₜᵡ * G' + H̃QH̃ᵀ + nudgeᴾ) 
                
                # divide by norm here to prevent trends/explosive growth.
                S̃ₜ₎ₜₛ = S̃ₜ₎ₜᵡ +  (_partial/norm(_partial)) * (ₙS̃ₜ - G*S̃ₜ₎ₜᵡ)
                P̃ₜ₎ₜₛ = Symmetric(P̃ₜ₎ₜᵡ - _partial * G * P̃ₜ₎ₜᵡ + nudgeᴾ)
                if isposdef(P̃ₜ₎ₜₛ) == false
                    P̃ₜ₎ₜₛ = make_psd(P̃ₜ₎ₜₛ)
                end
    
                # These states are being drawn conditional on error variances informed by Γ
                # dist=
                ᵥS̃ₜ .= rand(MvNormal(S̃ₜ₎ₜₛ, P̃ₜ₎ₜₛ))
                LX[t] = logpdf(MvNormal(S̃ₜ₎ₜₛ, P̃ₜ₎ₜₛ), ᵥS̃ₜ)
            end
        end
        Sₜ = KeyedArray(Sₜ', Variable = G_keyed.Variables, time=1:nObs)
    end

    return Sₜ, sum(LX[burnin+1:end])

end

"""
    get_sparse_solution(m::MacroModelling.ℳ, 𝐒)

Get the sparse state-space solution matrices (G, H) and steady-state (SS) from the model.

This function retrieves the solution of the model `m`, sets the state variables to the provided
steady-state `𝐒`, and extracts the state transition matrix `G`, shock impact matrix `H`,
and the steady-state vector `SS`. It also zeros out small values in the solution matrices.

# Arguments
- `m::MacroModelling.ℳ`: The model object.
- `𝐒`: The steady-state vector to impose on the model's solution.

# Returns
- `G`: KeyedArray representing the state transition matrix (mapping `VariablesPast` to `Variables`).
- `H`: KeyedArray representing the shock impact matrix.
- `SS`: KeyedArray representing the steady-state vector.
"""
function get_sparse_solution(m::MacroModelling.ℳ, 𝐒)
    GH = get_solution(m) #TODO: I only need this for the scaffolding, maybe store somewhere and pass along
    GH[2:end, :] = 𝐒'
    mask = abs.(GH) .< 1e-14
    GH[mask] .= 0.0
    nStates = length(GH.Variables)
    
    # make it square
    # prev_states will be the columns
    G = zeros(nStates, nStates)
    for (i,var) ∈ enumerate(GH.Variables)
        row_name = Symbol("$(var)₍₋₁₎")
        if row_name ∈ GH.Steady_state__States__Shocks
            G[:,i] =  GH(row_name,:)'
        end
    end
    
    shock_vars = Vector{Symbol}()
    for var in GH.Steady_state__States__Shocks
        if contains(String(var), "₍ₓ₎")
            push!(shock_vars, var)
        end
    end

    G = KeyedArray(G, Variables=GH.Variables, VariablesPast=GH.Variables)

    H = GH(shock_vars, :)'
    
    return G, H, GH(:Steady_state,:)
end

"""
    get_P₀(m, 𝐒)

Compute the initial covariance matrix `P̃₀` for the Kalman filter based on the model's
state-space representation, assuming an augmented state vector.

This function calculates the unconditional covariance matrix of an augmented state
vector (typically `[Sₜ', Sₜ₋₁']'`, where `Sₜ` is the original state vector).
It does this by solving the discrete Lyapunov equation:
`P̃₀ = G̃ P̃₀ G̃' + H̃ Q H̃'`
for `P̃₀`, where `G̃` is the transition matrix of the augmented state vector and
`H̃ Q H̃'` is its innovation covariance. The solution is found by vectorizing the
equation: `vec(P̃₀) = inv(I - G̃ ⊗ G̃) * vec(H̃ Q H̃')`.

This `P̃₀` represents the unconditional variance of the augmented state process and is
often used as the initial state covariance matrix when starting the Kalman filter,
assuming the process is stationary and no specific prior for the initial state is available.

# Arguments
- `m::MacroModelling.ℳ`: The model object, used to derive the state-space matrices.
- `𝐒`: The steady-state solution of the model, used to obtain `G_keyed` and `H_keyed` via `get_sparse_solution`.

# Returns
- `P̃₀::Matrix{Float64}`: The initial covariance matrix for the augmented state vector.
"""
function get_P₀(m, 𝐒)
    G_keyed, H_keyed, SS = get_sparse_solution(m, 𝐒)
    nStates = length(SS)
    zeros_nStates = zeros(nStates,nStates)
    G̃ = [Matrix(G_keyed)    zeros_nStates
        I(nStates)  zeros_nStates]
    H̃ = vcat(Matrix(H_keyed), zeros(size(H_keyed)))
    Q = I(size(H̃,2)) # Q is just the Identity matrix. H depends on θ
    H̃QH̃ᵀ = Symmetric(H̃*Q*H̃')
    P̃₀vec = (I(size(G̃,1)^2) - kron(G̃,G̃)) \ reshape(H̃QH̃ᵀ, size(H̃,1)^2,1);
    P̃₀ = reshape(P̃₀vec, size(G̃))
    return P̃₀
end

"""
    get_log_likelihood_θ(m::MacroModelling.ℳ, fixed_parameters=Vector{Symbol}())

Calculate the sum of log prior probabilities for the parameters of a model.

This function iterates through the parameters of the model `m`. For each parameter
not listed in `fixed_parameters`, it retrieves its prior distribution using
`prior_θ(Val(Symbol(m.model_name)), Val(param), value)` and calculates the
log probability density function (logpdf) of the parameter's current value
under this prior. The sum of these logpdf values is returned.

# Arguments
- `m::MacroModelling.ℳ`: The model object, containing parameter definitions and current values.
- `fixed_parameters=Vector{Symbol}()`: A vector of symbols representing parameters
  whose prior log-likelihood should not be included in the sum (e.g., if they are fixed
  or their priors are handled elsewhere).

# Returns
- `total_loglik::Float64`: The sum of the log prior probabilities for the non-fixed parameters.
"""
function get_log_likelihood_θ(m::MacroModelling.ℳ, fixed_parameters=Vector{Symbol}())
    total_loglik = 0
    for (i, param, value) in zip(1:length(m.parameters), m.parameters, m.parameter_values)
        if param ∉ fixed_parameters
            # println(i, param, value)
            dist = prior_θ(Val(Symbol(m.model_name)), Val(param), value)
            loglik = logpdf(dist, value)
            total_loglik += loglik
        end
    end
    return total_loglik
end

"""
    get_log_likelihood_X(m::MacroModelling.ℳ, X, S, Λ, 𝛙::Diagonal, R::Diagonal, burnin=0, c=1.0)

Calculate the log-likelihood of the observed data `X` given the states `S` and measurement equation parameters.

The measurement equation is assumed to be:
`Xₜ = ΛSₜ + eₜ`
`eₜ = 𝛙eₜ₋₁ + ϵₜ`
where `ϵₜ ~ NID(0, R)`.

The function calculates the measurement errors `eₜ` and then the structural measurement
innovations `ϵₜ`. The log-likelihood is computed as the sum of the log probability
density function (logpdf) values of these innovations `ϵₜ` under a multivariate
normal distribution `MvNormal(zeros(size(X,2)), R ./ c)`. A burn-in period can be
specified to exclude initial observations from the likelihood calculation. The scaling
factor `c` can be used to adjust the variance of the measurement error distribution.

# Arguments
- `m::MacroModelling.ℳ`: The model object (currently unused in the function body but part of the signature).
- `X`: A matrix of observed data, where rows are time periods and columns are variables.
- `S`: A matrix or KeyedArray of states, where columns are time periods and rows (or Variables) correspond to model states.
- `Λ`: The observation matrix linking states to observables.
- `𝛙::Diagonal`: A diagonal matrix representing the autoregressive coefficients for the measurement errors `eₜ`.
- `R::Diagonal`: A diagonal matrix representing the covariance matrix of the structural measurement innovations `ϵₜ`.
- `burnin=0`: The number of initial observations to exclude from the likelihood calculation.
- `c=1.0`: A scaling factor for the covariance matrix `R` in the likelihood calculation. `R_effective = R / c`.

# Returns
- `loglik::Float64`: The calculated log-likelihood of the data.
"""
function get_log_likelihood_X(m::MacroModelling.ℳ, X, S, Λ, 𝛙::Diagonal, R::Diagonal, burnin=0, c=1.0)
    """
    Xₜ = ΛSₜ + eₜ
    eₜ = 𝛙eₜ₋₁ + ϵₜ where ϵₜ ~ NID(0,R)
    """
    _X = X[burnin+1:end,:]
    _S = Matrix(S)[:, burnin+1:end]
    e = _X' - Λ * _S
    ϵ = e[:,2:end] - 𝛙*e[:, 1:end-1]

    errors_dist = MvNormal(zeros(size(_X,2)), R ./ c)
    loglik = sum(logpdf(errors_dist, ϵ))

    return loglik
end