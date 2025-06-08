include("$(pwd())/models/Gelfer_SW.jl")
m = deepcopy(SW);

fred_df_est, fred_mvts_est = load_and_treat_data(1986Q2:2008Q3, 1986Q2:2008Q3, fred_groups_sw);
data = get_data_obj(m, fred_df_est)

Σ⁻¹ = get_hessian_dfm(m)
Σ⁻¹ = make_psd(Σ⁻¹)
isposdef(Σ⁻¹)

disable_timer!(to)
@timeit to "MWG" m_final, res = metropolis_within_gibbs(m, data, Σ⁻¹; 
    niter=150000, 
    burnin=0,
    c=1e-2,
    w=1e-3,
    c_min=1e-7,
    target_acceptance_rate = 0.25,
    Λ_constraints = Λ_constraints_SW,
    verbose = false,
    fixed_parameters=[:β, :α, :τ, :Iy, :Gy, :λw],
    progress_plots = true,
    speedup=true
);

println("Done.")

store_results(res,m, "sw_dfm_2008Q3")

println("Saved.")


# θ_for_hessian = copy(SW.parameter_values)
# _, _, 𝐒_for_hessian, _, _ = MacroModelling.get_relevant_steady_state_and_state_update(Val(:first_order), θ_for_hessian, m)
# M = get_measurement_matrix(m, 𝐒_for_hessian, data, core=true)
# A = get_measurement_const(m, θ_for_hessian, data)
# SR = I(length(data.core_series));
# Q = I(length(m.exo))
# S, _ = kalman(m, 𝐒_for_hessian, data, M, SR, A, Q; core=true, backward=false)
# R_Θ, Λ_μ, 𝛙_μ, 𝛙_σ = get_initial_Γ(S, data, Λ_constraints_SW)
# for i = 1:10
#     S, loglikX = generate_states(m, θ_for_hessian, 𝐒_for_hessian, data, Λ_μ, 𝛙_μ, R_Θ) #should I remove the noise added here?
#     R_Θ, Λ_μ, 𝛙_μ, 𝛙_σ = get_initial_Γ(S, data, Λ_constraints_SW)
# end
# H = get_hessian(m, θ_for_hessian, Λ_μ, 𝛙_μ, R_Θ, data, [:β, :α, :τ, :Iy, :Gy, :λw,], 0.005, true)

