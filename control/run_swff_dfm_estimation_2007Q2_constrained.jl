include("$(pwd())/models/Gelfer_SWFF_COEF.jl")
m = deepcopy(SWFF);

fred_df_est, fred_mvts_est = load_and_treat_data(1986Q2:2007Q2, 1986Q2:2007Q2, fred_groups_swff);
data = get_data_obj(m, fred_df_est, false, ["NASDAQ", "WTI"])

Σ⁻¹ = get_hessian_dfm(m)
Σ⁻¹ = make_psd(Σ⁻¹)
isposdef(Σ⁻¹)

disable_timer!(to)
@timeit to "MWG" m_final, res = metropolis_within_gibbs(m, data, Σ⁻¹; 
    niter=150000, 
    burnin=0,
    c=1e-1,
    w=1e-3,
    c_min=1e-7,
    target_acceptance_rate = 0.25,
    Λ_constraints = Λ_constraints_SWFF,
    verbose = false,
    fixed_parameters=[:β, :α, :τ, :Iy, :Gy, :λw, :γ, :F, :S_ss],
    progress_plots = true,
    speedup=true,
    # restore_from="data/results/partial_gelfer_2025-05-10T2054_RESTART.jld2",
    # restore_from_current_iter=5000
);

println("Done.")

store_results(res,m, "swff_dfm_2007Q2_constrained")

println("Saved.")