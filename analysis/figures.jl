# please run load_results.jl first

# Figure 3.1
param_boxes_sw = param_box_plots(
    m_sw, 
    [bundle_08full.chain_sw_reg, gelfer_sw_reg_2008Q3, bundle_08full.chain_sw_dfm, gelfer_sw_dfm_2008Q3], 
    ["1. SW-Reg", "2. SW-Reg\n(Gelfer)", "3. SW-Rich", "4. SW-Rich\n(Gelfer)"], 
    [gelfer_means_sw_regular_estimation, gelfer_means_sw_regular_estimation, gelfer_means_sw_dfm_estimation, gelfer_means_sw_dfm_estimation]; 
    suptitle="Parameter posteriors, SW model vs Gelfer", 
    plot_size=(1600,900), param_range=1:10, burnin=75000)


# Figure 3.2
p_gdp_2007Q2 = make_comparison_chart(
    "RGDP", 2008Q2:2008Q3,
    m_sw, m_swff,
    bundle_07full, bundle_07part;
    n_parameter_draws=500, n_shock_draws = 100, hrz=16, 
    burnin=75000,
    data_X_start = 1986Q2,
    vartype=:growth,
    legend=:bottomright,
    ylim=(-3.05, 1.3))
# Plots.savefig(p_gdp_2007Q2, "graphs/recession_2007Q2.png")

# Figure 3.3
p_gdp_2008Q3 = make_comparison_chart(
    "RGDP", 2008Q2:2008Q3,
    m_sw, m_swff,
    bundle_08full, bundle_08part;
    n_parameter_draws=500, n_shock_draws = 100, hrz=16, 
    burnin=75000,
    data_X_start = 1986Q2,
    vartype=:growth,
    legend=:bottomright,
    ylim=(-3.05, 1.3))
# Plots.savefig(p_gdp_2008Q3, "graphs/recession_2008Q3.png")

# Figure 3.4
p_inf_2007Q2 = make_comparison_chart(
    "PGDP", 2008Q2:2008Q3,
    m_sw, m_swff,
    bundle_07full, bundle_07part;
    n_parameter_draws=500, n_shock_draws = 100, hrz=16, 
    burnin=75000,
    data_X_start = 1986Q2,
    vartype=:growth,
    legend=:bottomright,
    ylim=(-0.75, 0.5))
# Plots.savefig(p_inf_2007Q2, "graphs/recession_inf2007Q2.png")

# Figure 3.5
p_gdp_2019Q4 = generate_and_plot_forecasts(
    "RGDP", 
    [2020Q2, 2020Q4, 2021Q2], 
    m_sw, m_swff, 
    bundle_19full.df_sw, bundle_19full.df_swff, bundle_19full.mvts_sw, 
    bundle_19full.chain_sw_reg, bundle_19full.chain_swff_reg, bundle_19full.chain_sw_dfm, bundle_19full.chain_swff_dfm;
    n_parameter_draws=500, n_shock_draws=100, 
    zlb_level=bundle_19full.zlb_level, 
    ylim=(-9,7), size=(1000, 450), layout=(1,3))


# Figure 3.6
p_inf_2019Q4 = generate_and_plot_forecasts(
    "PGDP", 
    [2021Q3, 2022Q1, 2022Q3], 
    m_sw, m_swff, 
    bundle_19full.df_sw, bundle_19full.df_swff, bundle_19full.mvts_sw, 
    bundle_19full.chain_sw_reg, bundle_19full.chain_swff_reg, bundle_19full.chain_sw_dfm, bundle_19full.chain_swff_dfm;
    n_parameter_draws=500, n_shock_draws=100, 
    vartype=:level, zlb_level=bundle_19full.zlb_level, 
    legend=:topright,  ylim=(-2.5,2.5), size=(1000, 450), layout=(1,3))

# Figure 3.7
make_shock_decomp_charts(:RGDP, m_sw, m_swff, bundle_19full, 2007Q1, 2011Q4;
    suptitle="Shock decomposition of Output over the Great Recession", legend=:bottomleft) 

# Figure 3.8
make_shock_decomp_charts(:PGDP, m_sw, m_swff, bundle_24full, 2019Q4, 2024Q3;
    suptitle="Shock decomposition of Inflation over the Pandemic", legend=:topright)

# Figure B.1
 p_gdp_2024Q3 = generate_and_plot_forecasts(
    "RGDP", 
    2019Q4:2019Q4+5, 
    m_sw, m_swff, 
    bundle_24full.df_sw, bundle_24full.df_swff, bundle_24full.mvts_sw, 
    bundle_24full.chain_sw_reg, bundle_24full.chain_swff_reg, bundle_24full.chain_sw_dfm, bundle_24full.chain_swff_dfm;
    n_parameter_draws=500, n_shock_draws=100, 
    zlb_level=bundle_24full.zlb_level, 
    ylim=(-9,7), legend=:topright)

# Figure B.2
 p_inf_2024Q3 = generate_and_plot_forecasts(
    "PGDP", 
    2021Q2:2021Q2+5, 
    m_sw, m_swff, 
    bundle_24full.df_sw, bundle_24full.df_swff, bundle_24full.mvts_sw, 
    bundle_24full.chain_sw_reg, bundle_24full.chain_swff_reg, bundle_24full.chain_sw_dfm, bundle_24full.chain_swff_dfm;
    n_parameter_draws=500, n_shock_draws=100, 
    vartype=:level, zlb_level=bundle_24full.zlb_level, 
     ylim=(-0.7, 1.7), legend=:topright) 

# Figure B.3
make_shock_decomp_charts(:RGDP, m_sw, m_swff, bundle_07full, 2007Q1, 2011Q4;
    suptitle="Shock decomposition of Output over the Great Recession", legend=:bottomleft)

# Figure B.4
make_shock_decomp_charts(:RGDP, m_sw, m_swff, bundle_08full, 2007Q1, 2011Q4;
    suptitle="Shock decomposition of Output over the Great Recession", legend=:bottomleft)

# Figure B.5
make_shock_decomp_charts(:PGDP, m_sw, m_swff, bundle_08full, 2019Q4, 2024Q3;
    suptitle="Shock decomposition of Inflation over the Pandemic", legend=:topright)

# Figure B.6
make_shock_decomp_charts(:PGDP, m_sw, m_swff, bundle_19full, 2019Q4, 2024Q3;
    suptitle="Shock decomposition of Inflation over the Pandemic", legend=:bottomright) 