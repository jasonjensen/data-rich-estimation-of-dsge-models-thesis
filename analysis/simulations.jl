# please run load_results.jl first


gdp_forecasts_07full = generate_all_forecasts(
    "RGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_07full.chain_sw_reg, bundle_07full.chain_swff_reg, bundle_07full.chain_sw_dfm, bundle_07full.chain_swff_dfm, bundle_07full.chain_sw_mes, bundle_07full.chain_swff_mes,
    bundle_07full.df_sw, bundle_07full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_07full.zlb_level,
)
# jldsave("data/simulations/gdp_forecasts_full_2007Q2.jld2"; gdp_forecasts_07full)


gdp_forecast_medians_07full = Workspace(
    :sw_reg => get_growth_medians(gdp_forecasts_07full.forecasts_sw_reg),
    :swff_reg => get_growth_medians(gdp_forecasts_07full.forecasts_swff_reg),
    :sw_dfm => get_growth_medians(gdp_forecasts_07full.forecasts_sw_dfm),
    :swff_dfm => get_growth_medians(gdp_forecasts_07full.forecasts_swff_dfm),
    :sw_mes => get_growth_medians(gdp_forecasts_07full.forecasts_sw_mes),
    :swff_mes => get_growth_medians(gdp_forecasts_07full.forecasts_swff_mes),
)
# jldsave("data/simulations/gdp_forecasts_medians_2007Q2.jld2"; gdp_forecast_medians_07full)


inf_forecasts_07full = generate_all_forecasts(
    "PGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_07full.chain_sw_reg, bundle_07full.chain_swff_reg, bundle_07full.chain_sw_dfm, bundle_07full.chain_swff_dfm, bundle_07full.chain_sw_mes, bundle_07full.chain_swff_mes,
    bundle_07full.df_sw, bundle_07full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_07full.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/inf_forecasts_full_2007Q2.jld2"; inf_forecasts_07full)


inf_forecast_medians_07full = Workspace(
    :sw_reg => get_medians(inf_forecasts_07full.forecasts_sw_reg),
    :swff_reg => get_medians(inf_forecasts_07full.forecasts_swff_reg),
    :sw_dfm => get_medians(inf_forecasts_07full.forecasts_sw_dfm),
    :swff_dfm => get_medians(inf_forecasts_07full.forecasts_swff_dfm),
    :sw_mes => get_medians(inf_forecasts_07full.forecasts_sw_mes),
    :swff_mes => get_medians(inf_forecasts_07full.forecasts_swff_mes),
)
# jldsave("data/simulations/inf_forecasts_medians_2007Q2.jld2"; inf_forecast_medians_07full)



gdp_forecasts_07part = generate_all_forecasts(
    "RGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_07part.chain_sw_reg, bundle_07part.chain_swff_reg, bundle_07part.chain_sw_dfm, bundle_07part.chain_swff_dfm, bundle_07part.chain_sw_mes, bundle_07part.chain_swff_mes,
    bundle_07part.df_sw, bundle_07part.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_07part.zlb_level,
)
# jldsave("data/simulations/gdp_forecasts_full_2007Q2_constrained.jld2"; gdp_forecasts_07part)


gdp_forecast_medians_07part = Workspace(
    :sw_reg => get_growth_medians(gdp_forecasts_07part.forecasts_sw_reg),
    :swff_reg => get_growth_medians(gdp_forecasts_07part.forecasts_swff_reg),
    :sw_dfm => get_growth_medians(gdp_forecasts_07part.forecasts_sw_dfm),
    :swff_dfm => get_growth_medians(gdp_forecasts_07part.forecasts_swff_dfm),
    :sw_mes => get_growth_medians(gdp_forecasts_07part.forecasts_sw_mes),
    :swff_mes => get_growth_medians(gdp_forecasts_07part.forecasts_swff_mes),
)
# jldsave("data/simulations/gdp_forecasts_medians_2007Q2_constrained.jld2"; gdp_forecast_medians_07part)


inf_forecasts_07part = generate_all_forecasts(
    "PGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_07part.chain_sw_reg, bundle_07part.chain_swff_reg, bundle_07part.chain_sw_dfm, bundle_07part.chain_swff_dfm, bundle_07part.chain_sw_mes, bundle_07part.chain_swff_mes,
    bundle_07part.df_sw, bundle_07part.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_07part.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/inf_forecasts_full_2007Q2_constrained.jld2"; inf_forecasts_07part)


inf_forecast_medians_07part = Workspace(
    :sw_reg => get_medians(inf_forecasts_07part.forecasts_sw_reg),
    :swff_reg => get_medians(inf_forecasts_07part.forecasts_swff_reg),
    :sw_dfm => get_medians(inf_forecasts_07part.forecasts_sw_dfm),
    :swff_dfm => get_medians(inf_forecasts_07part.forecasts_swff_dfm),
    :sw_mes => get_medians(inf_forecasts_07part.forecasts_sw_mes),
    :swff_mes => get_medians(inf_forecasts_07part.forecasts_swff_mes),
)
# jldsave("data/simulations/inf_forecasts_medians_2007Q2_constrained.jld2"; inf_forecast_medians_07part)


gdp_forecasts_08full = generate_all_forecasts(
    "RGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_08full.chain_sw_reg, bundle_08full.chain_swff_reg, bundle_08full.chain_sw_dfm, bundle_08full.chain_swff_dfm, bundle_08full.chain_sw_mes, bundle_08full.chain_swff_mes,
    bundle_08full.df_sw, bundle_08full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_08full.zlb_level,
)
# jldsave("data/simulations/gdp_forecasts_full_2008Q3.jld2"; gdp_forecasts_08full)


gdp_forecast_medians_08full = Workspace(
    :sw_reg => get_growth_medians(gdp_forecasts_08full.forecasts_sw_reg),
    :swff_reg => get_growth_medians(gdp_forecasts_08full.forecasts_swff_reg),
    :sw_dfm => get_growth_medians(gdp_forecasts_08full.forecasts_sw_dfm),
    :swff_dfm => get_growth_medians(gdp_forecasts_08full.forecasts_swff_dfm),
    :sw_mes => get_growth_medians(gdp_forecasts_08full.forecasts_sw_mes),
    :swff_mes => get_growth_medians(gdp_forecasts_08full.forecasts_swff_mes),
)
# jldsave("data/simulations/gdp_forecasts_medians_2008Q3.jld2"; gdp_forecast_medians_08full)


inf_forecasts_08full = generate_all_forecasts(
    "PGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_08full.chain_sw_reg, bundle_08full.chain_swff_reg, bundle_08full.chain_sw_dfm, bundle_08full.chain_swff_dfm, bundle_08full.chain_sw_mes, bundle_08full.chain_swff_mes,
    bundle_08full.df_sw, bundle_08full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_08full.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/inf_forecasts_full_2008Q3.jld2"; inf_forecasts_08full)


inf_forecast_medians_08full = Workspace(
    :sw_reg => get_medians(inf_forecasts_08full.forecasts_sw_reg),
    :swff_reg => get_medians(inf_forecasts_08full.forecasts_swff_reg),
    :sw_dfm => get_medians(inf_forecasts_08full.forecasts_sw_dfm),
    :swff_dfm => get_medians(inf_forecasts_08full.forecasts_swff_dfm),
    :sw_mes => get_medians(inf_forecasts_08full.forecasts_sw_mes),
    :swff_mes => get_medians(inf_forecasts_08full.forecasts_swff_mes),
)
# jldsave("data/simulations/inf_forecasts_medians_2008Q3.jld2"; inf_forecast_medians_08full)



gdp_forecasts_08part = generate_all_forecasts(
    "RGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_08part.chain_sw_reg, bundle_08part.chain_swff_reg, bundle_08part.chain_sw_dfm, bundle_08part.chain_swff_dfm, bundle_08part.chain_sw_mes, bundle_08part.chain_swff_mes,
    bundle_08part.df_sw, bundle_08part.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_08part.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/gdp_forecasts_full_2008Q3_constrained.jld2"; gdp_forecasts_08part)


gdp_forecast_medians_08part = Workspace(
    :sw_reg => get_growth_medians(gdp_forecasts_08part.forecasts_sw_reg),
    :swff_reg => get_growth_medians(gdp_forecasts_08part.forecasts_swff_reg),
    :sw_dfm => get_growth_medians(gdp_forecasts_08part.forecasts_sw_dfm),
    :swff_dfm => get_growth_medians(gdp_forecasts_08part.forecasts_swff_dfm),
    :sw_mes => get_growth_medians(gdp_forecasts_08part.forecasts_sw_mes),
    :swff_mes => get_growth_medians(gdp_forecasts_08part.forecasts_swff_mes),
)
# jldsave("data/simulations/gdp_forecasts_medians_2008Q3_constrained.jld2"; gdp_forecast_medians_08part)


inf_forecasts_08part = generate_all_forecasts(
    "PGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_08part.chain_sw_reg, bundle_08part.chain_swff_reg, bundle_08part.chain_sw_dfm, bundle_08part.chain_swff_dfm, bundle_08part.chain_sw_mes, bundle_08part.chain_swff_mes,
    bundle_08part.df_sw, bundle_08part.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_08part.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/inf_forecasts_full_2008Q3_constrained.jld2"; inf_forecasts_08part)


inf_forecast_medians_08part = Workspace(
    :sw_reg => get_medians(inf_forecasts_08part.forecasts_sw_reg),
    :swff_reg => get_medians(inf_forecasts_08part.forecasts_swff_reg),
    :sw_dfm => get_medians(inf_forecasts_08part.forecasts_sw_dfm),
    :swff_dfm => get_medians(inf_forecasts_08part.forecasts_swff_dfm),
    :sw_mes => get_medians(inf_forecasts_08part.forecasts_sw_mes),
    :swff_mes => get_medians(inf_forecasts_08part.forecasts_swff_mes),
)
# jldsave("data/simulations/inf_forecasts_medians_2008Q3_constrained.jld2"; inf_forecast_medians_08part)



gdp_forecasts_19full = generate_all_forecasts(
    "RGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_19full.chain_sw_reg, bundle_19full.chain_swff_reg, bundle_19full.chain_sw_dfm, bundle_19full.chain_swff_dfm, bundle_19full.chain_sw_mes, bundle_19full.chain_swff_mes,
    bundle_19full.df_sw, bundle_19full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_19full.zlb_level,
)
# jldsave("data/simulations/gdp_forecasts_full_2019Q4.jld2"; gdp_forecasts_19full)


gdp_forecast_medians_19full = Workspace(
    :sw_reg => get_growth_medians(gdp_forecasts_19full.forecasts_sw_reg),
    :swff_reg => get_growth_medians(gdp_forecasts_19full.forecasts_swff_reg),
    :sw_dfm => get_growth_medians(gdp_forecasts_19full.forecasts_sw_dfm),
    :swff_dfm => get_growth_medians(gdp_forecasts_19full.forecasts_swff_dfm),
    :sw_mes => get_growth_medians(gdp_forecasts_19full.forecasts_sw_mes),
    :swff_mes => get_growth_medians(gdp_forecasts_19full.forecasts_swff_mes),
)
# jldsave("data/simulations/gdp_forecasts_medians_2019Q4.jld2"; gdp_forecast_medians_19full)


inf_forecasts_19full = generate_all_forecasts(
    "PGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_19full.chain_sw_reg, bundle_19full.chain_swff_reg, bundle_19full.chain_sw_dfm, bundle_19full.chain_swff_dfm, bundle_19full.chain_sw_mes, bundle_19full.chain_swff_mes,
    bundle_19full.df_sw, bundle_19full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_19full.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/inf_forecasts_full_2019Q4.jld2"; inf_forecasts_19full)


inf_forecast_medians_19full = Workspace(
    :sw_reg => get_medians(inf_forecasts_19full.forecasts_sw_reg),
    :swff_reg => get_medians(inf_forecasts_19full.forecasts_swff_reg),
    :sw_dfm => get_medians(inf_forecasts_19full.forecasts_sw_dfm),
    :swff_dfm => get_medians(inf_forecasts_19full.forecasts_swff_dfm),
    :sw_mes => get_medians(inf_forecasts_19full.forecasts_sw_mes),
    :swff_mes => get_medians(inf_forecasts_19full.forecasts_swff_mes),
)
# jldsave("data/simulations/inf_forecasts_medians_2019Q4.jld2"; inf_forecast_medians_19full)


gdp_forecasts_24full = generate_all_forecasts(
    "RGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_24full.chain_sw_reg, bundle_24full.chain_swff_reg, bundle_24full.chain_sw_dfm, bundle_24full.chain_swff_dfm, bundle_24full.chain_sw_mes, bundle_24full.chain_swff_mes,
    bundle_24full.df_sw, bundle_24full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_24full.zlb_level,
)
# jldsave("data/simulations/gdp_forecasts_full_2024Q3.jld2"; gdp_forecasts_24full)


gdp_forecast_medians_24full = Workspace(
    :sw_reg => get_growth_medians(gdp_forecasts_24full.forecasts_sw_reg),
    :swff_reg => get_growth_medians(gdp_forecasts_24full.forecasts_swff_reg),
    :sw_dfm => get_growth_medians(gdp_forecasts_24full.forecasts_sw_dfm),
    :swff_dfm => get_growth_medians(gdp_forecasts_24full.forecasts_swff_dfm),
    :sw_mes => get_growth_medians(gdp_forecasts_24full.forecasts_sw_mes),
    :swff_mes => get_growth_medians(gdp_forecasts_24full.forecasts_swff_mes),
)
# jldsave("data/simulations/gdp_forecasts_medians_2024Q3.jld2"; gdp_forecast_medians_24full)


inf_forecasts_24full = generate_all_forecasts(
    "PGDP", 1998Q1:2023Q3,
    m_sw, m_swff,
    bundle_24full.chain_sw_reg, bundle_24full.chain_swff_reg, bundle_24full.chain_sw_dfm, bundle_24full.chain_swff_dfm, bundle_24full.chain_sw_mes, bundle_24full.chain_swff_mes,
    bundle_24full.df_sw, bundle_24full.df_swff;
    n_parameter_draws=500, n_shock_draws= 100, horizon = 4,
    data_X_start = 1986Q2, zlb_level=bundle_24full.zlb_level,
    vartype=:level, 
)
# jldsave("data/simulations/inf_forecasts_full_2024Q3.jld2"; inf_forecasts_24full)


inf_forecast_medians_24full = Workspace(
    :sw_reg => get_medians(inf_forecasts_24full.forecasts_sw_reg),
    :swff_reg => get_medians(inf_forecasts_24full.forecasts_swff_reg),
    :sw_dfm => get_medians(inf_forecasts_24full.forecasts_sw_dfm),
    :swff_dfm => get_medians(inf_forecasts_24full.forecasts_swff_dfm),
    :sw_mes => get_medians(inf_forecasts_24full.forecasts_sw_mes),
    :swff_mes => get_medians(inf_forecasts_24full.forecasts_swff_mes),
)
# jldsave("data/simulations/inf_forecasts_medians_2024Q3.jld2"; inf_forecast_medians_24full)