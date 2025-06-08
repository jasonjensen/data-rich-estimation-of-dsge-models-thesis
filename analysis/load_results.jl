include("$(pwd())/models/Gelfer_SW.jl")
include("$(pwd())/models/Gelfer_SWFF_COEF.jl")

m_sw = deepcopy(SW)
m_swff = deepcopy(SWFF)


# 2007Q2
fred_df_sw_07full, fred_mvts_swff_07full = load_and_treat_data(1986Q2:2007Q2, 1986Q2:2007Q2, fred_groups_sw);
fred_df_swff_07full, fred_mvts_swff_07full= load_and_treat_data(1986Q2:2007Q2, 1986Q2:2007Q2, fred_groups_swff);
data_sw_07full = get_data_obj(m_sw, fred_df_sw_07full)
data_sw_core_07full = get_data_obj(m_sw, fred_df_sw_07full, true)
data_swff_07full = get_data_obj(m_swff, fred_df_swff_07full)
data_swff_core_07full = get_data_obj(m_swff, fred_df_swff_07full, true)
fred_df_sw_extended_07full, fred_mvts_sw_extended_07full = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2007Q2, fred_groups_sw);
fred_df_swff_extended_07full, fred_mvts_swff_extended_07full= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2007Q2, fred_groups_swff);
ZLB_LEVEL_07full = fred_mvts_sw_extended_07full.FedFunds[2008Q4]
chain_swff_dfm_07full = load_results("data/chains/2007Q2/swff_dfm_2007Q2_2025-05-11T0300.daec")
chain_swff_reg_07full = load_results("data/chains/2007Q2/swff_reg_2007Q2_2025-05-05T1530.daec")
chain_swff_mes_07full = load_results("data/chains/2007Q2/swff_mes_2007Q2_2025-05-05T1548.daec")
chain_sw_dfm_07full = load_results("data/chains/2007Q2/sw_dfm_2007Q2_2025-05-11T1014.daec")
chain_sw_reg_07full = load_results("data/chains/2007Q2/sw_reg_2007Q2_2025-05-05T0955.daec")
chain_sw_mes_07full = load_results("data/chains/2007Q2/sw_mes_2007Q2_2025-05-05T1240.daec")
gdp_forecast_medians_07full = load("data/simulations/gdp_forecasts_medians_2007Q2.jld2")["gdp_forecast_medians"]
inf_forecast_medians_07full = load("data/simulations/inf_forecasts_medians_2007Q2.jld2")["inf_forecast_medians"]


# 07 constrained
fred_df_sw_07part, fred_mvts_sw_07part = load_and_treat_data(1986Q2:2007Q2, 1986Q2:2007Q2, fred_groups_sw);
fred_df_swff_07part, fred_mvts_swff_07part= load_and_treat_data(1986Q2:2007Q2, 1986Q2:2007Q2, fred_groups_swff);
data_sw_07part = get_data_obj(m_sw, fred_df_sw_07part, false, ["NASDAQ", "WTI"])
data_sw_core_07part = get_data_obj(m_sw, fred_df_sw_07part, true, ["NASDAQ", "WTI"])
data_swff_07part = get_data_obj(m_swff, fred_df_swff_07part, false, ["NASDAQ", "WTI"])
data_swff_core_07part = get_data_obj(m_swff, fred_df_swff_07part, true, ["NASDAQ", "WTI"])
fred_df_sw_extended_07part, fred_mvts_sw_extended_07part = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2007Q2, fred_groups_sw);
fred_df_swff_extended_07part, fred_mvts_swff_extended_07part= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2007Q2, fred_groups_swff);
ZLB_LEVEL_07part = fred_mvts_sw_extended_07part.FedFunds[2008Q4]
chain_swff_dfm_07part = load_results("data/chains/2007Q2/swff_dfm_2007Q2_constrained_2025-05-12T0459.daec")
chain_swff_reg_07part = load_results("data/chains/2007Q2/swff_reg_2007Q2_2025-05-05T1530.daec")
chain_swff_mes_07part = load_results("data/chains/2007Q2/swff_mes_2007Q2_2025-05-05T1548.daec")
chain_sw_dfm_07part = load_results("data/chains/2007Q2/sw_dfm_2007Q2_constrained_2025-05-12T0359.daec")
chain_sw_reg_07part = load_results("data/chains/2007Q2/sw_reg_2007Q2_2025-05-05T0955.daec")
chain_sw_mes_07part = load_results("data/chains/2007Q2/sw_mes_2007Q2_2025-05-05T1240.daec")
gdp_forecast_medians_07part = load("data/simulations/gdp_forecasts_medians_2007Q2_constrained.jld2")["gdp_forecast_medians"]
inf_forecast_medians_07part = load("data/simulations/inf_forecasts_medians_2007Q2_constrained.jld2")["inf_forecast_medians"]

bundle_07full = Workspace(
    :chain_swff_dfm    => chain_swff_dfm_07full,  
    :chain_swff_reg    => chain_swff_reg_07full,  
    :chain_swff_mes    => chain_swff_mes_07full,  
    :chain_sw_dfm      => chain_sw_dfm_07full,  
    :chain_sw_reg      => chain_sw_reg_07full,  
    :chain_sw_mes      => chain_sw_mes_07full,
    :zlb_level       => ZLB_LEVEL_07full,
    :df_sw           => fred_df_sw_extended_07full,
    :df_swff         => fred_df_swff_extended_07full,
    :mvts_sw         => fred_mvts_sw_extended_07full,
    :gdp_forecast_medians => gdp_forecast_medians_07full,
    :inf_forecast_medians => inf_forecast_medians_07full,
)

bundle_07part = Workspace(
    :chain_swff_dfm    => chain_swff_dfm_07part,  
    :chain_swff_reg    => chain_swff_reg_07part,  
    :chain_swff_mes    => chain_swff_mes_07part,  
    :chain_sw_dfm      => chain_sw_dfm_07part,  
    :chain_sw_reg      => chain_sw_reg_07part,  
    :chain_sw_mes      => chain_sw_mes_07part,
    :zlb_level       => ZLB_LEVEL_07part,
    :df_sw           => fred_df_sw_extended_07part,
    :df_swff         => fred_df_swff_extended_07part,
    :mvts_sw         => fred_mvts_sw_extended_07part,
    :excluded_series  => ["NASDAQ", "WTI"],
    :gdp_forecast_medians => gdp_forecast_medians_07part,
    :inf_forecast_medians => inf_forecast_medians_07part,
)

# 08 full
fred_df_sw_08full, fred_mvts_swff_08full = load_and_treat_data(1986Q2:2008Q3, 1986Q2:2008Q3, fred_groups_sw);
fred_df_swff_08full, fred_mvts_swff_08full= load_and_treat_data(1986Q2:2008Q3, 1986Q2:2008Q3, fred_groups_swff);
data_sw_08full = get_data_obj(m_sw, fred_df_sw_08full)
data_sw_core_08full = get_data_obj(m_sw, fred_df_sw_08full, true)
data_swff_08full = get_data_obj(m_swff, fred_df_swff_08full)
data_swff_core_08full = get_data_obj(m_swff, fred_df_swff_08full, true)
fred_df_sw_extended_08full, fred_mvts_sw_extended_08full = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2008Q3, fred_groups_sw);
fred_df_swff_extended_08full, fred_mvts_swff_extended_08full= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2008Q3, fred_groups_swff);
ZLB_LEVEL_08full = fred_mvts_sw_extended_08full.FedFunds[2008Q4]
chain_swff_dfm_08full = load_results("data/chains/2008Q3/swff_dfm_2008Q3_2025-05-11T1337.daec")
chain_swff_reg_08full = load_results("data/chains/2008Q3/swff_reg_2008Q3_2025-05-03T1434.daec")
chain_swff_mes_08full = load_results("data/chains/2008Q3/swff_mes_2008Q3_2025-05-11T2022.daec")
chain_sw_dfm_08full = load_results("data/chains/2008Q3/sw_dfm_2008Q3_2025-05-10T1956.daec")
chain_sw_reg_08full = load_results("data/chains/2008Q3/sw_reg_2008Q3_2025-05-03T1122.daec")
chain_sw_mes_08full = load_results("data/chains/2008Q3/sw_mes_2008Q3_2025-04-10T0441.daec")
gdp_forecast_medians_08full = load("data/simulations/gdp_forecasts_medians_2008Q3.jld2")["gdp_forecast_medians"]
inf_forecast_medians_08full = load("data/simulations/inf_forecasts_medians_2008Q3.jld2")["inf_forecast_medians"]

# 08 constrained
fred_df_sw_08part, fred_mvts_sw_08part = load_and_treat_data(1986Q2:2008Q3, 1986Q2:2008Q3, fred_groups_sw);
fred_df_swff_08part, fred_mvts_swff_08part= load_and_treat_data(1986Q2:2008Q3, 1986Q2:2008Q3, fred_groups_swff);
data_sw_08part = get_data_obj(m_sw, fred_df_sw_08part, false, ["NASDAQ", "WTI"])
data_sw_core_08part = get_data_obj(m_sw, fred_df_sw_08part, true, ["NASDAQ", "WTI"])
data_swff_08part = get_data_obj(m_swff, fred_df_swff_08part, false, ["NASDAQ", "WTI"])
data_swff_core_08part = get_data_obj(m_swff, fred_df_swff_08part, true, ["NASDAQ", "WTI"])
fred_df_sw_extended_08part, fred_mvts_sw_extended_08part = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2008Q3, fred_groups_sw);
fred_df_swff_extended_08part, fred_mvts_swff_extended_08part= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2008Q3, fred_groups_swff);
ZLB_LEVEL_08part = fred_mvts_sw_extended_08part.FedFunds[2008Q4]
chain_swff_dfm_08part = load_results("data/chains/2008Q3/swff_dfm_2008Q3_constrained_2025-05-12T1423.daec")
chain_swff_reg_08part = load_results("data/chains/2008Q3/swff_reg_2008Q3_2025-05-03T1434.daec")
chain_swff_mes_08part = load_results("data/chains/2008Q3/swff_mes_2008Q3_2025-05-11T2022.daec")
chain_sw_dfm_08part = load_results("data/chains/2008Q3/sw_dfm_2008Q3_constrained_2025-05-12T1313.daec")
chain_sw_reg_08part = load_results("data/chains/2008Q3/sw_reg_2008Q3_2025-05-03T1122.daec")
chain_sw_mes_08part = load_results("data/chains/2008Q3/sw_mes_2008Q3_2025-04-10T0441.daec")
gdp_forecast_medians_08part = load("data/simulations/gdp_forecasts_medians_2008Q3_constrained.jld2")["gdp_forecast_medians"]
inf_forecast_medians_08part = load("data/simulations/inf_forecasts_medians_2008Q3_constrained.jld2")["inf_forecast_medians"]

bundle_08full = Workspace(
    :chain_swff_dfm    => chain_swff_dfm_08full,  
    :chain_swff_reg    => chain_swff_reg_08full,  
    :chain_swff_mes    => chain_swff_mes_08full,  
    :chain_sw_dfm      => chain_sw_dfm_08full,  
    :chain_sw_reg      => chain_sw_reg_08full,  
    :chain_sw_mes      => chain_sw_mes_08full,
    :zlb_level       => ZLB_LEVEL_08full,
    :df_sw           => fred_df_sw_extended_08full,
    :df_swff         => fred_df_swff_extended_08full,
    :mvts_sw         => fred_mvts_sw_extended_08full,
    :gdp_forecast_medians => gdp_forecast_medians_08full,
    :inf_forecast_medians => inf_forecast_medians_08full,
)

bundle_08part = Workspace(
    :chain_swff_dfm    => chain_swff_dfm_08part,  
    :chain_swff_reg    => chain_swff_reg_08part,  
    :chain_swff_mes    => chain_swff_mes_08part,  
    :chain_sw_dfm      => chain_sw_dfm_08part,  
    :chain_sw_reg      => chain_sw_reg_08part,  
    :chain_sw_mes      => chain_sw_mes_08part,
    :zlb_level       => ZLB_LEVEL_08part,
    :df_sw           => fred_df_sw_extended_08part,
    :df_swff         => fred_df_swff_extended_08part,
    :mvts_sw         => fred_mvts_sw_extended_08part,
    :excluded_series  => ["NASDAQ", "WTI"],
    :gdp_forecast_medians => gdp_forecast_medians_08part,
    :inf_forecast_medians => inf_forecast_medians_08part,
)


# 19 full
fred_df_sw_19full, fred_mvts_swff_19full = load_and_treat_data(1986Q2:2019Q4, 1986Q2:2019Q4, fred_groups_sw);
fred_df_swff_19full, fred_mvts_swff_19full= load_and_treat_data(1986Q2:2019Q4, 1986Q2:2019Q4, fred_groups_swff);
data_sw_19full = get_data_obj(m_sw, fred_df_sw_19full)
data_sw_core_19full = get_data_obj(m_sw, fred_df_sw_19full, true)
data_swff_19full = get_data_obj(m_swff, fred_df_swff_19full)
data_swff_core_19full = get_data_obj(m_swff, fred_df_swff_19full, true)
fred_df_sw_extended_19full, fred_mvts_sw_extended_19full = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2019Q4, fred_groups_sw);
fred_df_swff_extended_19full, fred_mvts_swff_extended_19full= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2019Q4, fred_groups_swff);
ZLB_LEVEL_19full = fred_mvts_sw_extended_19full.FedFunds[2008Q4]
chain_swff_dfm_19full = load_results("data/chains/2019Q4/swff_dfm_2019Q4_2025-05-30T1629.daec")
chain_swff_reg_19full = load_results("data/chains/2019Q4/swff_reg_2019Q4_2025-05-29T2225.daec")
chain_swff_mes_19full = load_results("data/chains/2019Q4/swff_mes_2019Q4_2025-05-30T2107.daec")
chain_sw_dfm_19full =   load_results("data/chains/2019Q4/sw_dfm_2019Q4_2025-05-30T1541.daec")
chain_sw_reg_19full =   load_results("data/chains/2019Q4/sw_reg_2019Q4Q3_2025-05-29T2209.daec")
chain_sw_mes_19full =   load_results("data/chains/2019Q4/sw_mes_2019Q4_2025-05-30T2055.daec")
gdp_forecast_medians_19full = load("data/simulations/gdp_forecasts_medians_2019Q4.jld2")["gdp_forecast_medians"]
inf_forecast_medians_19full = load("data/simulations/inf_forecasts_medians_2019Q4.jld2")["inf_forecast_medians"]

bundle_19full = Workspace(
    :chain_swff_dfm    => chain_swff_dfm_19full,  
    :chain_swff_reg    => chain_swff_reg_19full,  
    :chain_swff_mes    => chain_swff_mes_19full,  
    :chain_sw_dfm      => chain_sw_dfm_19full,  
    :chain_sw_reg      => chain_sw_reg_19full,  
    :chain_sw_mes      => chain_sw_mes_19full,
    :zlb_level       => ZLB_LEVEL_19full,
    :df_sw           => fred_df_sw_extended_19full,
    :df_swff         => fred_df_swff_extended_19full,
    :mvts_sw         => fred_mvts_sw_extended_19full,
    :gdp_forecast_medians => gdp_forecast_medians_19full,
    :inf_forecast_medians => inf_forecast_medians_19full,
)

# 24 full
fred_df_sw_24full, fred_mvts_swff_24full = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2024Q3, fred_groups_sw);
fred_df_swff_24full, fred_mvts_swff_24full= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2024Q3, fred_groups_swff);
data_sw_24full = get_data_obj(m_sw, fred_df_sw_24full)
data_sw_core_24full = get_data_obj(m_sw, fred_df_sw_24full, true)
data_swff_24full = get_data_obj(m_swff, fred_df_swff_24full)
data_swff_core_24full = get_data_obj(m_swff, fred_df_swff_24full, true)
fred_df_sw_extended_24full, fred_mvts_sw_extended_24full = load_and_treat_data(1986Q2:2024Q3, 1986Q2:2024Q3, fred_groups_sw);
fred_df_swff_extended_24full, fred_mvts_swff_extended_24full= load_and_treat_data(1986Q2:2024Q3, 1986Q2:2024Q3, fred_groups_swff);
ZLB_LEVEL_24full = fred_mvts_sw_extended_24full.FedFunds[2008Q4]
chain_swff_dfm_24full = load_results("data/chains/2024Q3/swff_dfm_2024_2025-05-11T1659.daec")
chain_swff_reg_24full = load_results("data/chains/2024Q3/swff_reg_2024Q3_2025-05-03T2111.daec")
chain_swff_mes_24full = load_results("data/chains/2024Q3/swff_mes_2024Q4_2025-05-01T2339.daec")
chain_sw_dfm_24full = load_results("data/chains/2024Q3/sw_dfm_2024Q3_2025-05-13T0912.daec")
chain_sw_reg_24full = load_results("data/chains/2024Q3/sw_reg_2024Q3_2025-05-03T1649.daec")
chain_sw_mes_24full = load_results("data/chains/2024Q3/sw_mes_2024Q4_2025-04-29T2341.daec")
gdp_forecast_medians_24full = load("data/simulations/gdp_forecasts_medians_2024Q3.jld2")["gdp_forecast_medians"]
inf_forecast_medians_24full = load("data/simulations/inf_forecasts_medians_2024Q3.jld2")["inf_forecast_medians"]


bundle_24full = Workspace(
    :chain_swff_dfm    => chain_swff_dfm_24full,  
    :chain_swff_reg    => chain_swff_reg_24full,  
    :chain_swff_mes    => chain_swff_mes_24full,  
    :chain_sw_dfm      => chain_sw_dfm_24full,  
    :chain_sw_reg      => chain_sw_reg_24full,  
    :chain_sw_mes      => chain_sw_mes_24full,
    :zlb_level       => ZLB_LEVEL_24full,
    :df_sw           => fred_df_sw_extended_24full,
    :df_swff         => fred_df_swff_extended_24full,
    :mvts_sw         => fred_mvts_sw_extended_24full,
    :gdp_forecast_medians => gdp_forecast_medians_24full,
    :inf_forecast_medians => inf_forecast_medians_24full,
)



# Gelfer rerun
gelfer_sw_reg_2008Q3_mat = MAT.matread("data/chains/gelfer_rerun/sw_reg.mat")["posterior"]
gelfer_swff_reg_2008Q3_mat = MAT.matread("data/chains/gelfer_rerun/swff_reg.mat")["posterior"]
gelfer_sw_dfm_2008Q3_mat = MAT.matread("data/chains/gelfer_rerun/DSGE_DFM_SW_estimates.mat")["posterior"]
gelfer_swff_dfm_2008Q3_mat = MAT.matread("data/chains/gelfer_rerun/DSGE_DFM_SWFF_estimates.mat")["posterior"]

gelfer_sw_reg_2008Q3 = process_gelfer_mat(m_sw, gelfer_sw_reg_2008Q3_mat, gelfer_map_sw, 150000)
gelfer_swff_reg_2008Q3 = process_gelfer_mat(m_swff, gelfer_swff_reg_2008Q3_mat, gelfer_map_swff, 150000)
gelfer_sw_dfm_2008Q3 = process_gelfer_mat(m_sw, gelfer_sw_dfm_2008Q3_mat, gelfer_map_sw, 150000)
gelfer_swff_dfm_2008Q3 = process_gelfer_mat(m_swff, gelfer_swff_dfm_2008Q3_mat, gelfer_map_swff, 150000)

