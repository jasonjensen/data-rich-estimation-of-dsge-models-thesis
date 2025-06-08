# please run load_results.jl first

# Table 3.1
print_cumulative_dm_table(
    bundle_07full.gdp_forecast_medians, 
    diff(bundle_07full.mvts_sw["RGDP"]) .* 100, 
    [1998Q1:2011Q4, 1998Q1:2007Q2, 2007Q3:2011Q4])

# Table 3.2
print_cumulative_dm_table(
    bundle_07full.inf_forecast_medians, 
    bundle_07full.mvts_sw["PGDP"],
    [1998Q1:2011Q4, 1998Q1:2007Q2, 2007Q3:2011Q4])

# Table 3.3
print_cumulative_dm_table(
    bundle_07part.gdp_forecast_medians, 
    diff(bundle_07part.mvts_sw["RGDP"]) .* 100, 
    [1998Q1:2011Q4, 1998Q1:2007Q2, 2007Q3:2011Q4])
print_cumulative_dm_table_full_vs_part(
    bundle_07full.gdp_forecast_medians, 
    bundle_07part.gdp_forecast_medians, 
    diff(bundle_07full.mvts_sw["RGDP"]) .* 100,
    [1998Q1:2011Q4, 1998Q1:2007Q2, 2007Q3:2011Q4])

# Table 3.4
print_cumulative_dm_table(
    bundle_08full.gdp_forecast_medians, 
    diff(bundle_08full.mvts_sw["RGDP"]) .* 100, 
    [2007Q3:2011Q4, 2012Q1:2019Q3, 2019Q4:2023Q3])  
    
# Table 3.5
print_cumulative_dm_table(
    bundle_08full.inf_forecast_medians, 
    bundle_08full.mvts_sw["PGDP"], 
    [2007Q3:2011Q4, 2012Q1:2019Q3, 2019Q4:2023Q3])  

# Table A.1
print_cumulative_dm_table(
    bundle_07part.inf_forecast_medians, 
    bundle_07part.mvts_sw["PGDP"], 
    [1998Q1:2011Q4, 1998Q1:2007Q2, 2007Q3:2011Q4])
print_cumulative_dm_table_full_vs_part(
    bundle_07full.inf_forecast_medians, 
    bundle_07part.inf_forecast_medians, 
    bundle_07full.mvts_sw["PGDP"],
    [1998Q1:2011Q4, 1998Q1:2007Q2, 2007Q3:2011Q4])

# Table A.2
print_cumulative_dm_table_multiple_targets(
    [bundle_19full.gdp_forecast_medians, bundle_19full.inf_forecast_medians], 
    [ diff(bundle_19full.mvts_sw["RGDP"]) .* 100,  bundle_19full.mvts_sw["PGDP"]], 
    2019Q4:2023Q3)


# Table A.3
print_cumulative_dm_table(
    bundle_08full.gdp_forecast_medians, 
    diff(bundle_08full.mvts_sw["RGDP"]) .* 100, 
    [2019Q4:2023Q3, 2020Q2:2023Q3, 2020Q3:2023Q3])

# Table A.4
print_cumulative_dm_table(
    bundle_08full.inf_forecast_medians, 
    bundle_08full.mvts_sw["PGDP"], 
    [2019Q4:2023Q3, 2020Q2:2023Q3, 2020Q3:2023Q3])

