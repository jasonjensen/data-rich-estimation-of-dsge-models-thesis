fred_codes = Dict{String,String}(
    "RGDP"=>"GDPC1",
    "IP_TOTAL"=>"INDPRO",
    "RGDI"=>"A261RX1Q020SBEA",
    "PGDP"=>"GDPDEF",
    "PCED"=>"PCECTPI",
    "CPI_ALL"=>"CPIAUCSL",
    "RCONS"=>"PCECC96",
    "RINV"=>"GPDIC1", # WAS GDPI
    "RWAGE"=>"AHETPI",
    "HOURS"=>"HOANBS",
    "EMP_CES"=>"PAYEMS+USGOVT",
    "EMP_CPS"=>"USPRIV", # changed from CE160V, is now "All Employeed, Total Private"
    "FedFunds"=>"FEDFUNDS",
    "Tbill_3m"=>"TB3MS",
    "AAABond"=>"AAA",
    "SFYBAAC"=>"BAA10Y",
    "SFYAAAC"=>"AAA10Y",
    "IP_FINAL"=>"IPFINAL", # changed from "IPS299",
    "IP_CONS_DBLE"=>"IPDCONGD",
    "IP_CONS_NONDBLE"=>"IPNCONGD",
    "IP_BUS_EQPT"=>"IPBUSEQ",
    "IP_DRBLE_MATS"=>"IPDMAT",
    "IP_NONDRBLE_MATS"=>"IPNMAT",
    "IP_MFG"=>"IPMAN",
    "IP_FUELS"=>"IPUTIL",
    "PMP"=>"NAPMPI",
    "RCONS_DRBLE"=>"DDURRA3Q086SBEA", # This series has the wrong code, it us the price series. Correct series is PCEDGC96, but that series is too short.
    "RCONS_NONDRBLE"=>"DNDGRA3Q086SBEA",  # This series has the wrong code, it us the price series. Correct series is PCENDC96, but that series is too short.
    "RCONS_SERV"=>"DSERRA3Q086SBEA", # This series has the wrong code, it us the price series. Correct series is PCESC96, but that series is too short.
    "REXPORTS"=>"B020RA3Q086SBEA",
    "RIMPORTS"=>"B255RA3Q086SBEA",
    "RGOV"=>"B823RA3Q086SBEA",
    "EMP_Mining"=>"USMINE",
    "EMP_CONST"=>"USCONS",
    "EMP_MFG"=>"MANEMP",
    "EMP_SERVICES"=>"SRVPRD",
    "EMP_TTU"=>"USTPU",
    "EMP_WHOLESALE"=>"USWTRADE",
    "EMP_RETAIL"=>"USTRADE",
    "EMP_FIN"=>"USFIRE",
    "EMP_GOVT"=>"USGOVT",
    "EMP_PROSERV"=>"USPBS",
    "EMP_LEISURE"=>"USLAH",
    "URATE"=>"UNRATE",
    "U_DURATION"=>"UEMPMEAN",
    "U_L5WKS"=>"UEMPLT5",
    "U_5_14WKS"=>"UEMP5TO14",
    "U_15_26WKS"=>"UEMP15T26",
    "U_M27WKS"=>"UEMP27OV",
    "HOURS_AVG"=>"CES0600000007",
    "HOURS_AVG_OT"=>"AWOTMAN",
    "HSTARTS_NE"=>"HOUSTNE",
    "HSTARTS_MW"=>"HOUSTMW",
    "HSTARTS_SOU"=>"HOUSTS",
    "HSTARTS_WST"=>"HOUSTW",
    # "RRRESINV"=>"B012RA3A086NBEA", # orig: B011RA3Q086SBEA (was quarterly)
    "SFYGM6"=>"TB6MS-TB3MS",
    "SFYGT1"=>"GS1-TB3MS",
    "SFYGT10"=>"GS10-TB3MS",
    "TOT_RES"=>"TOTRESNS",
    "TOT_RES_NB"=>"NONBORRES",
    "BUS_LOANS"=>"BUSLOANS",
    "CONS_CREDIT"=>"NONREVSL",
    # "SP500"=>"SP500",
    "NASDAQ"=>"NASDAQCOM",
    # "DJIA"=>"DJIA",
    "DJIA"=>"M1109BUSM293NNBR",
    # "EXR_US"=>"TWEXMMTH", # alt:  (needs splicing!)
    "EXR_SW"=>"EXSZUS",
    "EXR_JAN"=>"EXJPUS",
    "EXR_UK"=>"EXUSUK",
    "EXR_CAN"=>"EXCAUS",
    "NAPMI"=>"NAPM",
    "NAPM_NEW_ORDERS"=>"NAPMNOI",
    "NAPM_SUP_DEL"=>"NAPMSDI",
    "NAPM_INVENTORIES"=>"NAPMII",
    "RNONRESINV"=>"B009RA3Q086SBEA",
    "RAHE_CONST"=>"CES3000000008", # This series has the wrong code, it is a duplicate of RAHE_CONST. Correct code is CES2000000008
    "RAHE_MFG"=>"CES3000000008",
    "RCOMP_HR"=>"COMPRNFB",
    "ULC"=>"ULCNFB",
    "CPI_CORE"=>"CPILFESL",
    "PCED_DUR"=>"DDURRA3Q086SBEA",
    "PCED_NDUR"=>"DNDGRA3Q086SBEA",
    "PCED_SERV"=>"DSERRG3Q086SBEA",
    "PINV_GDP"=>"GPDICTPI",
    "PINV_NRES_STRUCT"=>"B009RG3Q086SBEA",
    # "PINV_NRES_EQP"=>"B010RG3Q086SBEA",
    "PINV_RES"=>"B011RG3Q086SBEA",
    "PEXPORTS"=>"B020RG3Q086SBEA",
    "PIMPORTS"=>"B021RG3Q086SBEA",
    "PGOV"=>"B822RG3Q086SBEA",
    "P_COM"=>"PPIACO",
    # "P_OIL"=>"PPICEM/PCEPILFE",
    "UTL11"=>"MCUMFN",
    "LABOR_PROD"=>"OPHNFB",
    "UMICH_CONS"=>"UMCSENT",
    "M_1"=>"M1SL",
    "M_2"=>"M2SL",
    # new series to make up for missing ones
    "WTI"=>"DCOILWTICO/PCEPILFE",
    "IP"=>"INDPRO",
    "MAN_NEW_ORDERS" => "AMTMNO",
    "INVENT" => "BUSINV",
)

fred_descriptions = Dict{String,String}(
    "RGDP"=>"Real GDP",
    "IP_TOTAL"=>"Industrial Production Index:total",
    "RGDI"=>"Real Domestic Income",
    "PGDP"=>"GDP Price deflator",
    "PCED"=>"PCE_ALL Price deflator",
    "CPI_ALL"=>"CPI_ALL Price index",
    "RCONS"=>"Real Personal Consumption Expenditures",
    "RINV"=>"Real Private Domestic Investment",
    "RWAGE"=>"Real Average Hourly wages:production:total private",
    "HOURS"=>"Hours Worked",
    "EMP_CES"=>"Employees:Total Nonfarm",
    "EMP_CPS"=>"Civilian Labor Force:Employed, Total",
    "FedFunds"=>"Federal Funds Rate (effective)",
    "Tbill_3m"=>"Interest Rate U.S. Treasury bill 3 month",
    "AAABond"=>"Bond Yield: Moody’s AAA corporate",
    "SFYBAAC"=>"Spread of BAA corporate yield to 10 year Treasury",
    "SFYAAAC"=>"Spread of AAA corporate yield to 10 year Treasury",
    "IP_FINAL"=>"Industrial Production Index:final products",
    "IP_CONS_DBLE"=>"Industrial Production Index:Durable Consumer Goods",
    "IP_CONS_NONDBLE"=>"Industrial Production Index:NonDurable Consumer Goods",
    "IP_BUS_EQPT"=>"Industrial Production Index:Business Equipment",
    "IP_DRBLE_MATS"=>"Industrial Production Index:Durable Goods Materials",
    "IP_NONDRBLE_MATS"=>"Industrial Production Index:NonDurable Goods Materials",
    "IP_MFG"=>"Industrial Production Index:Manufacturing",
    "IP_FUELS"=>"Industrial Production Index:Fuels",
    "PMP"=>"NAPM Production index",
    "RCONS_DRBLE"=>"Real Personal Consumption Expenditures index:Durables",
    "RCONS_NONDRBLE"=>"Real Personal Consumption Expenditures index:NonDurables",
    "RCONS_SERV"=>"Real Personal Consumption Expenditures index:Sevices",
    "REXPORTS"=>"Real Exports Quantity Index",
    "RIMPORTS"=>"Real Imports Quantity Index",
    "RGOV"=>"Real Government Consumption & Investment Quantity ",
    "EMP_Mining"=>"Employees:Mining & Logging",
    "EMP_CONST"=>"Employees:Construction",
    "EMP_MFG"=>"Employees:Manufacturing",
    "EMP_SERVICES"=>"Employees:Service Providing",
    "EMP_TTU"=>"Employees:Trade, Transportation, Utilities",
    "EMP_WHOLESALE"=>"Employees:Wholesale Trade",
    "EMP_RETAIL"=>"Employees:Retail Trade",
    "EMP_FIN"=>"Employees:Financial Activities",
    "EMP_GOVT"=>"Employees:Government",
    "EMP_PROSERV"=>"Employees:Professional Services",
    "EMP_LEISURE"=>"Employees:Leisure & Hospitality",
    "URATE"=>"Unemployment Rate",
    "U_DURATION"=>"Average Duration of Unemployment (weeks)",
    "U_L5WKS"=>"Unemployment Duration:Persons:Less than 5 Weeks",
    "U_5_14WKS"=>"Unemployment Duration:Persons:5–14 Weeks",
    "U_15_26WKS"=>"Unemployment Duration:Persons:15–26",
    "U_M27WKS"=>"Unemployment Duration:Persons:27 weeks +",
    "HOURS_AVG"=>"Average Weekly Hours:Goods Producing",
    "HOURS_AVG_OT"=>"Average Weekly Overtime Hours:Manufacturing",
    "HSTARTS_NE"=>"Housing Starts:Northeast",
    "HSTARTS_MW"=>"Housing Starts:Midwest",
    "HSTARTS_SOU"=>"Housing Starts:South",
    "HSTARTS_WST"=>"Housing Starts:West",
    "RRRESINV"=>"Real Private Domestic Investment:Residential Quantity ",
    "SFYGM6"=>"Spread of 6 month Tbill to 3 month Tbill",
    "SFYGT1"=>"Spread of 1 year Treasury to 3 month Tbill",
    "SFYGT10"=>"Spread of 10 year Treasury to 3 month Tbill",
    "TOT_RES"=>"Total Reserves of Depository Institutions",
    "TOT_RES_NB"=>"Total Reserves Of Depository Institutions, Nonborrowed",
    "BUS_LOANS"=>"Commercial and Industrial Loans at All Commercial Banks",
    "CONS_CREDIT"=>"Total Nonrevolving Credit Owned and Securitized, Outstanding",
    "SP500"=>"S&P 500 Stock Price Index",
    # "DJIA"=>"Dow Jones Industrial Average",
    "DJIA"=>"Dow-Jones Industrial Stock Price Index for United States",
    "EXR_US"=>"Trade Weighted U.S. Dollar Index: Major Currencies",
    "EXR_SW"=>"Switzerland / U.S. Foreign Exchange Rate",
    "EXR_JAN"=>"Japan / U.S. Foreign Exchange Rate",
    "EXR_UK"=>"U.S. / U.K. Foreign Exchange Rate",
    "EXR_CAN"=>"Canada / U.S. Foreign Exchange Rate",
    "NAPMI"=>"Purchasing Managers Index",
    "NAPM_NEW_ORDERS"=>"NAPM New Orders Index",
    "NAPM_SUP_DEL"=>"NAPM Supplier Deliveries",
    "NAPM_INVENTORIES"=>"NAPM Inventories Index",
    "RNONRESINV"=>"Real private fixed investment: Nonresidential quantity index",
    "RAHE_CONST"=>"",
    "RAHE_MFG"=>"Real Avg. Hourly wages:manufacturing (Deflated w/GDP ",
    "RCOMP_HR"=>"Real Compensation Per Hour (index)",
    "ULC"=>"Unit Labor Cost (index)",
    "CPI_CORE"=>"CPI:Less food and energy",
    "PCED_DUR"=>"PCE:Durable goods price index",
    "PCED_NDUR"=>"PCE:NonDurable goods price index",
    "PCED_SERV"=>"PCE:Services price index",
    "PINV_GDP"=>"Gross private domestic investment price index",
    "PINV_NRES_STRUCT"=>"GPDI:price index:structures",
    "PINV_NRES_EQP"=>"GPDI:price index:Equiptment and software",
    "PINV_RES"=>"GPDI:price index:Residential",
    "PEXPORTS"=>"GDP:Exports Price Index",
    "PIMPORTS"=>"GDP:Imports Price Index",
    "PGOV"=>"Government Consumption and gross investment price index",
    "P_COM"=>"PPI:All commodities price index",
    "P_OIL"=>"PPI:Crude (Divided by PCE Core)",
    "UTL11"=>"Capacity Utilization-Manufacturing",
    "LABOR_PROD"=>"Output per hour all persons:business sector index",
    "UMICH_CONS"=>"University of Michigan Consumer Expectations",
    "M_1"=>"M1 Money stock",
    "M_2"=>"M2 Money stock",
    # new series to make up for missing ones
    "NASDAQ"=>"NASDAQ Composite Index",
    "WTI"=>"Crude Oil Prices: West Texas Intermediate (WTI) - Cushing, Oklahoma", # replacing P_OIL
    "IP"=>"Industrial Production: Total Index", # replacing PMP
    "MAN_NEW_ORDERS"=>"Manufacturers' New Orders: Total Manufacturing", # replacing NAPM_NEW_ORDERS
    "INVENT"=>"Total Business Inventories", # replacing NAPM_INVENTORIES
)
# 0 Demeaned
# 1 Log() and demeaned
# 2 Linear detrended Log() per capita 
# 3 Log() differenced and demeaned 
# 4 Detrended Log()
# 5 Detrended per capita level
fred_transformations = Dict{String,Int64}(
    "RGDP"=>2, #OBS! use 8 for HS model
    "IP_TOTAL"=>2,
    "RGDI"=>2,
    "PGDP"=>3,
    "PCED"=>3,
    "CPI_ALL"=>3,
    "RCONS"=>2,
    "RINV"=>2,
    "RWAGE"=>4,
    "HOURS"=>2,
    "EMP_CES"=>2,
    "EMP_CPS"=>2,
    "FedFunds"=>7, #OBS! changed from 0
    "Tbill_3m"=>0, #OBS! changed from 0
    "AAABond"=>0, #OBS! changed from 0
    "SFYBAAC"=>0,
    "SFYAAAC"=>0,
    "IP_FINAL"=>2,
    "IP_CONS_DBLE"=>2,
    "IP_CONS_NONDBLE"=>2,
    "IP_BUS_EQPT"=>2,
    "IP_DRBLE_MATS"=>2,
    "IP_NONDRBLE_MATS"=>2,
    "IP_MFG"=>2,
    "IP_FUELS"=>2,
    "PMP"=>0,
    "RCONS_DRBLE"=>2,
    "RCONS_NONDRBLE"=>2,
    "RCONS_SERV"=>2,
    "REXPORTS"=>2,
    "RIMPORTS"=>2,
    "RGOV"=>2,
    "EMP_Mining"=>2,
    "EMP_CONST"=>2,
    "EMP_MFG"=>2,
    "EMP_SERVICES"=>2,
    "EMP_TTU"=>2,
    "EMP_WHOLESALE"=>2,
    "EMP_RETAIL"=>2,
    "EMP_FIN"=>2,
    "EMP_GOVT"=>2,
    "EMP_PROSERV"=>2,
    "EMP_LEISURE"=>2,
    "URATE"=>0,
    "U_DURATION"=>0,
    "U_L5WKS"=>2,
    "U_5_14WKS"=>2,
    "U_15_26WKS"=>2,
    "U_M27WKS"=>2,
    "HOURS_AVG"=>0,
    "HOURS_AVG_OT"=>0,
    "HSTARTS_NE"=>1,
    "HSTARTS_MW"=>1,
    "HSTARTS_SOU"=>1,
    "HSTARTS_WST"=>1,
    "RRRESINV"=>2,
    "SFYGM6"=>0,
    "SFYGT1"=>0,
    "SFYGT10"=>0,
    "TOT_RES"=>2,
    "TOT_RES_NB"=>5,
    "BUS_LOANS"=>2,
    "CONS_CREDIT"=>2,
    "SP500"=>3,
    "NASDAQ"=>3,
    "DJIA"=>3,
    "EXR_US"=>3,
    "EXR_SW"=>3,
    "EXR_JAN"=>3,
    "EXR_UK"=>3,
    "EXR_CAN"=>3,
    "NAPMI"=>0,
    "NAPM_NEW_ORDERS"=>0,
    "NAPM_SUP_DEL"=>0,
    "NAPM_INVENTORIES"=>0,
    "RNONRESINV"=>2,
    "RAHE_CONST"=>4,
    "RAHE_MFG"=>4,
    "RCOMP_HR"=>4,
    "ULC"=>4,
    "CPI_CORE"=>3,
    "PCED_DUR"=>3,
    "PCED_NDUR"=>3,
    "PCED_SERV"=>3,
    "PINV_GDP"=>3,
    "PINV_NRES_STRUCT"=>3,
    "PINV_NRES_EQP"=>3,
    "PINV_RES"=>3,
    "PEXPORTS"=>3,
    "PIMPORTS"=>3,
    "PGOV"=>3,
    "P_COM"=>3,
    "P_OIL"=>3,
    "UTL11"=>0,
    "LABOR_PROD"=>4,
    "UMICH_CONS"=>1,
    "M_1"=>2,
    "M_2"=>2,
    "POP"=>6,
    # new series to make up for missing ones
    "WTI"=>3,
    "IP"=>4, # replacing PMP
    "MAN_NEW_ORDERS"=>4, # replacing NAPM_NEW_ORDERS
    "INVENT"=>4, # replacing NAPM_INVENTORIES
)



fred_frequencies = Dict{String,Any}(
    "RGDP" => Quarterly,
    "PGDP" => Quarterly,
    "RCONS" => Quarterly,
    "RINV" => Quarterly,
    "RWAGE" => Monthly,
    "HOURS" => Quarterly,
    "FedFunds" => Monthly,
    "SFYBAAC" => BDaily,
    "NASDAQ" => BDaily,
    "DJIA" => Monthly,
    "POP" => Monthly,
)

fred_groups = Workspace()
fred_groups.core = Workspace()
fred_groups.core.output = ["RGDP", "IP_TOTAL", "RGDI"]
fred_groups.core.inflation = ["PGDP", "PCED", "CPI_ALL"]
fred_groups.core.consumption = ["RCONS"]
fred_groups.core.investment = ["RINV"]
fred_groups.core.wages = ["RWAGE"]
fred_groups.core.labour = ["HOURS", "EMP_CES", "EMP_CPS"]
fred_groups.core.rates = ["FedFunds", "Tbill_3m", "AAABond"]
fred_groups.core.spread = ["SFYBAAC", "SFYAAAC",]

fred_groups.sat = Workspace()
fred_groups.sat.output = ["IP_FINAL", "IP_CONS_DBLE", "IP_CONS_NONDBLE", "IP_BUS_EQPT", "IP_DRBLE_MATS", "IP_NONDRBLE_MATS", "IP_MFG", "IP_FUELS", "PMP", "RCONS_DRBLE", "RCONS_NONDRBLE", "RCONS_SERV", "REXPORTS", "RIMPORTS", "RGOV"]
fred_groups.sat.labour = ["EMP_Mining", "EMP_CONST", "EMP_MFG", "EMP_SERVICES", "EMP_TTU", "EMP_WHOLESALE", "EMP_RETAIL", "EMP_FIN", "EMP_GOVT", "EMP_PROSERV", "EMP_LEISURE", "URATE", "U_DURATION", "U_L5WKS", "U_5_14WKS", "U_15_26WKS", "U_M27WKS", "HOURS_AVG", "HOURS_AVG_OT"]
fred_groups.sat.housing = ["HSTARTS_NE", "HSTARTS_MW", "HSTARTS_SOU", "HSTARTS_WST", "RRRESINV"]
fred_groups.sat.financial = ["SFYGM6", "SFYGT1", "SFYGT10", "TOT_RES", "TOT_RES_NB", "BUS_LOANS", "CONS_CREDIT", "SP500", "DJIA"]
fred_groups.sat.fx = ["EXR_US", "EXR_SW", "EXR_JAP", "EXR_UK", "EXR_CAN"]
fred_groups.sat.investment = ["NAMPI", "NAMP_NEW_ORDERS", "NAMP_SUP_DEL", "NAPM_INVENTORIES", "RNONRESINV"]
fred_groups.sat.prices = ["RAHE_CONST", "RAHE_MFG", "RCOMP_HR", "ULC", "CPI_CORE", "PCED_DUR", "PCED_NDUR", "PCED_SERV", "PINV_GDP", "PINV_NRES_STRUCT", "PINV_NRES_EQP", "PINV_RES", "PEXPORTS", "PIMPORTS", "PGOV", "P_COM", "P_OIL", "WTI"]
fred_groups.sat.other = ["UTL1", "LABOR_PROD", "UMICH_CONS", "M_1", "M_2"]

fred_groups_sw = deepcopy(fred_groups)
fred_groups_sw.sat.financial = [fred_groups_sw.core.spread..., fred_groups_sw.sat.financial...]
fred_groups_sw.core.spread = []
# delete!(:spread, fred_groups_sw.core)
fred_groups_swff = fred_groups

data_series_core_map_sw = Dict{Symbol,Symbol}(
    :output => :Y,
    :inflation => :π,
    :consumption => :C,
    :investment => :I,
    :wages => :w,
    :labour => :L,
    :rates => :R,
)

data_series_core_map_HS = Dict{Symbol,Symbol}(
    :output => :y_obs,
    :inflation => :π_obs,
    :rates => :R_obs,
)

Λ_constraints_HS = Dict{String, Pair{Symbol,Symbol}}()
for group in keys(fred_groups)
    for core_group in keys(fred_groups.core)
        for (i, var) ∈ enumerate(fred_groups.core[core_group])
            if core_group ∈ keys(data_series_core_map_HS)
                if i == 1
                    Λ_constraints_HS[var] = data_series_core_map_HS[core_group] => :one
                else
                    Λ_constraints_HS[var] = data_series_core_map_HS[core_group] => :any
                end
            else
                # skip, no constraints
            end
        end
    end
end
Λ_constraints_HS["y_obs"] = :y_obs=>:one
Λ_constraints_HS["π_obs"] = :π_obs=>:one
Λ_constraints_HS["R_obs"] = :R_obs=>:one
Λ_constraints_HS


Λ_constraints_SW = Dict{String, Pair{Symbol,Symbol}}()
for group in keys(fred_groups)
    for core_group in keys(fred_groups.core)
        for (i, var) ∈ enumerate(fred_groups.core[core_group])
            if core_group ∈ keys(data_series_core_map_sw)
                if i == 1
                    Λ_constraints_SW[var] = data_series_core_map_sw[core_group] => :one
                    if var == "FedFunds"
                        Λ_constraints_SW[var] = data_series_core_map_sw[core_group] => :four
                    end
                else
                    Λ_constraints_SW[var] = data_series_core_map_sw[core_group] => :any
                    if var ∈ ("Tbill_3m", "AAABond")
                        Λ_constraints_SW[var] = data_series_core_map_sw[core_group] => :anyfour
                    end
                end
            else
                # skip, no constraints
            end
        end
    end
end

data_series_core_map_swff = Dict{Symbol,Symbol}(
    :output => :Y,
    :inflation => :π,
    :consumption => :C,
    :investment => :I,
    :wages => :w,
    :labour => :L,
    :rates => :R,
    :spread => :S,
)

Λ_constraints_SWFF = Dict{String, Pair{Symbol,Symbol}}()
for group in keys(fred_groups)
    for core_group in keys(fred_groups.core)
        for (i, var) ∈ enumerate(fred_groups.core[core_group])
            if core_group ∈ keys(data_series_core_map_swff)
                if i == 1
                    Λ_constraints_SWFF[var] = data_series_core_map_swff[core_group] => :one
                    if var == "FedFunds"
                        Λ_constraints_SWFF[var] = data_series_core_map_swff[core_group] => :four
                    end
                else
                    Λ_constraints_SWFF[var] = data_series_core_map_swff[core_group] => :any
                    if var ∈ ("Tbill_3m", "AAABond")
                        Λ_constraints_SWFF[var] = data_series_core_map_swff[core_group] => :anyfour
                    end
                end
            else
                # skip, no constraints
            end
        end
    end
end