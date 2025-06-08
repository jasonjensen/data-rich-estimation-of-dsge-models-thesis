# This file contains code used to download data from FRED
# issues may be encountered when running this code.

using TimeSeriesEcon
using FredData
using JLD2

fred = Fred()


#########################################################
# SECTION 1: GET data
#########################################################
fred_data = Dict{String,FredSeries}()
for (key,val) in pairs(fred_codes)
    try
        if contains(val, '/')
            val1, val2 = split(val, '/')
            fred1 = get_data(fred, val1)
            fred2 = get_data(fred, val2)
            fred_data[key] = Workspace(
                op => '/',
                arg1 => fred1,
                arg2 => fred2
            )
        elseif contains(val, '-')
            val1, val2 = split(val, '-')
            fred1 = get_data(fred, val1)
            fred2 = get_data(fred, val2)
            fred_data[key] = Workspace(
                op => '-',
                arg1 => fred1,
                arg2 => fred2
            )
        elseif contains(val, '+')
            val1, val2 = split(val, '+')
            fred1 = get_data(fred, val1)
            fred2 = get_data(fred, val2)
            fred_data[key] = Workspace(
                op => '+',
                arg1 => fred1,
                arg2 => fred2
            )
        else
            fred_data[key] = get_data(fred, val)
        end
    catch err
        println("ISSUE WITH: $val")
    end
end
get_data(fred, "SP500")
fred_data["POP"] = get_data(fred, "CNP16OV")
fred_data["SFYBAAC"] = get_data(fred, fred_codes["SFYBAAC"])

# can use NASDAQ instead of SP500
NASDAQCOM
fred_data["SP500"] = get_data(fred, fred_codes["SP500"])
fred_data["NASDAQ"] = get_data(fred, fred_codes["NASDAQ"])
fred_data["DJIA"] = get_data(fred, fred_codes["DJIA"])
fred_data["IP"] = get_data(fred, fred_codes["IP"])
fred_data["MAN_NEW_ORDERS"] = get_data(fred, fred_codes["MAN_NEW_ORDERS"])
fred_data["INVENT"] = get_data(fred, fred_codes["INVENT"])
fred_data["DCOILWTICO"] = get_data(fred, "DCOILWTICO")


# corrected series (Jun 1st, 2025)
fred_data["RCONS_DRBLE"] = get_data(fred, fred_codes["RCONS_DRBLE"])
fred_data["RCONS_NONDRBLE"] = get_data(fred, fred_codes["RCONS_NONDRBLE"])
fred_data["RCONS_SERV"] = get_data(fred, fred_codes["RCONS_SERV"])
fred_data["RAHE_CONST"] = get_data(fred, fred_codes["RAHE_CONST"])
# can use M1109BUSM293NNBR for DJI

# key = "SFYGM6"
# val = "TB6MS-TB3MS"
# val1, val2 = split(val, '-')


# fred1 = get_data(fred, val1)
# fred1 = get_data(fred, val; realtime_start=string(Dates.today()), realtime_end=realtime_start=string(Dates.today()))

# fred2 = get_data(fred, val2)
# series = []
#     op => '/',
#     arg1 => fred1,
#     arg2 => fred2
# )

# dt = string(Dates.today())


# fred_data[val1] = fred1
# fred_data[val2] = fred2
# TODO: PPICEM/PCEPILFE
# TODO: GS1-TB3MS, GS10-TB3MS, AAA-GS10, PAYEMS+USGOVT, BAA-GS10

# Institute for Supply Management Data To Be Removed from FRED
# NAPM, NAPMNOI, NAPMII, NAPMPI, 

# Just missing:
# GDPI, IPS299, DNDG, RA3Q086SBEA, CE160V, 

# start_date issue?
# SP500, DJIA

# # issues with: NAPM, GS1-TB3MS,  PPICEM/PCEPILFE, NAPMNOI, GDPI
# # IPS299,  GS10-TB3MS, DNDG RA3Q086SBEA, AAA-GS10, NAPMII, PAYEMS+USGOVT
# # SP500, CE160V, BAA-GS10, MAPMSDI, DJIA, NAPMPI

save("data/input/fred_data.jld2", fred_data)

fred_data = load("data/input/fred_data.jld2")

regular_estimation_series = ["RGDP", "PGDP", "RCONS", "RINV", "RWAGE", "HOURS", "FedFunds", "SFYBAAC", "CNP16OV"]

for key ∈ regular_estimation_series
    @assert key ∈ keys(fred_data) "$key is not in saved data"
end

for key ∈ keys(fred_codes)
    if !haskey(fred_data, key)
        println(key)
    end
end