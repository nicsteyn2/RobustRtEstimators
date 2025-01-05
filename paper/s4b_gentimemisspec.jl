
using Plots, Random, CSV, DataFrames

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Pre-allocate outputs
out_all = DataFrame()
ee_parampost_all = DataFrame()
ef_parampost_all= DataFrame()

# Load data
df = CSV.read("data/exampledata.csv", DataFrame)
Ct = df.Cases[df.simulation.=="Random walk"]
Rt = df.TrueRt[df.simulation.=="Random walk"]

# Create a Gamma serial interval distribution with specified mean and std dev
function serial_interval(; mean=6.5, stddev=4.2)

    shape0 = (mean / stddev)^2
    scale0 = (stddev^2) / mean

    return(Gamma(shape0, scale0))
    
end

# # Fit models
w = pdf.(serial_interval(), 1:100)
w = w/sum(w)
plot(w[1:21])
(out, ee_parampost, ef_parampost) = fitModels(Ct, Rt, "Standard", w=w)
out_all = vcat(out_all, out)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost)

w = pdf.(serial_interval(mean=4.2), 1:100)
w = w/sum(w)
plot!(w[1:21])
(out, ee_parampost, ef_parampost) = fitModels(Ct, Rt, "Exponential(4.2)", w=w)
out_all = vcat(out_all, out)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost)

w = pdf.(serial_interval(mean=6.5*1.5), 1:100)
w = w/sum(w)
plot!(w[1:21])
(out, ee_parampost, ef_parampost) = fitModels(Ct, Rt, "Mean x1.5", w=w)
out_all = vcat(out_all, out)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost)

CSV.write("paper/outputs/s4b_SImisspec.csv", out_all)
CSV.write("paper/outputs/s4b_SImisspec_epiestimposteriors.csv", ee_parampost_all)
CSV.write("paper/outputs/s4b_SImisspec_epifilterposteriors.csv", ef_parampost_all)