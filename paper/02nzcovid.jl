
using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Pre-allocate output dataframes
out_all = DataFrame()
ee_parampost_all = DataFrame()
ef_parampost_all= DataFrame()

# Load data
(Ct, w) = loadData("NZCOVID_AUG2021")
Rt = similar(Ct) * NaN

# Fit models
(out_tmp_rw, ee_parampost_tmp_rw, ef_parampost_tmp_rw) = fitModels(Ct, Rt, "NZ COVID Aug 2021"; w=w, windin=3)
out_all = vcat(out_all, out_tmp_rw)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_rw)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_rw)

# Save outputs
CSV.write("paper/outputs/02nzcovid.csv", out_all)
CSV.write("paper/outputs/02nzcovid_epiestimposteriors.csv", ee_parampost_all)
CSV.write("paper/outputs/02nzcovid_epifilterposteriors.csv", ef_parampost_all)

# Fetch credible intervals on eta
ef_parampost_all.eta[findall(ef_parampost_all.p .== maximum(ef_parampost_all.p))]
ef_parampost_all.eta[findlast(cumsum(ef_parampost_all.p) .<= 0.025)]
ef_parampost_all.eta[findfirst(cumsum(ef_parampost_all.p) .>= 0.975)]