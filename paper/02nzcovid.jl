
using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Save outputs
out_all = DataFrame()
ee_parampost_all = DataFrame()
ef_parampost_all= DataFrame()

(Ct, w) = loadData("NZCOVID_AUG2021")
Rt = similar(Ct) * NaN

# 7-day smooth Ct data
# Ct_smooth = Int.(round.([mean(Ct[max(1, ii-2):min(length(Ct), ii+2)]) for ii in 1:length(Ct)]))
# display(plot(Ct_smooth))

(out_tmp_rw, ee_parampost_tmp_rw, ef_parampost_tmp_rw) = fitModels(Ct_smooth, Rt, "Random walk"; w=w, windin=3)
out_all = vcat(out_all, out_tmp_rw)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_rw)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_rw)

CSV.write("paperv2/outputs/02nzcovid.csv", out_all)
CSV.write("paperv2/outputs/02nzcovid_epiestimposteriors.csv", ee_parampost_all)
CSV.write("paperv2/outputs/02nzcovid_epifilterposteriors.csv", ef_parampost_all)
