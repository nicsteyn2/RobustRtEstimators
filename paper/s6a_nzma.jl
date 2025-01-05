

using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Pre-allocate outputs
out_all = DataFrame()
ee_parampost_all = DataFrame()
ef_parampost_all= DataFrame()

# Load data manually
nzdata = CSV.read("data/nzcovid_moh.csv", DataFrame)
sort!(nzdata, :date)
Ctraw = nzdata.local[539:720]

# Apply 5-day and 10-day smoothers
Ct5 = Int.(round.([mean(Ctraw[max(1, ii-2):min(length(Ctraw), ii+2)]) for ii in 1:length(Ctraw)]))
Ct10 = Int.(round.([mean(Ctraw[max(1, ii-4):min(length(Ctraw), ii+4)]) for ii in 1:length(Ctraw)]))

# Fetch standard generation time distribution
w = pdf.(Weibull(2.826, 5.664), 1:length(Ctraw))
w = w/sum(w)
Rt = similar(Ctraw) * NaN # As we don't know what the true Rt is

# Fit model to raw data
(out_tmp_rw, ee_parampost_tmp_rw, ef_parampost_tmp_rw) = fitModels(Ctraw, Rt, "Raw"; w=w, windin=3)
out_all = vcat(out_all, out_tmp_rw)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_rw)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_rw)

# Fit model to 5-day MA
(out_tmp_rw, ee_parampost_tmp_rw, ef_parampost_tmp_rw) = fitModels(Ct5, Rt, "5-day"; w=w, windin=3)
out_all = vcat(out_all, out_tmp_rw)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_rw)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_rw)

# Fit model to 10-day MA
(out_tmp_rw, ee_parampost_tmp_rw, ef_parampost_tmp_rw) = fitModels(Ct10, Rt, "10-day"; w=w, windin=3)
out_all = vcat(out_all, out_tmp_rw)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_rw)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_rw)

rename!(out_all, :simulation => :movingaverage)
rename!(ee_parampost_all, :simulation => :movingaverage)
rename!(ef_parampost_all, :simulation => :movingaverage)

CSV.write("paper/outputs/s6a_nzmovingaverage.csv", out_all)
CSV.write("paper/outputs/s6a_nzmovingaverage_epiestimposteriors.csv", ee_parampost_all)
CSV.write("paper/outputs/s6a_nzmovingaverage_epifilterposteriors.csv", ef_parampost_all)
