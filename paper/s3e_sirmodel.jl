

using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Save outputs
out_all = DataFrame()
ee_parampost_all = DataFrame()
ef_parampost_all= DataFrame()

# Specify the generation time distribution
γ = 1/6.5
w = γ * (1 - γ).^(0:99)
w = w/sum(w)

# Simulate and fit basic model
(Ct, Rt) = simulateStochasticSIR()
(out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct, Rt, "SIR"; w=w)
out_all = vcat(out_all, out_tmp)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)

# Simulate and fit model with extra noise
ρ = 0.5
Ct2 = simulateBinomialObservationNoise(Int.(round.(Ct/ρ)), p=ρ)
(out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct2, Rt, "SIR (with noise)"; w=w)
out_all = vcat(out_all, out_tmp)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)

CSV.write("paper/outputs/s3e_sirmodel.csv", out_all)
CSV.write("paper/outputs/s3e_sirmodel_epiestimposteriors.csv", ee_parampost_all)
CSV.write("paper/outputs/s3e_sirmodel_epifilterposteriors.csv", ef_parampost_all)