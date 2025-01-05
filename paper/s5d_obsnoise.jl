
using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Set options
Random.seed!(42)
p = 0.5

# Save outputs
out_all = DataFrame()
ee_parampost_all = DataFrame()
ef_parampost_all= DataFrame()

println("Starting random walk example...")
(Ct, Rt) = simulateRandomWalkRt(; minCt=5, maxCt=500)
Ct = simulateBetaBinomialObservationNoise(Int.(round.(Ct/p)), p = p)
(out_tmp_rw, ee_parampost_tmp_rw, ef_parampost_tmp_rw) = fitModels(Ct, Rt, "Random walk")
out_all = vcat(out_all, out_tmp_rw)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_rw)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_rw)

println("Starting sinusoidal example...")
(Ct, Rt) = simulateSinusoidalRt(; ω=21)
Ct = simulateBetaBinomialObservationNoise(Int.(round.(Ct/p)), p = p)
(out_tmp_sin, ee_parampost_tmp_sin, ef_parampost_tmp_sin) = fitModels(Ct, Rt, "Sinusoidal")
out_all = vcat(out_all, out_tmp_sin)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_sin)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_sin)

println("Starting step-change example...")
(Ct, Rt) = simulateStepChangeRt( ; tChange=28)
Ct = simulateBetaBinomialObservationNoise(Int.(round.(Ct/p)), p = p)
(out_tmp_step, ee_parampost_tmp_step, ef_parampost_tmp_step) = fitModels(Ct, Rt, "Step change")
out_all = vcat(out_all, out_tmp_step)
ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp_step)
ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp_step)

CSV.write("paper/outputs/s5d_obsnoise.csv", out_all)
CSV.write("paper/outputs/s5d_obsnoise_epiestimposteriors.csv", ee_parampost_all)
CSV.write("paper/outputs/s5d_obsnoise_epifilterposteriors.csv", ef_parampost_all)