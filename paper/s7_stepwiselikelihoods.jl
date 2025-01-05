
using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Set options
Random.seed!(42)
windin=2
w = defaultGenTimeDist()
kvals = collect(1:30)
ηvals = 0.05:0.01:0.3
Rgrid = LinRange(0.05, 10, 1000)

# Save outputs
data_out = DataFrame()
epifilter_out = DataFrame()
epiestim_out = DataFrame()

# Run on random walk example data
println("Starting random walk example...")
(Ct, Rt) = simulateRandomWalkRt(; minCt=5, maxCt=500)

(_, ee_swll) = EpiEstimLogLik(kvals, w, Ct; windin=windin)
(_, ef_swll) = EpiFilterLogLik(ηvals, w, Ct, Rgrid)

data_tmp = DataFrame(t=1:100, Ct=Ct, Rt=Rt, simulation="Random walk")
epiestim_tmp = DataFrame()
for ii = 1:length(kvals)
    global epiestim_tmp = vcat(epiestim_tmp, DataFrame(t=1:100, swll=ee_swll[ii,:], k = kvals[ii], simulation="Random walk"))
end
epifilter_tmp = DataFrame()
for ii = 1:length(ηvals)
    global epifilter_tmp = vcat(epifilter_tmp, DataFrame(t=1:100, swll=ef_swll[ii,:], eta = ηvals[ii], simulation="Random walk"))
end

data_out = vcat(data_out, data_tmp)
epiestim_out = vcat(epiestim_out, epiestim_tmp)
epifilter_out = vcat(epifilter_out, epifilter_tmp)


# Run on sinusoidal example data
println("Starting sinusoidal example...")
(Ct, Rt) = simulateSinusoidalRt(; ω=35)

(_, ee_swll) = EpiEstimLogLik(kvals, w, Ct; windin=windin)
(_, ef_swll) = EpiFilterLogLik(ηvals, w, Ct, Rgrid)

data_tmp = DataFrame(t=1:100, Ct=Ct, Rt=Rt, simulation="Sinusoidal")
epiestim_tmp = DataFrame()
for ii = 1:length(kvals)
    global epiestim_tmp = vcat(epiestim_tmp, DataFrame(t=1:100, swll=ee_swll[ii,:], k = kvals[ii], simulation="Sinusoidal"))
end
epifilter_tmp = DataFrame()
for ii = 1:length(ηvals)
    global epifilter_tmp = vcat(epifilter_tmp, DataFrame(t=1:100, swll=ef_swll[ii,:], eta = ηvals[ii], simulation="Sinusoidal"))
end

data_out = vcat(data_out, data_tmp)
epiestim_out = vcat(epiestim_out, epiestim_tmp)
epifilter_out = vcat(epifilter_out, epifilter_tmp)


# Run on step-change example data
println("Starting step-change example...")
(Ct, Rt) = simulateStepChangeRt( ; tChange=28)

(_, ee_swll) = EpiEstimLogLik(kvals, w, Ct; windin=windin)
(_, ef_swll) = EpiFilterLogLik(ηvals, w, Ct, Rgrid)

data_tmp = DataFrame(t=1:100, Ct=Ct, Rt=Rt, simulation="Step-change")
epiestim_tmp = DataFrame()
for ii = 1:length(kvals)
    global epiestim_tmp = vcat(epiestim_tmp, DataFrame(t=1:100, swll=ee_swll[ii,:], k = kvals[ii], simulation="Step-change"))
end
epifilter_tmp = DataFrame()
for ii = 1:length(ηvals)
    global epifilter_tmp = vcat(epifilter_tmp, DataFrame(t=1:100, swll=ef_swll[ii,:], eta = ηvals[ii], simulation="Step-change"))
end

data_out = vcat(data_out, data_tmp)
epiestim_out = vcat(epiestim_out, epiestim_tmp)
epifilter_out = vcat(epifilter_out, epifilter_tmp)



CSV.write("paper/outputs/s7_swll_data.csv", data_out)
CSV.write("paper/outputs/s7_swll_epiestim.csv", epiestim_out)
CSV.write("paper/outputs/s7_swll_epifilter.csv", epifilter_out)


