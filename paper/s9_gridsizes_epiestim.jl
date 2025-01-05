
using Plots, Random, Measures

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Generate simulated data
Random.seed!(42)
w = defaultGenTimeDist()
(Ct, Rt) = simulateRandomWalkRt(; w=w, η=0.1, C0=100, minCt=5, maxCt=500)
windin = 10


# First start by examining impact on EpiEstim marginal posterior on Rt (this is the only time that a grid impacts the results of EpiEstim)
plotR = plot(label="EpiEstim", xlabel="Day", ylabel="Rt (EpiEstim)", ylims=(0, 1.7))

Rgrid = LinRange(0.01, 10, 500)
pRmarginalised = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plot!(plotR, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Estimated Rt (coarser grid)")

Rgrid = LinRange(0.01, 10, 1000)
pRmarginalised = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plot!(plotR, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Estimated Rt (standard grid)")

Rgrid = LinRange(0.01, 10, 2000)
pRmarginalised = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plot!(plotR, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Estimated Rt (finer grid)")

plot!(plotR, 1:100, Rt, label="True Rt", linestyle=:dash, color=:black)

# And repeat on a much coarser grid
plotR2 = plot(label="EpiEstim", xlabel="Day", ylabel="Rt (EpiEstim)", ylims=(0, 1.7))

Rgrid = LinRange(0.01, 10, 100)
pRmarginalised = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plot!(plotR2, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Estimated Rt (extremely coarse grid)")

Rgrid = LinRange(0.01, 10, 1000)
pRmarginalised = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plot!(plotR2, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Estimated Rt (standard grid)")

plot!(plotR2, 1:100, Rt, label="True Rt", linestyle=:dash, color=:black)

plt = plot(plotR, plotR2, layout=(2, 1), size=(800, 500), dpi=300, margin=3mm)
png("paper/figures/s9_gridsize_epiestim.png")