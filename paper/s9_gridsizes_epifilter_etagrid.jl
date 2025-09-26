
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

# Specify default grids
Rgrid = LinRange(0.01, 10, 1000)
ηgrid = LinRange(0.001, 1, 1000)
Cgrid = 0:Int(10*maximum(Ct))

# Setup empty plots
plotEtaOverTime = plot(label="EpiFilter", xlabel="Day", ylabel="η")
plotEtaFinal = plot(label="EpiFilter", xlabel="η", ylabel=" Density*")
plotRcond = plot(label="EpiFilter", xlabel="Day", ylabel="Rt\n(EpiFilter, default)", ylims=(0, 1.7))
plotRmarg = plot(label="EpiFilter", xlabel="Day", ylabel="Rt\n(EpiFilter, marginalised)", ylims=(0, 1.7))
plotCcond = plot(label="EpiFilter", xlabel="Day", ylabel="Cases\n(EpiFilter, conditional)", ylims=(0, 500))
plotCmarg = plot(label="EpiFilter", xlabel="Day", ylabel="Cases\n(EpiFilter, marginalised)", ylims=(0, 500))


# Run EpiFilter at a coarser grid
ηgrid = LinRange(0.001, 1, 500)
pη0 = ones(length(ηgrid))/length(ηgrid)

(pR, pObs, pRup) = EpiFilterForwards(0.1, w, Ct, Rgrid)
(m, med, l, u) = calculateResults(pR, Rgrid)
plotRcond = plot!(plotRcond, (windin+1):100, m[(windin+1):100], ribbon=(m[(windin+1):100]-l[(windin+1):100], u[(windin+1):100]-m[(windin+1):100]), fillalpha=0.2, label="Coarser grid")

pC = EpiFilterConditionalPredictive(0.1, w, Ct, Rgrid, Cgrid)
(mC, medC, lC, uC) = calculateResults(pC, Cgrid)
plotCcond = plot!(plotCcond, (windin+1):100, mC[(windin+1):100], ribbon=(mC[(windin+1):100]-lC[(windin+1):100], uC[(windin+1):100]-mC[(windin+1):100]), fillalpha=0.2, label="Coarser grid")

(pη, pRgivenη, pRupgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=windin)
(mη, medη, lη, uη) = calculateResults(pη, ηgrid)
plotEtaOverTime = plot!(plotEtaOverTime, 2:100, mη[2:100], ribbon=(mη[2:100]-lη[2:100], uη[2:100]-mη[2:100]), fillalpha=0.2, label="Coarser grid")
plotEtaFinal = plot!(plotEtaFinal, ηgrid[1:findfirst(ηgrid .>= 0.5)], 0.5 * pη[1:findfirst(ηgrid .>= 0.5),end], fillalpha=0.2, label="Coarser grid")

pRmarginalised = EpiFilterMarginalPosterior(pη, pRgivenη)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plotRmarg = plot!(plotRmarg, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Coarser grid")

pCmarginalised = EpiFilterMarginalPredictive(pη, pRupgivenη, w, Ct, Rgrid, Cgrid)
(mC, medC, lC, uC) = calculateResults(pCmarginalised, Cgrid)
plotCmarg = plot!(plotCmarg, (windin+1):100, mC[(windin+1):100], ribbon=(mC[(windin+1):100]-lC[(windin+1):100], uC[(windin+1):100]-mC[(windin+1):100]), fillalpha=0.2, label="Coarser grid")


# Run EpiFilter at default grid
ηgrid = LinRange(0.001, 1, 1000)
pη0 = ones(length(ηgrid))/length(ηgrid)

(pR, pObs, pRup) = EpiFilterForwards(0.1, w, Ct, Rgrid)
(m, med, l, u) = calculateResults(pR, Rgrid)
plotRcond = plot!(plotRcond, (windin+1):100, m[(windin+1):100], ribbon=(m[(windin+1):100]-l[(windin+1):100], u[(windin+1):100]-m[(windin+1):100]), fillalpha=0.2, label="Default grid")

pC = EpiFilterConditionalPredictive(0.1, w, Ct, Rgrid, Cgrid)
(mC, medC, lC, uC) = calculateResults(pC, Cgrid)
plotCcond = plot!(plotCcond, (windin+1):100, mC[(windin+1):100], ribbon=(mC[(windin+1):100]-lC[(windin+1):100], uC[(windin+1):100]-mC[(windin+1):100]), fillalpha=0.2, label="Default grid")

(pη, pRgivenη, pRupgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=windin)
(mη, medη, lη, uη) = calculateResults(pη, ηgrid)
plotEtaOverTime = plot!(plotEtaOverTime, 2:100, mη[2:100], ribbon=(mη[2:100]-lη[2:100], uη[2:100]-mη[2:100]), fillalpha=0.2, label="Default grid")
plotEtaFinal = plot!(plotEtaFinal, ηgrid[1:findfirst(ηgrid .>= 0.5)], pη[1:findfirst(ηgrid .>= 0.5),end], fillalpha=0.2, label="Default grid")

pRmarginalised = EpiFilterMarginalPosterior(pη, pRgivenη)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plotRmarg = plot!(plotRmarg, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Default grid")

pCmarginalised = EpiFilterMarginalPredictive(pη, pRupgivenη, w, Ct, Rgrid, Cgrid)
(mC, medC, lC, uC) = calculateResults(pCmarginalised, Cgrid)
plotCmarg = plot!(plotCmarg, (windin+1):100, mC[(windin+1):100], ribbon=(mC[(windin+1):100]-lC[(windin+1):100], uC[(windin+1):100]-mC[(windin+1):100]), fillalpha=0.2, label="Default grid")


# Run EpiFilter at finer grid
ηgrid = LinRange(0.001, 1, 2000)
pη0 = ones(length(ηgrid))/length(ηgrid)

(pR, pObs, pRup) = EpiFilterForwards(0.1, w, Ct, Rgrid)
(m, med, l, u) = calculateResults(pR, Rgrid)
plotRcond = plot!(plotRcond, (windin+1):100, m[(windin+1):100], ribbon=(m[(windin+1):100]-l[(windin+1):100], u[(windin+1):100]-m[(windin+1):100]), fillalpha=0.2, label="Finer grid")

pC = EpiFilterConditionalPredictive(0.1, w, Ct, Rgrid, Cgrid)
(mC, medC, lC, uC) = calculateResults(pC, Cgrid)
plotCcond = plot!(plotCcond, (windin+1):100, mC[(windin+1):100], ribbon=(mC[(windin+1):100]-lC[(windin+1):100], uC[(windin+1):100]-mC[(windin+1):100]), fillalpha=0.2, label="Finer grid")

(pη, pRgivenη, pRupgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=windin)
(mη, medη, lη, uη) = calculateResults(pη, ηgrid)
plotEtaOverTime = plot!(plotEtaOverTime, 2:100, mη[2:100], ribbon=(mη[2:100]-lη[2:100], uη[2:100]-mη[2:100]), fillalpha=0.2, label="Finer grid")
plotEtaFinal = plot!(plotEtaFinal, ηgrid[1:findfirst(ηgrid .>= 0.5)], 2 * pη[1:findfirst(ηgrid .>= 0.5),end], fillalpha=0.2, label="Finer grid")

pRmarginalised = EpiFilterMarginalPosterior(pη, pRgivenη)
(mR, medR, lR, uR) = calculateResults(pRmarginalised, Rgrid)
plotRmarg = plot!(plotRmarg, (windin+1):100, mR[(windin+1):100], ribbon=(mR[(windin+1):100]-lR[(windin+1):100], uR[(windin+1):100]-mR[(windin+1):100]), fillalpha=0.2, label="Finer grid")

pCmarginalised = EpiFilterMarginalPredictive(pη, pRupgivenη, w, Ct, Rgrid, Cgrid)
(mC, medC, lC, uC) = calculateResults(pCmarginalised, Cgrid)
plotCmarg = plot!(plotCmarg, (windin+1):100, mC[(windin+1):100], ribbon=(mC[(windin+1):100]-lC[(windin+1):100], uC[(windin+1):100]-mC[(windin+1):100]), fillalpha=0.2, label="Finer grid")


# Make the giant plot
plt = plot(plotEtaOverTime, plotEtaFinal, plotRcond, plotRmarg, plotCcond, plotCmarg, layout=(3, 2), size=(800, 700), dpi=300, margin=3mm)
# png("paper/figures/s9_gridsize_epifilter_etagrid.png")
savefig(plt, "paper/figures/s9_gridsize_epifilter_etagrid.pdf")