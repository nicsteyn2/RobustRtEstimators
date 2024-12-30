

using Plots, Random, CSV, DataFrames

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Load data
df = CSV.read("data/exampledata.csv", DataFrame)
Ct = df.Cases[df.simulation.=="Random walk"]
Rt = df.TrueRt[df.simulation.=="Random walk"]

# Sample a realisation of an uncertain serial interval
function serial_interval()

    x = sqrt(rand(Normal(1, 0.5/3)))
    shape1 = 2.3669 * x
    scale1 = 2.7463 * x
    return(Gamma(shape1, scale1))
    
end


# N = 1000
# X = zeros(N,100)
# for ii = 1:N
#     X[ii,:] = pdf.(serial_interval(), 1:100)
# end

# X = X ./ sum(X, dims=2)

# x = mean(X,dims=1)'
# sum(x' .* collect(1:100))

# sh = (6.5 / 4.2)^2
# sc = (4.2^2) / 6.5
# x2 = pdf.(Gamma(sh, sc), 1:100)

# plot(x[1:21])
# display(plot!(x2[1:21]))
# display(plot!(pdf.(serial_interval(), 1:21), color=:black, label=false))


# Ffit EpiFilter with an uncertain serial interval
function FitEpiFilter_SIUncertainty(SI, Ct, Rt, simulationName; N=50, windin=10, showProgress=true)
    
    # Fit on a finer grid to speed things up
    Rgrid = LinRange(0.01, 10, 400)
    ηgrid = LinRange(0.001, 1, 400)  
    Cgrid = 0:(10*maximum(Ct))
    pη0 = ones(length(ηgrid))/length(ηgrid)
    
    # Sample a serial interval
    w = pdf.(SI(), 1:100)
    w = w/sum(w)
    
    # Fit EpiFilter N times
    pη = zeros(length(ηgrid), length(Ct), N)
    pR = zeros(length(Rgrid), length(Ct), N)
    pC = zeros(length(Cgrid), length(Ct), N)
    ProgBar = Progress(N+1, dt=1, barlen=50, desc="Fitting EpiFilter at sampled serial intervals...", enabled=showProgress)
    Threads.@threads for ii = 1:N
        (pη_new, pRgivenη, pRupgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=windin, showProgress=false)
        pη[:,:,ii] = pη_new
        pR[:,:,ii] = EpiFilterMarginalPosterior(pη_new, pRgivenη)
        pC[:,:,ii] = EpiFilterMarginalPredictive(pη_new, pRupgivenη, w, Ct, Rgrid, Cgrid; windin=windin)
        next!(ProgBar)
    end
    pη = dropdims(sum(pη, dims=3), dims=3)/N
    pR = dropdims(sum(pR, dims=3), dims=3)/N
    pC = dropdims(sum(pC, dims=3), dims=3)/N

    # Calculate results
    (mη, _, lη, uη) = calculateResults(pη, ηgrid)
    (mR, _, lR, uR) = calculateResults(pR, Rgrid)
    (mC, _, lC, uC) = calculateResults(pC, Cgrid)
    score = CRPS(cumsum(pC, dims=1), Cgrid, Ct; stepwise=true)

    out = DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mR, lower=lR, upper=uR, mean_param=mη, lower_param=lη, upper_param=uη, mean_cases=mC, lower_cases=lC, upper_cases=uC, crps=score, simulation=simulationName, model="EpiFilter", fit="SI Uncertainty")
    paramposterior = DataFrame(eta=ηgrid, peta=pη[:,end], simulation=simulationName, model="EpiFilter", fit="SI Uncertainty")

    # Fit EpiFilter at mean SI
    w = pdf.(Gamma(2.3669, 2.7463), 1:100)
    w = w/sum(w)
    (pη, pRgivenη, pRupgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=windin, showProgress=false)
    pR = EpiFilterMarginalPosterior(pη, pRgivenη)
    pC = EpiFilterMarginalPredictive(pη, pRupgivenη, w, Ct, Rgrid, Cgrid; windin=windin)
    next!(ProgBar)
    (mη, _, lη, uη) = calculateResults(pη, ηgrid)
    (mR, _, lR, uR) = calculateResults(pR, Rgrid)
    (mC, _, lC, uC) = calculateResults(pC, Cgrid)
    score = CRPS(cumsum(pC, dims=1), Cgrid, Ct; stepwise=true)

    out = vcat(out, DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mR, lower=lR, upper=uR, mean_param=mη, lower_param=lη, upper_param=uη, mean_cases=mC, lower_cases=lC, upper_cases=uC, crps=score, simulation=simulationName, model="EpiFilter", fit="No Uncertainty"))
    paramposterior = vcat(paramposterior, DataFrame(eta=ηgrid, peta=pη[:,end], simulation=simulationName, model="EpiFilter", fit="No Uncertainty"))

    return(out, paramposterior)
    
end

function FitEpiEstim_SIUncertainty(SI, Ct, Rt, simulationName; N=50, windin=10, showProgress=true)

    # Set grids and options
    kvals = collect(1:30)
    Rgrid = LinRange(0.01, 10, 500)
    Cgrid = 0:(20*maximum(Ct))

    # Sample a serial interval
    w = pdf.(SI(), 1:100)
    w = w/sum(w)
    
    # Fit EpiFilter N times
    pK = zeros(length(kvals), length(Ct))
    pR = zeros(length(Rgrid), length(Ct))
    pC = zeros(length(Cgrid), length(Ct))
    ProgBar = Progress(N+1, dt=1, barlen=50, desc="Fitting EpiEstim at sampled serial intervals...", enabled=showProgress)
    for ii = 1:N
        (loglik, _) = EpiEstimLogLik(kvals, w, Ct; windin=windin)
        pK_unnorm = exp.(loglik .- maximum(loglik, dims=1))
        pK = pK + (pK_unnorm ./ sum(pK_unnorm, dims=1))
        pR = pR + EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
        pC = pC + EpiEstimMarginalPredictive(w, Ct, Cgrid; windin=windin)
        next!(ProgBar)
    end
    pK = pK/N
    pR = pR/N
    pC = pC/N

    # Calculate results
    (mk, _, lk, uk) = calculateResults(pK, kvals)
    (mR, _, lR, uR) = calculateResults(pR, Rgrid)
    (mC, _, lC, uC) = calculateResults(pC, Cgrid)
    score = CRPS(cumsum(pC, dims=1), Cgrid, Ct; stepwise=true)

    out = DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mR, lower=lR, upper=uR, mean_param=mk, lower_param=lk, upper_param=uk, mean_cases=mC, lower_cases=lC, upper_cases=uC, crps=score, simulation=simulationName, model="EpiEstim", fit="SI Uncertainty")
    paramposterior = DataFrame(k=kvals, pk=pK[:,end], simulation=simulationName, model="EpiEstim", fit="SI Uncertainty")

    # Fit EpiEstim at mean SI
    w = pdf.(Gamma(2.3669, 2.7463), 1:100)
    w = w/sum(w)
    (loglik, _) = EpiEstimLogLik(kvals, w, Ct; windin=windin)
    pK = exp.(loglik .- maximum(loglik, dims=1))
    pK = pK./sum(pK, dims=1)
    pR = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin)
    pC = EpiEstimMarginalPredictive(w, Ct, Cgrid; windin=windin)
    (mk, _, lk, uk) = calculateResults(pK, kvals)
    (mR, _, lR, uR) = calculateResults(pR, Rgrid)
    (mC, _, lC, uC) = calculateResults(pC, Cgrid)
    score = CRPS(cumsum(pC, dims=1), Cgrid, Ct; stepwise=true)
    next!(ProgBar)

    out = vcat(out, DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mR, lower=lR, upper=uR, mean_param=mk, lower_param=lk, upper_param=uk, mean_cases=mC, lower_cases=lC, upper_cases=uC, crps=score, simulation=simulationName, model="EpiEstim", fit="No Uncertainty"))
    paramposterior = vcat(paramposterior, DataFrame(k=kvals, pk=pK[:,end], simulation=simulationName, model="EpiEstim", fit="No Uncertainty"))

    return(out, paramposterior)

end


(EpiFilterSI, EpiFilterSI_params) = FitEpiFilter_SIUncertainty(serial_interval, Ct, Rt, "Random walk"; N=50)
CSV.write("paper/outputs/s9epifilter_siuncertainty.csv", EpiFilterSI)
CSV.write("paper/outputs/s9epifilter_siuncertainty_params.csv", EpiFilterSI_params)

(EpiEstimSI, EpiEstimSI_params) = FitEpiEstim_SIUncertainty(serial_interval, Ct, Rt, "Random walk"; N=50)
CSV.write("paper/outputs/s9epiestim_siuncertainty.csv", EpiEstimSI)
CSV.write("paper/outputs/s9epiestim_siuncertainty_params.csv", EpiEstimSI_params)




