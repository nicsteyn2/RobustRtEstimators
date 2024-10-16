 
 using CSV, DataFrames, Distributions, Dates
 
 
 function calculateResults(X::Vector, Rgrid; α=0.05)
    
    m = sum(X .* Rgrid)
    med = Rgrid[findfirst(cumsum(X) .> 0.5)]
    l = Rgrid[findfirst(cumsum(X) .> α/2)]
    u = Rgrid[findfirst(cumsum(X) .> 1 - α/2)]
    
    return(m, med, l, u)
    
end

function calculateResults(X::Matrix, Rgrid; α=0.05)
    
    if length(size(X)) > 1
        
        nDay = size(X)[2]
        
        valid = vec(sum(X, dims=1) .>= 0.99)
        
        if sum(.!valid) > 0
            println("Warning: $(sum(.!valid)) element(s) of X do not appear to be normalised, or may contain missing/invalid data.")
        end
        
        m = [valid[tt] ? sum(X[:,tt] .* Rgrid) : NaN for tt = 1:nDay]
        med = [valid[tt] ? Rgrid[findfirst(cumsum(X[:,tt]) .> 0.5)] : NaN for tt = 1:nDay]
        l = [valid[tt] ? Rgrid[findfirst(cumsum(X[:,tt]) .> α/2)] : NaN for tt = 1:nDay]
        u = [valid[tt] ? Rgrid[findfirst(cumsum(X[:,tt]) .> 1 - α/2)] : NaN for tt = 1:nDay]
        
        return(m, med, l, u)
        
    else
        
        m = sum(X .* Rgrid)
        med = Rgrid[findfirst(cumsum(X) .> 0.5)]
        l = Rgrid[findfirst(cumsum(X) .> α/2)]
        u = Rgrid[findfirst(cumsum(X) .> 1 - α/2)]
        
        return(m, med, l, u)
        
    end
    
    
end


function calculateResults(X::Vector{<:Distribution}; α=0.05)
    
    m = mean.(X)
    med = quantile.(X, 0.5)
    l = quantile.(X, α/2)
    u = quantile.(X, 1 - α/2)
    return(m, med, l, u)
    
end


function defaultGenTimeDist(tMax=100)
    
    w = pdf.(Gamma(2.3669, 2.7463), 1:tMax)
    return(w/sum(w))
    
end


function loadData(source)    
    
    if source == "Simulated"
        
        data = CSV.read("data/sim.csv", DataFrame)
        w = CSV.read("data/gensim.csv", DataFrame)
        w = w[:, 1]
        
        return(data.C, w)
        
    elseif source == "NZCOVID" # Data from here: https://github.com/minhealthnz/nz-covid-data
        
        nzdata = CSV.read("data/nzcovid_moh.csv", DataFrame)
        sort!(nzdata, :date)

        cases = nzdata.border[1:100] .+ nzdata.local[1:100]
        w = defaultGenTimeDist(length(cases))

        return(cases, w)

        
    elseif source == "NZCOVID_AUG2021" # Data from here: https://github.com/minhealthnz/nz-covid-data
        
        nzdata = CSV.read("data/nzcovid_moh.csv", DataFrame)
        sort!(nzdata, :date)
        
        cases = nzdata.local[539:720]
        w = defaultGenTimeDist(length(cases))

        # Apply 5-day smoother
        cases_smooth = Int.(round.([mean(cases[max(1, ii-2):min(length(cases), ii+2)]) for ii in 1:length(cases)]))
        
        return(cases_smooth, w)
        
    elseif source == "flu"
        
        flu = CSV.read("data/Iflu.csv", DataFrame, header=false)
        flu = flu[:, 1]
        
        genflu = CSV.read("data/genflu.csv", DataFrame, header=false)
        genflu = genflu[2:end, 1]
        
        w = zeros(length(flu))
        w[1:length(genflu)] = genflu
        
        return(flu, w)
        
    elseif source == "flufilt"
        
        flu = CSV.read("data/IfluFilt.csv", DataFrame, header=false)
        flu = flu[:, 1]
        
        genflu = CSV.read("data/genflu.csv", DataFrame, header=false)
        genflu = genflu[2:end, 1]
        
        w = zeros(length(flu))
        w[1:length(genflu)] = genflu
        
        return(flu, w)
        
    elseif source == "sars"
        
        sars = CSV.read("data/Isars.csv", DataFrame, header=false)
        sars = sars[:, 1]
        gensars = CSV.read("data/gensars.csv", DataFrame, header=false)
        gensars = gensars[2:end, 1]
        w = zeros(length(sars))
        w[1:length(gensars)] = gensars
        
        return(sars, w)
        
    elseif source == "sarsfilt"
        
        sars = CSV.read("data/IsarsFilt.csv", DataFrame, header=false)
        sars = sars[:, 1]
        gensars = CSV.read("data/gensars.csv", DataFrame, header=false)
        gensars = gensars[2:end, 1]
        w = zeros(length(sars))
        w[1:length(gensars)] = gensars
        
        return(sars, w)
        
    elseif source == "nz2"
        
        nz2 = CSV.read("data/nz2.csv", DataFrame)
        cases = nz2.L[555:650]
        w = CSV.read("data/gensim.csv", DataFrame)
        w = w[:, 1]
        return(cases, w)
        
    else
        
        @error("Invalid choice of data souce. See support.jl::loadData() for options.")
        
    end
    
end




# This function fits both EpiEstim and EpiFilter models to data and returns a bunch of results.
# While we use this a lot for the paper, we recommend using the example code to start with.
function fitModels(Ct, Rt, simulationName; windin=10, informativePriors=false, kPrior=missing, ηPrior=missing, w=missing, showProgress=true)
    
    out = DataFrame()
    ee_parampost = DataFrame()
    ef_parampost = DataFrame()
    
    if ismissing(w)
        w = defaultGenTimeDist()
    end
    Rgrid = LinRange(0.01, 10, 1000)
    ηgrid = LinRange(0.001, 1, 1000)
    Cgrid = 0:(20*maximum(Ct))
    pη0 = missing
    if !informativePriors
        pη0 = ones(length(ηgrid))/length(ηgrid)
    else
        if ismissing(ηPrior)
            error("You must specify a prior distribution for η if informativePriors = true.")
        end
        pη0 = ηPrior
    end
    
    # Fit EpiEstim (both default and conditional)
    kvals = collect(1:30)
    (loglik, _) = EpiEstimLogLik(kvals, w, Ct; windin=windin)
    loglik = loglik .- maximum(loglik, dims=1)
    pK = missing
    if !informativePriors
        pK = exp.(loglik .- maximum(loglik, dims=1))
        pK = pK./sum(pK, dims=1) # This assumes a Uniform prior distribution on 1:kMax
    else
        if ismissing(kPrior)
            error("You must specify a prior distribution for k if informativePriors = true.")
        end
        pK = exp.(loglik .- maximum(loglik, dims=1)) .* kPrior
        pK = pK ./ sum(pK, dims=1)
    end
    (mk, _, lk, uk) = calculateResults(pK, kvals)
    ee_parampost = vcat(ee_parampost, DataFrame(k=kvals, p=pK[:,end], simulation=simulationName, model="EpiEstim"))
    
    # Marginal results
    pR = EpiEstimMarginalPosterior(w, Ct, Rgrid; windin=windin, kPrior=kPrior)
    (mR, _, lR, uR) = calculateResults(pR, Rgrid)
    
    pC = EpiEstimMarginalPredictive(w, Ct, Cgrid; windin=windin, kPrior=kPrior)
    (mC, _, lC, uC) = calculateResults(pC, Cgrid)
    
    scoreEEMarginal = CRPS(cumsum(pC, dims=1), Cgrid, Ct; stepwise=true)
    
    out = vcat(out, DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mR, lower=lR, upper=uR, mean_param=mk, lower_param=lk, upper_param=uk, mean_cases=mC, lower_cases=lC, upper_cases=uC, crps=scoreEEMarginal, simulation=simulationName, model="EpiEstim", fit="Marginalised"))
    
    # (Conditional results)
    pRcond = EpiEstimConditionalPosterior(7, w, Ct)
    (mRcond, lRcond, uRcond) = [mean.(pRcond), quantile.(pRcond, 0.025), quantile.(pRcond, 0.975)]
    
    pCcond = EpiEstimConditionalPredictive(7, w, Ct)
    (mCcond, lCcond, uCcond) = [mean.(pCcond), quantile.(pCcond, 0.025), quantile.(pCcond, 0.975)]
    
    scoreEECond = CRPS(pCcond, Ct; stepwise=true)
    
    out = vcat(out, DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mRcond, lower=lRcond, upper=uRcond, mean_param=7, lower_param=NaN, upper_param=NaN, mean_cases=mCcond, lower_cases=lCcond, upper_cases=uCcond, crps=scoreEECond, simulation=simulationName, model="EpiEstim", fit="Conditional (k = 7)"))
    
    
    # Fit EpiFilter (both default and conditional)
    (pη, pRgivenη, pRupgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=windin, showProgress=showProgress)
    (mη, _, lη, uη) = calculateResults(pη, ηgrid)
    ef_parampost = vcat(ef_parampost, DataFrame(eta=ηgrid, p=pη[:,end], simulation=simulationName, model="EpiFilter"))
    
    # (Marginal results)    
    pR = EpiFilterMarginalPosterior(pη, pRgivenη)
    (mR, _, lR, uR) = calculateResults(pR, Rgrid)
    
    pC = EpiFilterMarginalPredictive(pη, pRupgivenη, w, Ct, Rgrid, Cgrid; windin=windin)
    (mC, _, lC, uC) = calculateResults(pC, Cgrid)
    
    scoreEFMarginal = CRPS(cumsum(pC, dims=1), Cgrid, Ct; stepwise=true)
    
    out = vcat(out, DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mR, lower=lR, upper=uR, mean_param=mη, lower_param=lη, upper_param=uη, mean_cases=mC, lower_cases=lC, upper_cases=uC, crps=scoreEFMarginal, simulation=simulationName, model="EpiFilter", fit="Marginalised"))
    
    # (Conditional results)
    pRcond = EpiFilterConditionalPosterior(0.1, w, Ct, Rgrid)
    (mRcond, _, lRcond, uRcond) = calculateResults(pRcond, Rgrid)
    
    pCcond = EpiFilterConditionalPredictive(0.1, w, Ct, Rgrid, Cgrid)
    (mCcond, _, lCcond, uCcond) = calculateResults(pCcond, Cgrid)
    
    scoreEFConditional = CRPS(cumsum(pCcond, dims=1), Cgrid, Ct; stepwise=true)
    
    out = vcat(out, DataFrame(t=1:length(Ct), TrueRt=Rt, Cases=Ct, mean=mRcond, lower=lRcond, upper=uRcond, mean_param=0.1, lower_param=NaN, upper_param=NaN, mean_cases=mCcond, lower_cases=lCcond, upper_cases=uCcond, crps=scoreEFConditional, simulation=simulationName, model="EpiFilter", fit="Conditional (eta = 0.1)"))
    
    return(out, ee_parampost, ef_parampost)
    
end