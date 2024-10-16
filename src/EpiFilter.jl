
using ProgressMeter

include("fastDists.jl")

# --------------------------------------------------
# ------------ Core EpiFilter functions ------------
# --------------------------------------------------


function CalculateInfectionPressure(w::Vector, C::Vector; normalised=true)
    
    Λ = missing
    if normalised
        Λ = [sum(C[(tt-1):-1:1] .* w[1:(tt-1)])/sum(w[1:max(tt-1, 1)]) for tt = 1:length(C)]
    else
        Λ = [sum(C[(tt-1):-1:1] .* w[1:(tt-1)]) for tt = 1:length(C)]
    end
    return(Λ)
    
end


function EpiFilterForwards(η, w, C, Rgrid)
    
    # Extract frequently used parameters
    nPts = length(Rgrid)
    nDay = length(C)
    
    # Pre-allocate memory
    pR = zeros(nPts, nDay) # Filtering distribution for Rt
    pObs = zeros(nPts, nDay) # Observation distribution
    pRup = zeros(nPts, nDay) # Update distribution for Rt
    
    # Set initial state distribution
    pR[:,1] = ones(nPts)/nPts # Start the filter at a uniform distribution
    
    # Precompute state transition distribution
    pstate = zeros(nPts, nPts)
    for ii = 1:nPts
        pstate[:,ii] = pdf.(Normal(Rgrid[ii], η .* sqrt(Rgrid[ii])), Rgrid)
        pstate[:,ii] = pstate[:,ii]/sum(pstate[:,ii])
    end
    
    # Precompute infection pressure terms
    Λ = CalculateInfectionPressure(w, C)
    
    # Iterate over time
    for tt = 2:nDay
        
        # Predict, calculate observation weights, and update
        pRup[:,tt] = pstate * pR[:,tt-1]
        pObs[:,tt] = fastPoissonPDF(C[tt], Rgrid .* Λ[tt])
        pR[:,tt] = pRup[:,tt] .* pObs[:,tt]
        
        # Normalise
        s = sum(pR[:,tt])
        if s > 0
            pR[:,tt] = pR[:,tt]/s
        else
            @warn("EpiFilter failed. This usually occurs when η is too small, resulting in zero observation probability on the support of the update distribution.")
            return(zeros(nPts, nDay), zeros(nPts, nDay), zeros(nPts, nDay))
        end
        pRup[:,tt] = pRup[:,tt] ./ sum(pRup[:,tt])
        
    end
    
    return(pR, pObs, pRup)
    
end


function EpiFilterBackwards(η, w, pR, pRup, Rgrid)
    
    nPts = length(Rgrid)
    nDay = size(pR)[2]
    
    # Allocate memory
    qR = zeros(nPts, nDay)
    qR[:,end] = pR[:,end]
    
    # Copy pRup (so we don't edit the input) and remove zeeros
    pRup = copy(pRup)
    pRup[pRup.==0] .= 1e-8
    
    # Precompute state transition distribution
    pstate = zeros(nPts, nPts)
    for ii = 1:nPts
        pstate[:,ii] = pdf.(Normal(Rgrid[ii], η .* sqrt(Rgrid[ii])), Rgrid)
        pstate[:,ii] = pstate[:,ii]/sum(pstate[:,ii])
    end
    
    # Iterate backwards over time
    for tt = nDay-1:-1:1
        
        # Integral term in smoother
        integ = qR[:,tt+1]./pRup[:,tt+1]
        integ = pstate*integ
        
        # Smoothed posterior over Rgrid
        qR[:,tt] = pR[:,tt] .* integ
        qR[:,tt] = qR[:,tt]/sum(qR[:,tt])
        
    end
    
    return(qR)
    
end



# --------------------------------------------------
# ------------- Likelihood calculators -------------
# --------------------------------------------------

function EpiFilterLogLik(η::Float64, w, C, Rgrid; windin=1)
    
    (_, pObs, pRup) = EpiFilterForwards(η, w, C, Rgrid)
    stepwiseloglik = log.(vec(sum(pRup .* pObs, dims=1)))
    if windin > 0
        stepwiseloglik[1:windin] .= 0
    end
    return(cumsum(stepwiseloglik), stepwiseloglik)
    
end


function EpiFilterLogLik(η::AbstractVector, w, C, Rgrid; windin=1)
    
    # Pre-allocate output
    loglik = zeros(length(η), length(C))
    stepwiseloglik = zeros(length(η), length(C))
    
    # Iterate over all values of η
    for (ii, ηi) in enumerate(η)
        (loglik[ii,:], stepwiseloglik[ii,:]) = EpiFilterLogLik(ηi, w, C, Rgrid; windin=windin)
    end
    
    return(loglik, stepwiseloglik)
    
end



# --------------------------------------------------
# -------- EpiFilter Conditional Posteriors ---------
# --------------------------------------------------

function EpiFilterConditionalPosterior(η, w, C, Rgrid; smoothing=false)
    
    # Get filtering posterior on Rt
    (pR, _, pRup) = EpiFilterForwards(η, w, C, Rgrid)
    
    # Run backwards smoother if enabled (only to be used for retrospective analysis)
    if smoothing
        pR = EpiFilterBackwards(η, w, pR, pRup, Rgrid)
    end
    
    return(pR)
    
end


function EpiFilterConditionalPredictive(η, w, C, Rgrid, Cgrid)
    
    # Pre-allocate
    pC = zeros(length(Cgrid), length(C))
    
    # Precompute infection pressure terms
    Λ = CalculateInfectionPressure(w, C)
    
    # Run EpiFilter to get P(R_t|C_{1:t-1}, η)
    (_, _, pRup) = EpiFilterForwards(η, w, C, Rgrid)
    
    # Calculate the predictive distribution over Cgrid
    pC[:,1] .= 1/length(pC[:,1])
    #TODO: Add a version for non-UnitRange Cgrid
    for tt = 2:length(C)
        λt = Rgrid .* Λ[tt]
        # pC[:,tt] = [sum(fastPoissonPDF(Ct, λt) .* pRup[:,tt]) for Ct in Cgrid]
        pC[:,tt] = sum(veryFastPoissonPDF(Cgrid, λt).* pRup[:,tt], dims=1)
    end
    
    # Normalise the grid and check that the probability on the upper bound is negligible
    pC = pC ./ sum(pC, dims=1)
    if any(pC[end, 2:end] .> 1e-3)
        error("Non-negligible probability on upper-bound of Cgrid. Increase max value.")
    end
    if any(pC[end, 2:end] .> 1e-5)
        @warn("Approaching non-negligible probability on upper bound of Cgrid. This may cause an error in future.")
    end
    
    return(pC)
    
end




# --------------------------------------------------
# -------- EpiFilter Marginal Posteriors ---------
# --------------------------------------------------

# As we need to run EpiFilter on many η values, it is helpful to compute everything once
function EpiFilterRunAllη(w, C, Rgrid, pη0, ηgrid; showProgress=true, windin=1)
    
    # Pre-allocate output. The "cond" suffixes indicate results are conditional on η
    pRgivenη = zeros(length(Rgrid), length(C), length(ηgrid)) 
    pRupgivenη = zeros(length(Rgrid), length(C), length(ηgrid))
    stepwiseloglik = zeros(length(ηgrid), length(C))
    
    # Iterate over all values of η
    ProgBar = Progress(length(ηgrid), dt=1, barlen=50, desc="Running EpiFilter on all η values...", enabled=showProgress)
    for (ii, ηi) in enumerate(ηgrid)
        (pRgivenη[:,:,ii], pObs, pRupgivenη[:,:,ii]) = EpiFilterForwards(ηi, w, C, Rgrid)
        stepwiseloglik[ii,:] = vec(log.(sum(pRupgivenη[:,:,ii] .* pObs, dims=1)))
        next!(ProgBar)
    end

    # Allow for windin in calculating log-likelihood
    if windin > 0
        stepwiseloglik[:,1:windin] .= 0
    end
    
    # Convert stepwise loglik into cumulative loglik
    loglik = cumsum(stepwiseloglik, dims=2)
    loglik = loglik .- maximum(loglik, dims=1) # Normalise loglik to avoid numerical issues
    
    # Convert to posterior and normalise
    pη = exp.(loglik) .* pη0 # Likelihood * prior
    pη = pη ./ sum(pη, dims=1)
    
    return(pη, pRgivenη, pRupgivenη)
    
end


function EpiFilterMarginalPosterior(pη, pRgivenη)
    
    # Pre-allocate
    pR = zeros(size(pRgivenη[:,:,1]))
    
    # Iterate over time
    for tt = 1:size(pR)[2]
        pR[:,tt] = sum(pRgivenη[:,tt,:] .* pη[:,tt]', dims=2)
    end
    
    return(pR)
    
end


function EpiFilterMarginalPredictive(pη, pRupgivenη, w, C, Rgrid, Cgrid; retrospective=false, windin=1)
    
    # Pre-allocate
    pRup = zeros(size(pRupgivenη[:,:,1]))
    pC = zeros(length(Cgrid), length(C))
        
    # If retrospective we use only the final entry in pη
    if retrospective
        for ii = 1:(size(pη)[2]-1)
            pη[:,ii] = pη[:,end]
        end
    end
    
    # Precompute force of infection
    Λ = CalculateInfectionPressure(w, C)
    
    # Calculate marginal pRup dist
    for tt = 1:length(C)
        pRup[:,tt] = sum(pRupgivenη[:,tt,:] .* pη[:,tt]', dims=2)
    end
    pRup[:,1] .= 1/length(pRup[:,1])
    
    # Calculate the predictive distribution over Cgrid
    pC[:,1] .= 1/length(pC[:,1])
    for tt = 2:length(C)
        λt = Rgrid .* Λ[tt]
        # pC[:,tt] = [sum(fastPoissonPDF(Ct, λt) .* pRup[:,tt]) for Ct in Cgrid]
        pC[:,tt] = sum(veryFastPoissonPDF(Cgrid, λt).* pRup[:,tt], dims=1)
    end

    
    # Normalise the grid and check that the probability on the upper bound is negligible
    pC = pC ./ sum(pC, dims=1)
    if any(pC[end, (windin+1):end] .> 1e-3)
        error("Non-negligible probability on upper-bound of Cgrid. Increase max value.")
    end
    if any(pC[end, (windin+1):end] .> 1e-5)
        @warn("Approaching non-negligible probability on upper bound of Cgrid. This may cause an error in future.")
    end
    
    return(pC)
    
end

