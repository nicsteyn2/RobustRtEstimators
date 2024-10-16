
using ProgressMeter


# --------------------------------------------------
# ------------- Core EpiEstim functions ------------
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


function EpiEstim(k::Int, w::Vector, C::Vector; a0=1, b0=1/5)
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w::Vector, C::Vector)
    
    # Fetch posterior parameters and posterior dist
    (a, b) = EpiEstimPosteriorParams(k, w, C; a0=a0, b0=b0)
    pRgivenk = Gamma.(a, 1 ./ b)
    
    return(pRgivenk)
    
end


function EpiEstimPosteriorParams(k::Int, w::Vector, C::Vector; a0=1, b0=1/5)
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w, C)
    
    # Preallocate output
    a = zeros(length(C))
    b = zeros(length(C))
    
    # Run for each step (careful about parameters for tt<k)
    if k == 0
        a[1:length(C)] .= a0
        b[1:length(C)] .= b0
    else
        a[1] = a0
        b[1] = b0 # At t = 1 there is no past data to inform b
        for tt = 2:length(C)
            a[tt] = a0 + sum(C[max(tt-k+1,1):tt])
            b[tt] = b0 + sum(Λ[max(tt-k+1,1):tt])
        end
    end
    
    return(a, b)
    
end


function EpiEstimPosteriorParams(k::Vector, w::Vector, C::Vector; a0=1, b0=1/5)
    
    # Pre-allocate output
    a = zeros(length(k), length(C))
    b = zeros(length(k), length(C))
    
    # Iterate over k values
    for (ii, kk) in enumerate(k)
        (a[ii,:], b[ii,:]) = EpiEstimPosteriorParams(kk, w, C; a0=a0, b0=b0)
    end
    
    return(a,b)
    
end



# --------------------------------------------------
# ------------- Likelihood calculators -------------
# --------------------------------------------------

function EpiEstimLogLik(k::Int, w, C; a0=1, b0=1/5, windin=10)
    
    # Pre-allocate output
    stepwiseloglik = zeros(length(C))
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w, C)
    
    # Fetch posterior params for Rt
    (a, b) = EpiEstimPosteriorParams(k-1, w, C; a0=a0, b0=b0) # Evaluate at k-1 
    
    # And calculate stepwise log-likelihood
    r = a[1:end-1]
    p = b[1:end-1] ./ (Λ[2:end] .+ b[1:end-1])
    stepwiseloglik[2:end] = log.(pdf.(NegativeBinomial.(r, p), C[2:end]))
    if windin > 0
        stepwiseloglik[1:windin] .= 0
    end
    
    # Return cumulative sum and stepwise
    return(cumsum(stepwiseloglik), stepwiseloglik)
    
end


function EpiEstimLogLik(k::AbstractVector, w, C; a0=1, b0=1/5, windin=10)
    
    # Pre-allocate output
    loglik = zeros(length(k), length(C))
    stepwiseloglik = zeros(length(k), length(C))
    
    # Iterate over all values of k
    for (ii, kk) in enumerate(k)
        (loglik[ii,:], stepwiseloglik[ii,:]) = EpiEstimLogLik(kk, w, C; a0=a0, b0=b0, windin=windin)
    end
    
    return(loglik, stepwiseloglik)
    
end


# --------------------------------------------------
# -------- EpiEstim Conditional Posteriors ---------
# --------------------------------------------------

function EpiEstimConditionalPosterior(k::Int, w, C; a0=1, b0=1/5)
    
    return(EpiEstim(k, w, C; a0=a0, b0=b0))
    
end


function EpiEstimConditionalPredictive(k::Int, w, C; a0=1, b0=1/5)
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w, C)
    
    # Fetch posterior parameters
    (a, b) = EpiEstimPosteriorParams(k .- 1, w, C; a0=a0, b0=b0)
    
    # Get NegBin params
    r = a[1:end-1]
    p = b[1:end-1] ./ (Λ[2:end] .+ b[1:end-1])
    pC = vcat(NegativeBinomial.(r[1], p[1]), NegativeBinomial.(r, p))
    
    return(pC)
    
end


# --------------------------------------------------
# ---------- EpiEstim Marginal Posteriors ----------
# --------------------------------------------------

function EpiEstimMarginalPosterior(w, C, Rgrid; a0=1, b0=1/5, kMax=30, windin=10, kPrior=missing)
    
    # Setup and pre-allocate
    kvals = collect(1:kMax)
    pR = zeros(length(Rgrid), length(C))
    
    # Calculate log-lik for all k, t, convert to likelihood, and normalise
    (loglik, _) = EpiEstimLogLik(kvals, w, C; a0=a0, b0=b0, windin=windin)
    loglik = loglik .- maximum(loglik, dims=1)
    pK = missing
    if ismissing(kPrior)
        pK = exp.(loglik)./sum(exp.(loglik), dims=1) # This assumes a Uniform prior distribution on 1:kMax
    else
        pK = exp.(loglik .- maximum(loglik, dims=1)) .* kPrior
        pK = pK ./ sum(pK, dims=1)
    end
    
    # Find the marginal posterior density on our grid
    (a, b) = EpiEstimPosteriorParams(kvals, w, C; a0=a0, b0=b0)
    for (ii, Rt) in enumerate(Rgrid)
        pR[ii,:] = sum(pdf.(Gamma.(a, 1 ./ b), Rt) .* pK, dims=1)
    end
    
    # Normalise columns and set wind-in period to prior distribution
    pR = pR ./ sum(pR, dims=1)
    pR[:,1:windin] .= 1/length(pR[:,1])
    
    return(pR)
    
end


function EpiEstimMarginalPredictive(w, C, Cgrid; a0=1, b0=1/5, kMax=30, retrospective=false, windin=10, kPrior=missing)
    
    # Setup and pre-allocate
    kvals = collect(1:kMax)
    pC = zeros(length(Cgrid), length(C))
    
    # Precompute force of infection
    Λ = CalculateInfectionPressure(w, C)
    
    # Calculate log-lik for all k, t, convert to likelihood, and normalise
    (loglik, _) = EpiEstimLogLik(kvals, w, C; a0=a0, b0=b0, windin=windin)
    loglik = loglik .- maximum(loglik, dims=1)
    pK = missing
    if ismissing(kPrior)
        pK = exp.(loglik)./sum(exp.(loglik), dims=1) # This assumes a Uniform prior distribution on 1:kMax
    else
        pK = exp.(loglik .- maximum(loglik, dims=1)) .* kPrior
        pK = pK ./ sum(pK, dims=1)
    end
    
    # If we are analysing it retrospectively, we use the final estimate of K
    if retrospective
        for ii = 1:(size(pK)[2]-1)
            pK[:,ii] = pK[:,end]
        end
    end
    
    # Find the marginal posterior density on our grid
    (a, b) = EpiEstimPosteriorParams(kvals .- 1, w, C; a0=a0, b0=b0)
    for tt = 2:length(C)
        r = a[:,tt-1]
        p = b[:,tt-1] ./ (Λ[tt] .+ b[:,tt-1])
        pC[:,tt] = [sum(pdf.(NegativeBinomial.(r, p), Ct) .* pK[:,tt]) for Ct in Cgrid]
    end
    
    # Normalise the grid and check that the probability on the upper bound is negligible
    pC[:,1:windin] .= 1
    pC = pC ./ sum(pC, dims=1)
    if any(pC[end, (windin+1):end] .> 1e-3)
        error("Non-negligible probability on upper-bound of Cgrid. Increase max value.")
    end
    if any(pC[end, (windin+1):end] .> 1e-5)
        @warn("Approaching non-negligible probability on upper bound of Cgrid. This may cause an error in future.")
    end
    
    return(pC)
    
end


