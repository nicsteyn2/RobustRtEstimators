
using Random, Distributions

# ---------------------------------------------------------------------
# ---------------- Epidemic simulator functions -----------------------
# ---------------------------------------------------------------------

function simulateRandomWalkRt( ; η = 0.1, w=missing, tMax = 100, R0 = 1, C0=100, windin=20, Rt_min=0.1, Rt_max=5, minCt=1, maxCt=5000, forceNoElimination=true)
    
    isFinished = false
    maxAttempts = 1e6
    
    # Load generation time distribution if it is not otherwise specified
    if ismissing(w)
        w = defaultGenTimeDist(tMax)
    end

    # Pre-allocate output
    Rt = zeros(tMax + windin)
    Rt[1:windin] .= R0
    Ct = Int.(zeros(tMax + windin))
    Ct[1:windin] .= C0
    
    while !isFinished && maxAttempts >= 0
        
        # Pre-allocate output
        Rt = zeros(tMax + windin)
        Rt[1:windin] .= R0
        Ct = Int.(zeros(tMax + windin))
        Ct[1:windin] .= C0
        
        # Simulate the epidemic
        for tt = (windin+1):(tMax+windin)
            Rt[tt] = rand(Truncated(Normal(Rt[tt-1], η * sqrt(Rt[tt-1])), Rt_min, Rt_max))
            Λ = sum(Ct[(tt-1):-1:(tt-windin)] .* w[1:windin])
            Ct[tt] = rand(Poisson(Rt[tt] * Λ))
        end

        # Check finishing conditions
        if (minimum(Ct) > minCt) && (maximum(Ct) < maxCt)
            isFinished = true
        end

        maxAttempts = maxAttempts - 1
        
    end
    
    return(Ct[(windin+1):end], Rt[(windin+1):end])
    
end



function simulateSinusoidalRt( ; A=0.5, ω=21, w=missing, tMax=100, R0=1, C0=100, windin=50)
    
    # Load generation time distribution if it is not otherwise specified
    if ismissing(w)
        w = defaultGenTimeDist(tMax)
    end
    
    # Pre-allocate output
    Rt = zeros(tMax + windin)
    Ct = Int.(zeros(tMax+windin))
    Ct[1:windin] .= C0
    
    # Simulate the epidemic
    for tt = (windin+1):(tMax+windin)
        Rt[tt] = R0 + A * sin(((2*pi)/ω) * (tt-(windin+1)))
        Λ = sum(Ct[(tt-1):-1:(tt-windin)] .* w[1:windin])
        Ct[tt] = rand(Poisson(Rt[tt] * Λ))
    end
    
    (Ctout, Rtout) = Ct[(windin+1):end], Rt[(windin+1):end]
    Ctout[1] = max(Ctout[1], 1) # Ensure there is at least one case at initialisation
    
    return(Ctout, Rtout)
    
end


function simulateStepChangeRt( ; Ra = 2/3, Rb=3/2, tChange=28, w=missing, tMax=100, C0=100, windin=50)
    
    # Load generation time distribution if it is not otherwise specified
    if ismissing(w)
        w = defaultGenTimeDist(tMax)
    end
    
    # Pre-allocate output
    Rt = zeros(tMax + windin)
    Ct = Int.(zeros(tMax+windin))
    Rt[windin] = Ra
    Ct[1:windin] .= C0
    
    # Simulate the epidemic
    for tt = (windin+1):(tMax+windin)
        
        # Start at Rt = Ra, then swap between Ra and Rb every tChange days
        if mod(tt - windin, tChange) == 0
            if Rt[tt-1] == Ra
                Rt[tt] = Rb
            else
                Rt[tt] = Ra
            end
        else
            Rt[tt] = Rt[tt-1]
        end
        
        Λ = sum(Ct[(tt-1):-1:(tt-windin)] .* w[1:windin])
        Ct[tt] = rand(Poisson(Rt[tt] * Λ))
    end
    
    return(Ct[(windin+1):end], Rt[(windin+1):end])
    
end



# ---------------------------------------------------------------------
# ----------------- Simulate observation noise ------------------------
# ---------------------------------------------------------------------

function simulateBinomialObservationNoise(Ct; p=0.5)

    Cobs = rand.(Binomial.(Ct, p))
    return(Cobs)

end

function simulateBetaBinomialObservationNoise(Ct; p=0.5, N=100)

    a = Int(round(p * N))
    b = N - a
    Cobs = rand.(BetaBinomial.(Ct, a, b))
    return(Cobs)

end