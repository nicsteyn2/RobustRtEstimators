


function CRPS(CDF::Matrix, X::AbstractVector, Obs::AbstractVector; stepwise=false, windin=0)
    
    DRPS = [sum([ (CDF[ii,tt] - (X[ii] >= Obs[tt]))^2 for ii = 1:length(X) ]) for tt = 1:length(Obs)]
    
    if stepwise
        return(DRPS)
    else
        return(mean(DRPS[(windin+1):end]))    
    end

end


function CRPS(D::Vector{<:Distribution}, Obs::AbstractVector; stepwise=false, windin=0)

    maxValue = maximum([quantile(D[ii], 0.999) for ii = 1:length(D)])
    X = collect(0:maxValue)
    DRPS = [sum( [ (cdf(D[tt], X[ii]) - (X[ii] >= Obs[tt]))^2 for ii=1:length(X) ] ) for tt=1:length(D)]

    if stepwise
        return(DRPS)
    else
        return(mean(DRPS[(windin+1):end]))
    end

end