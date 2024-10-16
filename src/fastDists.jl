
import SpecialFunctions: loggamma
import LogExpFunctions: xlogy

using BenchmarkTools


# If we are evaluating the Poisson PDF at many values of x and a single value of λ
function fastPoissonPDF(x::AbstractVector, λ::Real)
    logλ = log(λ)
    logpdf = @. x * logλ - λ - loggamma(x + 1)
    return exp.(logpdf)
end

# If we are evaluating the Poisson PDF at a single value of x and many values of λ
function fastPoissonPDF(x::Real, λ::AbstractVector)
    logxkfact = loggamma(x + 1)
    logpdf = @. (x * log.(λ)) - λ - logxkfact
    return exp.(logpdf)
end

# If we are evaluating the Poisson PDF at multiple values of x and λ. Returns results as a matrix of size length(λ) x length(x).
function veryFastPoissonPDF(x::UnitRange, λ::AbstractVector)
    
    if x[1] != 0
        error("The veryFastPoissonPDF() function is designed to work with x = 0:Cmax for some Cmax.")
    end
    
    loglambda = log.(λ)
    X = zeros(length(λ), length(x))
    X[:,1] = -λ
    for kk = 2:length(x)
        X[:,kk] = X[:,kk-1] .+ loglambda .- log(kk-1)
    end
    return(exp.(X))
    
end

