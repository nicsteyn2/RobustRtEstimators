
using Plots, Random, Measures

include("../src/EpiEstim.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

# Load data
(Ct, w) = loadData("NZCOVID")

# Specify grids
Rgrid = LinRange(0.01, 10, 1000)
Cgrid = 0:1:(10*maximum(Ct))

# Specify windin period
windin = 10


function EpiEstimAPE(k::Int, w, C; a0=1, b0=1/5, windin=10)
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w, C)
    
    # Find parameters of P_APE(Rt|C_{1:t})
    (a, b) = EpiEstimPosteriorParams(k, w, C; a0=a0, b0=b0)

    # Shift by one day and calculate NegBin params
    r = a[1:end-1]
    p = b[1:end-1] ./ (Λ[2:end] .+ b[1:end-1])

    # And then calculate the APE metric, allowing for a windin period of 10 days
    StepwiseAPE = -log.(pdf.(NegativeBinomial.(r, p), C[2:end]))
    return(sum(StepwiseAPE[windin+1:end]))
    
end


kvals = 1:30
APE = [EpiEstimAPE(k, w, Ct) for k in kvals]
(loglik, _) = EpiEstimLogLik(kvals, w, Ct)

p = plot(xlabel="k", ylabel="Value of metric", dpi=300, size=(600,300), margin=3mm, minorgrid=true)
plot!(p, kvals, APE, label="APE metric", marker=:circle)
plot!(p, kvals, -1 * loglik[:,end], label="Negative log-likelihood", marker=:circle, legend=:bottomright)
png(p, "paper/figures/s8_APEestim.png")