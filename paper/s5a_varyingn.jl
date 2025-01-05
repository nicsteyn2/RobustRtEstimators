
using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

function runVaryingC0_sinusoidal(; nIter=10)
    
    C0vals = [25, 50, 100, 200, 400, 800, 1600, 3200]

    # Save outputs
    out_all = DataFrame()
    ee_parampost_all = DataFrame()
    ef_parampost_all= DataFrame()
    
    # Iterate over C0 values
    ProgBar = Progress(length(C0vals)*nIter, dt=1, barlen=50, desc="Running multiple sinusoidal simulations")
    for jj in 1:nIter
        
        for (ii, C0) in enumerate(C0vals)
            
            println()
            println()
            println()
            println("Starting new simulation with C0 = $C0...")
            
            (Ct, Rt) = simulateSinusoidalRt(ω=21, C0=C0)
            
            println("Starting model fitting...")
            (out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct, Rt, "Sinusoidal"; windin=10, showProgress=false)
            
            out_tmp[!,:TrueParam] .= C0
            ee_parampost_tmp[!,:TrueParam] .= C0
            ef_parampost_tmp[!,:TrueParam] .= C0
            
            out_tmp[!,:iter] .= jj
            ee_parampost_tmp[!,:iter] .= jj
            ef_parampost_tmp[!,:iter] .= jj
            
            out_all = vcat(out_all, out_tmp)
            ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
            ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)
            
            CSV.write("paper/outputs/s5a_varying_C0_sinusoidal.csv", out_all)
            CSV.write("paper/outputs/s5a_varying_C0_sinusoidal_eepost.csv", ee_parampost_all)
            CSV.write("paper/outputs/s5a_varying_C0_sinusoidal_efpost.csv", ef_parampost_all)
            
            next!(ProgBar)
            
        end
    end
    
end


runVaryingC0_sinusoidal(nIter=2)

