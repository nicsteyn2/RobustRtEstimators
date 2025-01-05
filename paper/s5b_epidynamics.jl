
using Plots, Random


include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/Simulations.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")


function runVaryingRW( ; nIter=10)
    
    ηvals = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    
    # Save outputs
    out_all = DataFrame()
    ee_parampost_all = DataFrame()
    ef_parampost_all= DataFrame()
    
    # Iterate over eta values
    ProgBar = Progress(length(ηvals)*nIter, dt=1, barlen=50, desc="Running multiple RW simulations")
    for jj in 1:nIter
        
        for (ii, η) in enumerate(ηvals)

            println()
            println()
            println()
            println("Starting new simulation with η = $η...")
            
            (Ct, Rt) = simulateRandomWalkRt(η=η; minCt=5, maxCt=5000)

            println("Starting model fitting...")
            (out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct, Rt, "Random walk"; windin=10, showProgress=false)
            
            out_tmp[!,:TrueParam] .= η
            ee_parampost_tmp[!,:TrueParam] .= η
            ef_parampost_tmp[!,:TrueParam] .= η
            
            out_tmp[!,:iter] .= jj
            ee_parampost_tmp[!,:iter] .= jj
            ef_parampost_tmp[!,:iter] .= jj
            
            out_all = vcat(out_all, out_tmp)
            ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
            ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)
            
            CSV.write("paper/outputs/s5b_varying_rw.csv", out_all)
            CSV.write("paper/outputs/s5b_varying_rw_eepost.csv", ee_parampost_all)
            CSV.write("paper/outputs/s5b_varying_rw_efpost.csv", ef_parampost_all)
            
            next!(ProgBar)
            
        end
    end
    
end



function runVaryingSinusoidal( ; nIter=10)
    
    ωvals = [10, 15, 20, 25, 30, 35, 40, 45, 50]
    
    # Save outputs
    out_all = DataFrame()
    ee_parampost_all = DataFrame()
    ef_parampost_all= DataFrame()
    
    # Iterate over eta values
    ProgBar = Progress(length(ωvals)*nIter, dt=1, barlen=50, desc="Running multiple sinusoidal simulations")
    for jj in 1:nIter
        
        for (ii, ω) in enumerate(ωvals)

            println()
            println()
            println()
            println("Starting new simulation with ω = $ω...")
            
            (Ct, Rt) = simulateSinusoidalRt(ω=ω)

            println("Starting model fitting...")
            (out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct, Rt, "Sinusoidal"; windin=10, showProgress=false)
            
            out_tmp[!,:TrueParam] .= ω
            ee_parampost_tmp[!,:TrueParam] .= ω
            ef_parampost_tmp[!,:TrueParam] .= ω
            
            out_tmp[!,:iter] .= jj
            ee_parampost_tmp[!,:iter] .= jj
            ef_parampost_tmp[!,:iter] .= jj
            
            out_all = vcat(out_all, out_tmp)
            ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
            ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)
            
            CSV.write("paper/outputs/s5b_varying_sinusoidal.csv", out_all)
            CSV.write("paper/outputs/s5b_varying_sinusoidal_eepost.csv", ee_parampost_all)
            CSV.write("paper/outputs/s5b_varying_sinusoidal_efpost.csv", ef_parampost_all)
            
            next!(ProgBar)
            
        end
    end
    
end




function runVaryingStepChange( ; nIter=10)
    
    numberOfChanges = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    changePoints = Int.(round.(100 ./ (numberOfChanges .+ 1) .+ 1))
    
    # Save outputs
    out_all = DataFrame()
    ee_parampost_all = DataFrame()
    ef_parampost_all= DataFrame()
    
    # Iterate over eta values
    ProgBar = Progress(length(changePoints)*nIter, dt=1, barlen=50, desc="Running multiple step change sims...")
    for jj in 1:nIter
        
        for (ii, changePoint) in enumerate(changePoints)

            println()
            println()
            println()
            println("Starting new simulation with changePoints every = $changePoint...")
            
            (Ct, Rt) = simulateStepChangeRt(; tChange = changePoint)

            println("Starting model fitting...")
            (out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct, Rt, "Step change"; windin=10, showProgress=false)
            
            out_tmp[!,:TrueParam] .= changePoint
            ee_parampost_tmp[!,:TrueParam] .= changePoint
            ef_parampost_tmp[!,:TrueParam] .= changePoint
            
            out_tmp[!,:iter] .= jj
            ee_parampost_tmp[!,:iter] .= jj
            ef_parampost_tmp[!,:iter] .= jj
            
            out_all = vcat(out_all, out_tmp)
            ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
            ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)
            
            CSV.write("paper/outputs/s5b_varying_stepchanges.csv", out_all)
            CSV.write("paper/outputs/s5b_varying_stepchanges_eepost.csv", ee_parampost_all)
            CSV.write("paper/outputs/s5b_varying_stepchanges_efpost.csv", ef_parampost_all)
            
            next!(ProgBar)
            
        end
    end
    
end



# runVaryingRW()
# runVaryingSinusoidal()
# runVaryingStepChange()