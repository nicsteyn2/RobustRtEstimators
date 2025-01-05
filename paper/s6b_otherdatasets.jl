
using Plots, Random

include("../src/EpiEstim.jl")
include("../src/EpiFilter.jl")
include("../src/support.jl")
include("../src/scoringRules.jl")

function fitModelsToRealWorldData()
    
    # Save outputs
    out_all = DataFrame()
    ee_parampost_all = DataFrame()
    ef_parampost_all= DataFrame()
    
    # Iterate over all datasources
    sources = ["flu", "flufilt", "sars", "sarsfilt"]
    
    for source in sources
        
        (Ct, w) = loadData(source)
        Rt = similar(Ct) * NaN
        
        println("Starting model on $source...")
        (out_tmp, ee_parampost_tmp, ef_parampost_tmp) = fitModels(Ct, Rt, source; w=w)
        out_all = vcat(out_all, out_tmp)
        ee_parampost_all = vcat(ee_parampost_all, ee_parampost_tmp)
        ef_parampost_all = vcat(ef_parampost_all, ef_parampost_tmp)
        
        CSV.write("paper/outputs/s6b_otherdatasets.csv", out_all)
        CSV.write("paper/outputs/s6b_otherdatasets_epiestimposteriors.csv", ee_parampost_all)
        CSV.write("paper/outputs/s6b_otherdatasets_epifilterposteriors.csv", ef_parampost_all)
        
    end
    
end

fitModelsToRealWorldData()