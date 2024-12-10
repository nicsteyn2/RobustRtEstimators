
using DataFrames, CSV, Plots, Distributions

include("../../src/EpiFilter.jl")
include("../../src/scoringRules.jl")
include("../../src/support.jl")

function fitEpiFilterSmoothing(simulation)

    println("Starting simulation $simulation...")
    
    # Load example data
    df = sort(CSV.read("othermethods/exampledata.csv", DataFrame), :t)
    Ct = df.Cases[df.simulation .== simulation]
    TrueRt = df.TrueRt[df.simulation .== simulation]
    
    # Set grids and default serial interval
    w = defaultGenTimeDist()
    Rgrid = LinRange(0.01, 10, 1000)
    ηgrid = LinRange(0.001, 1, 1000)
    Cgrid = 0:(20*maximum(Ct))
    pη0 = ones(length(ηgrid))/length(ηgrid)
    
    # Fit EpiFilter on the grid of η values
    (pη, pRgivenη, pRupgivenη, qRgivenη) = EpiFilterRunAllη(w, Ct, Rgrid, pη0, ηgrid; windin=10, showProgress=true, smoothingposterior=true)

    # Print mean and quantilers of pη[:,end]
    (mη, medη, lη, uη) = calculateResults(pη[:,end], ηgrid)
    modeη = ηgrid[argmax(pη[:,end])]
    println("Posterior η = $modeη (95% CI: $lη, $uη)")
    
    # Fetch filtering posterior distribution for Rt
    pRfilt = EpiFilterMarginalPosterior(pη, pRgivenη)
    (mF, medF, lF, uF) = calculateResults(pRfilt, Rgrid)
    
    # Fetch filtering predictive distribution for Ct
    pCfilt = EpiFilterMarginalPredictive(pη, pRupgivenη, w, Ct, Rgrid, Cgrid; retrospective=false)
    (mCF, medCF, lCF, uCF) = calculateResults(pCfilt, Cgrid)
    
    # Fetch smoothing posterior distribution for Rt
    pRsmooth = EpiFilterMarginalPosterior(pη[:,end], qRgivenη)
    (mS, medS, lS, uS) = calculateResults(pRsmooth, Rgrid)
    
    # Fetch smoothing predictive distribution for Ct
    pCsmooth = EpiFilterMarginalPredictive(pη, qRgivenη, w, Ct, Rgrid, Cgrid; retrospective=true)
    (mCS, medCS, lCS, uCS) = calculateResults(pCsmooth, Cgrid)
    
    # Plot
    p = plot(2:100, mF[2:end], ribbon=(mF[2:end]-lF[2:end], uF[2:end]-mF[2:end]), label="Filtering dist.")
    p = plot!(p, 1:100, mS, ribbon=(mS-lS, uS-mS), label="Smoothing dist.")
    p = plot!(p, 1:100, TrueRt, color=:black, label="True Rt")
    
    p2 = plot(3:100, mCF[3:end], ribbon=(mCF[3:end]-lCF[3:end], uCF[3:end]-mCF[3:end]), label="Filtering dist.")
    p2 = plot!(p2, 2:100, mCS[2:end], ribbon=(mCS[2:end]-lCS[2:end], uCS[2:end]-mCS[2:end]), label="Smoothing dist.")
    p2 = plot!(p2, 1:100, Ct, color=:black, label="True incidence")
    
    display(plot(p, p2, layout=(2,1)))
    
    # Calculate CRPS
    CRPSfilt = CRPS(cumsum(pCfilt, dims=1), Cgrid, Ct; windin=10, stepwise=true)
    CRPSsmooth = CRPS(cumsum(pCsmooth, dims=1), Cgrid, Ct; windin=10, stepwise=true)
    
    # Save results
    df_out_filt_Rt = DataFrame(t=1:100, Mean=mF, Lower=lF, Upper=uF, variable="R", CRPS=NaN, simulation=simulation, TrueValue=TrueRt, posterior="Filtering")
    df_out_smooth_Rt = DataFrame(t=1:100, Mean=mS, Lower=lS, Upper=uS, variable="R", CRPS=NaN, simulation=simulation, TrueValue=TrueRt, posterior="Smoothing")
    df_out_filt_Ct = DataFrame(t=1:100, Mean=mCF, Lower=lCF, Upper=uCF, variable="I", CRPS=CRPSfilt, simulation=simulation, TrueValue=Ct, posterior="Filtering")
    df_out_smooth_Ct = DataFrame(t=1:100, Mean=mCS, Lower=lCS, Upper=uCS, variable="I", CRPS=CRPSsmooth, simulation=simulation, TrueValue=Ct, posterior="Smoothing")
    df_out = vcat(df_out_filt_Rt, df_out_smooth_Rt, df_out_filt_Ct, df_out_smooth_Ct)

    df_params_out = DataFrame(simulation=simulation, parameter="eta", mean=mη, mode=modeη, lower=lη, upper=uη)
    
    return(df_out, df_params_out)
    
end

(df_randomwalk, df_randomwalk_params) = fitEpiFilterSmoothing("Random walk")
(df_sinusoidal, df_sinusoidal_params) = fitEpiFilterSmoothing("Sinusoidal")
(df_stepchange, df_stepchange_params) = fitEpiFilterSmoothing("Step change")

df_all = vcat(df_randomwalk, df_sinusoidal, df_stepchange)
df_all_params = vcat(df_randomwalk_params, df_sinusoidal_params, df_stepchange_params)

CSV.write("othermethods/EpiFilterSmoothing/results.csv", df_all)
CSV.write("othermethods/EpiFilterSmoothing/params.csv", df_all_params)