
# EpiLPS (MAP)

Simple worked example fitting the MAP version of [EpiLPS](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010618).

In addition to producing posterior means and credible intervals of $R_t$, we also provide methods to find the smoothing predictive distribution for reported cases $C_t$, as well as an estimator of the CRPS for this smoothing predictive distribution.

Full methods are included in the supplementary material of [our paper](https://doi.org/10.1093/aje/kwaf165).

This folder contains the following files:

|File|Description|
|---|---|
|``EpiLPSMAP-example.R``| **Start here.** R script to fit EpiLPS (MAP) to the example data.|
|``LPSMAP-incidence.R``| R script providing the function ``observedIncidenceCurve(LPSfit)`` which computes the posterior mean and credible intervals of the posterior predictive distribution for reported cases.|
|``LPSMAP-CRPS.R``| R script providing the function ``computeCRPS(LPSfit)`` which computes the CRPS of the smoothing predictive distribution of reported cases.|
|``KercubicBspline.cpp``| C++ code for calculating the B-spline basis functions. **Copyright restrictions do not permit us to include the code in this repository. Please download it from the [original source](https://github.com/oswaldogressani/EpiLPS/tree/main/src) to this folder.**|
|``results.csv``| An output. Results of time-varying states.|
|``params.csv``| An output. Posterior values of selected model parameters.|