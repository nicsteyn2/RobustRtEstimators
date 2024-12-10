
# EpiLPS (MALA)

Simple worked example fitting the MALA version of [EpiLPS](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010618).

In addition to producing posterior means and credible intervals of $R_t$, we also provide methods to find the smoothing predictive distribution for reported cases $C_t$, as well as an estimator of the CRPS for this smoothing predictive distribution.

Full methods are included in the supplementary material of [our paper](#TODO) and on the [corresponding website](#TODO).

This folder contains the following files:

|File|Description|
|---|---|
|``EpiLPSMALA-example.R``| **Start here.** R script to fit EpiLPS (MAP) to the example data.|
|``LPSMALA-incidence.R``| R script providing the function ``observedIncidenceCurve(LPSfit)`` which computes the posterior mean and credible intervals of the posterior predictive distribution for reported cases.|
|``LPSMALA-CRPS.R``| R script providing the function ``computeCRPS(LPSfit)`` which computes the CRPS of the smoothing predictive distribution of reported cases.|
|``/packagefiles/``| Additional scripts required to run the example. Please see the details below.|
|``results.csv``| An output. Results of time-varying states.|
|``params.csv``| An output. Posterior values of selected model parameters.|

In order to extract the MCMC samples, we modify the ``estimRmcmc()`` function provided in the original package. The modified function is called ``estimRmcmc_custom()`` and is provided in ``packagefiles/estimRmcmc-custom.R``. Other required scripts from the [github repository for EpiLPS](https://github.com/oswaldogressani/EpiLPS) are also provided/required in the ``packagefiles`` folder:

|File|Description|Included?|
|---|---|---|
|``estimRmcmc-custom.R``| Modified version of the original function ``estimRmcmc()``.|Yes|
|``KercubicBspline.cpp``| C++ code for the cubic B-spline kernel.|No. [Download here](https://github.com/oswaldogressani/EpiLPS/tree/main/src).|
|``KerIncidCheck.R``| Function to check incidence data.|Yes|
|``KerLaplace.cpp``| C++ code for the Laplace kernel.|No. [Download here](https://github.com/oswaldogressani/EpiLPS/tree/main/src).|
|``KerLikelihood.R``| Function to compute the likelihood.|Yes|
|``KerMCMC.R``| Function to run the MCMC/MALA.|Yes|
|``KerMVN.cpp``| C++ code for the multivariate normal kernel.| Yes|
|``KerPtheta.R``| Function to compute the posterior of the B-spline parameter vector.|Yes|
|``KerRpostmcmc.cpp``| C++ code to calculate R from B-spline posterior.|No. [Download here](https://github.com/oswaldogressani/EpiLPS/tree/main/src).|