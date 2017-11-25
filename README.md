# BNPInferenceForMRL
This repository contains the R code for the model presented in the following paper

## ***Nonparametric Bayesian inference for mean residual life functions in survival analysis***

## Authors:
* **Valerie Poynor** California State University, Fullerton <vpoynor@fullerton.edu>  
* **Athanasios Kottas** Univsersity of California, Santa Cruz <thanos@soe.ucsc.edu>
        

        
In this paper (see references for earlier version on arXiv), we consider a Dirichlet Process (DP) Mixture Model:

![DP Mixture Model](/images/DPMM.gif)

with the primary aim of infering on the Mean Residual Life (MRL) function. For this purpose, we assume a Gamma distribution for the kernel, k(), which is shown to be weakly dense in the MRL space.  We mix over both the log(shape) and log(rate) parameters of the Gamma kernel.  We assume a bivariate Normal baseline distribution with conjugate priors for the mean vector (Normal) and covariance matrix (Inverse Wishart). Finally, we place the conjugate Gamma prior on the precision parameter. 

![Priors](/images/Model.gif)

Under the truncated stick-breaking constructive definition of the DP, the MRL function has the following form:

![MRL](/images/MRL.gif)

Note that the Gamma MRL function can be simplified to 

![MRL2](/images/MRL2.gif)

In this repository, we provide the code and data used in the analysis of a dataset describing the survival times
(in days) of patients with small cell lung cancer.  These data are described in (Ying et al., 1995) and originally collected in (Maksymuik et al., 1994). The patients were randomly assigned to
one of two treatments referred to as Arm A and Arm B. Arm A patients received cisplatin (P)
followed by etoposide (E), while Arm B patients received (E) followed by (P). There were a total
of 62 patients in Arm A with 15 right censored survival times, while Arm B consisted of 59
patients with 8 right censored survival times.  

The data can be found in the files: *ArmA.txt* and *ArmB.txt*

The variables recorded in the dataset(s) are:
* age: age of the patient upon entry (in years)
* time: survival time of the patient (in days)
* cens: indicator for survival time - 0 if observed, 1 if right censored

The code is written in the R language and is broken up into two parts:

*ArmAArmBGammaDPMMPart1.R* provides the MCMC algorithm to obtain posterior samples of the model parameters.

*ArmAArmBGammaDPMMPart2.R* provides the algorithm to obtain posterior samples of the survival functionals (density, Survival, hazard rate, and MRL)

The code notation is described in detail of Part 1, but there are two main notational changes from the paper to the code - 

1) L -> N

2) t -> y




# References

Jett, R, J D Earle, J Q Su, F A Diegert, J A Mailliard, C G Kardinal, J E Krook, M H Veeder, and M Wiesenfeld. Sequencing and schedule effects of cisplatin plus etoposide in small-cell lung cancer: results of a North Central Cancer Treatment Group randomized clinical trial. A W Maksymiuk, J Journal of Clinical Oncology 1994 12:1, 70-76

Poynor V. and Kottas A. (2017). Nonparametric Bayesian Inference for Mean Residual Life Functions in Survival Analysis. Submitted to Biostatistics. arXiv:1411.7481 [stat.ME].

Ying, Z., S. H. Jung and L. J. Wei. Survival Analysis with Median Regression Models. Journal of the American Statistical Association, Vol. 90, No. 429 (Mar., 1995), pp.178-184





