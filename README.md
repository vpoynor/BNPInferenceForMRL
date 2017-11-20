# BNPInferenceForMRL
This repository contains the R code for the model presented in the following paper

## ***Nonparametric Bayesian inference for mean residual life functions in survival analysis***

## Authors:
* **Valerie Poynor** California State University, Fullerton <vpoynor@fullerton.edu>  
* **Athanasios Kottas** Univsersity of California, Santa Cruz <thanos@soe.ucsc.edu>
        

        
The model proposed in this paper is a Dirichlet Process Mixture Model using a Gamma kernel distribution.  The primary aim of the paper is to infer on the Mean Residual Life (MRL) function.  

Here, we provide the code and data used in the analysis of a dataset describing the survival times
(in days) of patients with small cell lung cancer.  These data are described in (Ying et al., 1995) and originally collected in (Maksymuik et al., 1994). The patients were randomly assigned to
one of two treatments referred to as Arm A and Arm B. Arm A patients received cisplatin (P)
followed by etoposide (E), while Arm B patients received (E) followed by (P). There were a total
of 62 patients in Arm A with 15 right censored survival times, while Arm B consisted of 59
patients with 8 right censored survival times.  

The data can be found in the files: *ArmA.txt* and *ArmB.txt*

The code is written in the R language andis broken up into two parts:

*ArmAArmBGammaDPMMPart1.R* provides the MCMC algorithm to obtain posterior samples of the model parameters

*ArmAArmBGammaDPMMPart2.R* provides the algorithm to obtain posterior samples of the survival functionals



