# Spatial soft sweeps reveals distinct continent-wide demographic-dispersal evolutionary strategies of human and mosquito populations across Africa

This repository contains R scripts to perform Approximate Bayesian Computation (ABC) analyses on spatial soft sweeps simulations for *Anopheles* and humans.

## **Scripts**

### **Anopheles Folder**

    - `Rejection_Wassertstein.R`: R script to perform ABC using a simple rejection approach on Wasserstein distances for different sampling regimes.

    - `Rejection_RandomSampling.R`: R script to perform ABC using a simple rejection approach for different sampling regimes.

    - `LocLinear_RandomSampling.R`: R script to perform ABC using a local linear regression approach for different sampling regimes.

    - `NN_RandomSampling.R`: R script to perform ABC using a non-linear regression correction approach using neural networks for different sampling regimes.

    - `Emp_NN_RandomSampling_eta.R`: R script to perform ABC with empirical data using a non-linear regression correction approach using neural networks.


### **Humans Folder**

    - `Rejection_Wassertstein_mu12e07.R`: R script to perform ABC with empirical data using a simple rejection approach on Wasserstein distances.

    - `LCT_AfricaRange_mu12e07.R`: R script to perform ABC with empirical data using a non-linear regression correction approach using neural networks.
   
---
