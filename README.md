# R functions for Semiparametric Partial Common PrincipalComponent Analysis (PCPCA) for Covariance Matrices

## Overview

This repo stores R functions that can estimate the common principle components (CPCs) in PCPCA with or without knowing the number of CPCs.

For details, please see the manuscript [here](https://www.biorxiv.org/content/10.1101/808527v1.abstract). Currently we are under revision and the functions in this repo are subject to change.

## Functions

- __R/pcpc_est.R__ A function that can estimate CPCs with or without knowing the number of CPCs. When the number of CPCs is unknown, this function calls __R/pcpc_test.R__ to estimate it. 

- We caution that, when the dimension of matrices is large (e.g., more than 100), estimating the number of CPCs takes __long__ time. In this case, we recommend parallel programming to reduce the run time. We give an example of parallism in __data-analysis/all-region-analysis.R__.

## Code for data application

The R scripts under the "data-analysis" folder contain the code of data application in the [paper](https://arxiv.org/abs/1910.13954).

## Code for simulation studies

The R scripts under the "simulations" folder contain the code of simulations in the [paper](https://arxiv.org/abs/1910.13954).

