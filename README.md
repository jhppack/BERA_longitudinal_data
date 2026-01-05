@ -0,0 +1,40 @@
# Bayesian Extended Redundancy Analysis for Longitudinal Data

This repository contains the R code and JAGS model implementation for the **Bayesian Extended Redundancy Analysis (ERA)** adapted for **longitudinal (panel) data**.

The code provided here **demonstrates** the simulation design and model estimation procedures presented in the manuscript tentatively titled:
> *"Bayesian Extended Redundancy Analysis for Longitudinal Data"*

## Overview

This project extends the original Bayesian ERA model (Park et al., 2019) to handle longitudinal data structures. The model is designed to analyze relationships between high-dimensional predictors and multiple outcomes over time, featuring:

* **Time-varying intercepts** to capture temporal trends.
* **Time-invariant regression coefficients** for component stability.
* **Subject-specific random effects** to account for heterogeneity.
* **Factor analysis structure** for residual covariance at each time point.
* **Handling of missing data** under the Missing At Random (MAR) assumption.

## Repository Contents

* `BERA_Panel_model.txt`: The JAGS model specification file defining the likelihood, priors, and structural equations.
* `simulation_script.R`: An R script that:
    1.  Generates synthetic longitudinal data based on the simulation design.
    2.  Prepares the data for JAGS.
    3.  Fits the model using the `R2jags` interface.
    4.  Summarizes posterior estimates.

## Prerequisites

To run the code, you need to have **R** and **JAGS** installed on your system.

### 1. Install JAGS
This project requires **JAGS (Just Another Gibbs Sampler)**.
* Download and install it from [SourceForge](https://sourceforge.net/projects/mcmc-jags/).

### 2. Install R Packages
The following R packages are required:

```r
install.packages("R2jags")
install.packages("mvtnorm")
