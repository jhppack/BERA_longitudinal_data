## -----------------------------------------------------------------------------
## Simulation Study: Bayesian Extended Redundancy Analysis for Longitudinal Data
## -----------------------------------------------------------------------------
## Description:
##   This R code demonstrates the simulation study presented in the manuscript 
##   tentatively titled "Bayesian Extended Redundancy Analysis for Longitudinal 
##   Data". The model estimation is implemented using JAGS.
## -----------------------------------------------------------------------------

rm(list=ls())
library(R2jags)
library(mvtnorm)

# Set working directory to your project path
setwd("D:/Project/BERA/Panel/")

################################################################################
#### 1. Data Generation (Simulation)
################################################################################

## -----------------------------------------------------------------------------
## Simulation Settings
## -----------------------------------------------------------------------------
ns <- 100             # Sample size (N)
K  <- 2               # Number of latent components
Q  <- 3               # Number of outcome variables
T  <- 3               # Number of time points (waves)
nX <- 4               # Number of predictors (P)
prop.missing <- 0.1   # Proportion of missing values in outcomes

## -----------------------------------------------------------------------------
## True Parameter Specification
## -----------------------------------------------------------------------------

# 1. Component Weights (Time-invariant)
#    - Maps predictors to components
tr.W <- matrix(c(1, 0.5, 0, 0, 0, 0, 0.5, 0.5), nX, K)

# 2. Component Regression Coefficients (Time-invariant)
#    - Effect of components on outcomes (Slopes)
tr.A <- matrix(c(1, 1.2, 2, 0, 0, 2), K, Q)

# 3. Intercepts (Time-varying)
#    - Rows: Outcomes (Q), Columns: Time points (T)
tr.bt0 <- cbind(c(1, -1, 0), c(0, 0, 0), c(1, 1, 0))

# 4. Residual Covariance Matrix for Outcomes (Time-invariant)
tr.Sig <- matrix(0.3, Q, Q)
diag(tr.Sig) <- c(1.2, 0.5, 1)

# 5. Covariance Matrix for Predictors
tr.Sig.Xs <- matrix(c(1, 0.3, 0.1, 0.1, 
                      0.3, 1, 0.1, 0.1, 
                      0.1, 0.1, 2, 0.3, 
                      0.1, 0.1, 0.3, 2), nX, nX)

# 6. Covariance Matrix for Random Effects
#    - Random effects account for subject-specific heterogeneity
n_re <- nrow(tr.bt0) + length(tr.A) # 3 + 6 = 9
tr.Sig.u <- matrix(1, n_re, n_re) * 0.1
diag(tr.Sig.u) <- 0.5


## -----------------------------------------------------------------------------
## Data Generation Loop
## -----------------------------------------------------------------------------

## Generate Predictors (X): Array [N, nX, T]
Xs <- array(NA, dim = c(ns, nX, T)) 
for (tm in 1:T) {
  Xs[, , tm] <- rmvnorm(ns, mean = rep(0, nX), sigma = tr.Sig.Xs)
}

## Generate Outcomes (Y): Array [N, Q, T]
ys <- array(0, dim = c(ns, Q, T))

for (i in 1:ns) {
  # Sample subject-specific random effects
  ui <- rmvnorm(1, mean = rep(0, nrow(tr.Sig.u)), sigma = tr.Sig.u)
  
  # Parse random effects:
  # - First Q elements: Random intercepts for each outcome
  # - Remaining elements: Random slopes for components (reshaped to K x Q)
  ui_intercepts <- as.vector(ui[1:Q])
  ui_slopes     <- matrix(ui[-(1:Q)], K, Q) 
  
  for (tm in 1:T) {
    # Mean Structure:
    # Intercept(t) + Random_Intercept + (X(t) * W) * (Slope + Random_Slope)
    mu_i_tm <- tr.bt0[, tm] + 
               ui_intercepts + 
               Xs[i, , tm] %*% tr.W %*% (tr.A + ui_slopes)
    
    # Sample outcome from multivariate normal
    ys[i, , tm] <- rmvnorm(1, mean = mu_i_tm, sigma = tr.Sig)
  }
}

## Introduce Missing Values (MAR/MCAR)
for (j in 1:2) { 
  # Select random indices for missingness
  idx_missing <- sample.int(ns)[1:floor(ns * prop.missing)]
  ys[idx_missing, j, T] <- NA
}


################################################################################
#### 2. Bayesian Analysis (JAGS Fitting)
################################################################################

## -----------------------------------------------------------------------------
## Model Setup & Hyperparameters
## -----------------------------------------------------------------------------

# Indicator vectors mapping predictors to components
ind.W1 <- c(1:2)
ind.W2 <- c(3:4)

# Hyperparameters for Priors
m_b <- rep(0, T + K)      # Prior mean for regression coefs
S_b <- diag(T + K)        # Prior covariance for regression coefs
m_w <- rep(0, T + K)      # Prior mean for weights/loadings
S_w <- diag(T + K)        # Prior covariance for weights/loadings

## -----------------------------------------------------------------------------
## Run JAGS
## -----------------------------------------------------------------------------

# Initialize Model
jags.fit <- jags.model(
  file = "BERA_longitudinal_model.txt",
  data = list(
    Y = ys, 
    X = Xs, 
    N = ns, 
    Q = dim(ys)[2], 
    K = K, 
    T = T,
    m_b = m_b, S_b = S_b, 
    m_w = m_w, S_w = S_w, 
    ind.W1 = ind.W1, 
    ind.W2 = ind.W2, 
    n.W1 = length(ind.W1), 
    n.W2 = length(ind.W2), 
    W.a = c(1, 0.5) 
  ),
  n.chains = 3, 
  n.adapt = 100
)

# Parameters to monitor
parm.list <- c("Wss", "Betas", "sig2.u", "Cov.y", "sig2.y")

# Sampling (Posterior Generation)
BERA_fit <- coda.samples(
  jags.fit, 
  variable.names = parm.list, 
  n.iter = 50000, 
  thin = 5
)

# Post-processing (Burn-in & Thinning)
BERA_fit <- window(BERA_fit, burnin = 10000, thin = 5)

## -----------------------------------------------------------------------------
## Summary of Results
## -----------------------------------------------------------------------------

fit_summary <- summary(BERA_fit)
print(cbind(fit_summary[[1]][, 1:2], fit_summary[[2]][, c(1, 5)]))