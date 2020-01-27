# Theo AO Rashid -- January 2020

# ----- Linear model -----
# Simple model with only "common" terms

# Numbers of each quantity
n_age_groups <- 19
n_sex        <- 2
n_LSOA       <- 4833
n_year       <- 17 # 2001 - 2017

code <- nimbleCode({
    # Priors
    alpha ~ dnorm(0, 0.00001)
    beta  ~ dnorm(0, 0.00001)

    # Likelihood
    for (i in 1:N) {
        y[i] ~ dpois(mu[i])
        log(mu[i]) <- log(n[i]) + epsilon[age[i], LSOA[i], yr[i]]
    }
})

constants <- list(N = nregions, L = length(nbInfo$adj), 
               adj = nbInfo$adj, weights = nbInfo$weights, num = nbInfo$num,
                        x = x, expected = respiratorydata_spatial$expected)
data <- list(y = respiratorydata_spatial$observed)
inits <- list(beta = 0, sigma = 1, s = rnorm(nregions))

model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Compile model into C++
# Cmodel <- compileNimble(model)

# # Run the MCMC
# MCMC <- buildMCMC(Cmodel)
# CMCMC <- compileNimble(MCMC)