# Theo AO Rashid -- March 2020

# ----- NIMBLE website examples -----
# Disease mapping CAR model example with
# a BYM model

library(nimble)

# ----- Disease Mapping example -----
# HES for a respiratory disease in 2010
# 134 Intermediate Geographies near Glasgow

library(CARBayesdata, quietly = TRUE)
library(sp, quietly = TRUE)
library(spdep, quietly = TRUE)

# ----- Unpack the data -----

data(GGHB.IG)
data(respiratorydata)

respiratorydata_spatial <- merge(x = GGHB.IG, y = respiratorydata, by = "IG", all.x = FALSE)
W.nb <- poly2nb(respiratorydata_spatial, row.names =  rownames(respiratorydata_spatial@data))
# Determine neighbourhood/adjacency information needed for neighbourhood-based CAR model
nbInfo <- nb2WB(W.nb)

# A vector of indices indicating which regions are neighbours of which.
# nbInfo$adj
# A vector of weights. In this case, all weights are 1.
# nbInfo$weights

# ----- Build the model -----
# For the CAR model, we need to specify:
# adj -- vectors indicating neighbouring regions
# weights -- weights of each pair of neighbours
# num -- number of neighbours for each region

nregions <- nrow(respiratorydata_spatial)

code <- nimbleCode({
    # priors
    beta ~ dnorm(0, sd = 100)
    sigma_u ~ dunif(0,2)
    tau_u <- pow(sigma_u,-2)
    sigma_v ~ dunif(0,2)
    tau_v <- pow(sigma_v,-2)
    # latent spatial process using BYM
    u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau_u, zero_mean = 0)
    for (j in 1:N) {
        v[j] ~ dnorm(u[j], tau_v) # centred on CAR term
    }
    # likelihood
    for(i in 1:N) {
        # Poisson includes an offset for the expected count
        log(lambda[i]) <- log(expected[i]) + beta*x[i] + v[i]
        y[i] ~ dpois(lambda[i])
    }
})

# Covariate -- number of people income deprived in spatial unit
x <- respiratorydata_spatial$incomedep
x <- x - mean(x)  # centre for improved MCMC performance

set.seed(1)

constants <- list(N = nregions, L = length(nbInfo$adj), 
               adj = nbInfo$adj, weights = nbInfo$weights, num = nbInfo$num,
                        x = x, expected = respiratorydata_spatial$expected)
data <- list(y = respiratorydata_spatial$observed)
inits <- list(beta = 0, sigma_u = 1, sigma_v = 1, u = rnorm(nregions))

# ----- Create the model -----
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# ----- Compile the model -----
cModel <- compileNimble(model)

# ---- MCMC integration -----
conf <- configureMCMC(model, monitors = c('beta', 'sigma_u', 'sigma_v', 'lambda'))
# conf$printSamplers()

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

mcmc.out <- runMCMC(cMCMC, niter = 10000, nburnin = 1000, summary = TRUE, samples = TRUE)

# Density plots of posterior samples
# plot(density(mcmc.out$samples[,"sigma"]))
# plot(density(mcmc.out$samples[,"lambda[1]"]))
# plot(density(mcmc.out$samples[,"beta"]))