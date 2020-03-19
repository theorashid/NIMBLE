# Theo AO Rashid -- March 2020

# ----- NIMBLE website examples -----
# Disease mapping CAR model example with
# a BYM2 model (Riebler et al. 2016)

library(nimble)
library(INLA)

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
# Determine neighbourhood/adjacency
nbInfo <- nb2WB(W.nb)

# Find scale parameter of Q for BYM2
nb2INLA("W.graph", W.nb) # create the adjacency matrix in INLA format
W.adj <- paste(getwd(),"/W.graph",sep="") # name the object

# usage from help(inla.scale.model)
g <- inla.read.graph(filename = "W.graph")
print(paste("Number of connected components", g$cc$n))
Q = -inla.graph2matrix(g)
diag(Q) = 0
diag(Q) = -rowSums(Q)
n = dim(Q)[1]
Q.scaled = inla.scale.model(Q, constr = list(A = matrix(1, 1, n), e=0))

scale = Q.scaled[1,1]/Q[1,1] # scale is 0.4340388

# ----- Build the model -----
# For the CAR model, we need to specify:
# adj -- vectors indicating neighbouring regions
# weights -- weights of each pair of neighbours
# num -- number of neighbours for each region

nregions <- nrow(respiratorydata_spatial)

code <- nimbleCode({
    # priors
    beta ~ dnorm(0, sd = 100)
    sigma ~ dunif(0,2)
    rho ~ dbeta(1,1)

    # latent spatial process using BYM2
    u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean = 1) # structured
    for (j in 1:N) {
        v[j] ~ dnorm(0,1) # unstructured
        Corr[j] <- sigma * u[j] * sqrt(rho/scale)
        UCorr[j] <- sigma * v[j] * sqrt(1 - rho)
        bym2[j] <- Corr[j] + UCorr[j]
    }

    # likelihood
    for(i in 1:N) {
        # Poisson includes an offset for the expected count
        log(lambda[i]) <- log(expected[i]) + beta*x[i] + bym2[i]
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
inits <- list(beta = 0, sigma = 1, rho = 0.5, u = rnorm(nregions))

# ----- Create the model -----
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# ----- Compile the model -----
cModel <- compileNimble(model)

# ---- MCMC integration -----
conf <- configureMCMC(model, monitors = c('beta', 'sigma', 'rho', 'lambda'))
# conf$printSamplers()

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

mcmc.out <- runMCMC(cMCMC, niter = 10000, nburnin = 1000, summary = TRUE, samples = TRUE)

# Density plots of posterior samples
# plot(density(mcmc.out$samples[,"sigma"]))
# plot(density(mcmc.out$samples[,"lambda[1]"]))
# plot(density(mcmc.out$samples[,"beta"]))