# Theo AO Rashid -- January 2020

# ----- NIMBLE User Manual -----
# Pump failure model -- lightning introduction for Chapter 2

library(nimble)
library(igraph)

# ----- 2.2 Creating the model -----

# time duration length -- t[i]
# pump number -- i
pumpConsts <- list(N = 10,
t = c(94.3, 15.7, 62.9, 126, 5.24,
31.4, 1.05, 1.05, 2.1, 10.5))

# number of failures for pump -- x[i]
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

# failure rate -- theta
# Goal is to estimate the value of parameters alpha and beta
# which are used to specify the arrival time -- lambda
pumpCode <- nimbleCode(
  {
for (i in 1:N){
  theta[i] ~ dgamma(alpha,beta)
  lambda[i] <- theta[i]*t[i]
  x[i] ~ dpois(lambda[i])
  }
alpha ~ dexp(1.0)
beta ~ dgamma(0.1,1.0)
})

pumpInits <- list(alpha = 1, beta = 1,
theta = rep(0.1, pumpConsts$N))

# Create the model object using nimbleModel()
pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

# model$getNodeNames() lets us look at the nodes that comprise the model's
# directed acyclic graph (DAG), and model$plotGraph() lets us plot this DAG

# To look at individual nodes, use the syntax model$nodename
# e.g pump$x, pump$theta, pump$lambda

# Show all dependencies of alpha and beta terminating in stochastic nodes
# pump$getDependencies(c("alpha", "beta"))
# and only deterministic dependencies
# pump$getDependencies(c("alpha", "beta"), determOnly = TRUE)

# To get the log probabilities of the likelihood, pump$logProb_x
# and pump$getLogProb("x") for the sum

# Simulate new theta values
set.seed(1) # This makes the simulations here reproducible
pump$simulate("theta")
# pump$theta # shows the new theta values

# Now calculate the new values of the dependencies of theta
pump$calculate(pump$getDependencies(c("theta")))
# pump$lambda # shows the new lambda values
# pump$logProb_x # shows the new log probabilities of the likelihood

# ----- 2.3 Compiling the model -----
Cpump <- compileNimble(pump) # complies in C++
# Cpump$theta
# the compiled model needs to be used before any compilation of algorithms

# ----- 2.4 One-line invocation of MCMC -----
# default NIMBLE MCMC algorithm is used to generate posterior samples
mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                       data = pumpData, inits = pumpInits,
                       monitors=c("alpha","beta","theta","lambda"),
                       nchains = 2, niter = 10000,
                       summary = TRUE, WAIC = TRUE)

# names(mcmc.out) # output, each called by mcmc.out$name

# --- 2.5 Creating, compiling and running a basic MCMC configuration ---
pumpConf <- configureMCMC(pump, print = TRUE) # create an MCMC configuration
pumpConf$addMonitors(c("alpha", "beta", "theta"))

pumpMCMC <- buildMCMC(pumpConf)
CpumpMCMC <- compileNimble(pumpMCMC, project = pump) # compile into C++

niter <- 1000
set.seed(1)
samples <- runMCMC(CpumpMCMC, niter = niter) # run the MCMC sampler

# ----- 2.6 Customizing the MCMC -----
# Add an adaptive block sampler on alpha and beta jointly
pumpConf$addSampler(target = c("alpha","beta"), type = "RW_block",
                    control = list(adaptInterval = 100))

pumpMCMC2 <- buildMCMC(pumpConf)
# compile into C++ (need to reset nimbleFunctions to add new MCMC)
CpumpNewMCMC <- compileNimble(pumpMCMC2, project = pump,
                              resetFunctions = TRUE)

set.seed(1)
CpumpNewMCMC$run(niter)
# samples <- as.matrix(CpumpNewMCMC$mvSamples)

# ----- 2.7 Running MCEM -----
# # Maximise the maginal likelihood for alpha and beta
# pump2 <- pump$newModel()

# box = list(list(c("alpha","beta"), c(0, Inf)))

# pumpMCEM <- buildMCMC(model = pump2, latentNodes = "theta[1:10]",
#                       boxConstraints = box)
# pumpMLE <- pumpMCEM$run()
# # pumpMLE # return very accurate values of alpha and beta