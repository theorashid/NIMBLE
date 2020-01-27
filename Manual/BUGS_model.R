# Theo AO Rashid -- January 2020

# ----- NIMBLE website examples -----
# Building a model from BUGS code in R using NIMBLE

library(nimble)
library(igraph)

# ----- Dyes Example -----
# Normal hierarchical model
dyesCode <- nimbleCode({
# Model
   for (i in 1:BATCHES) {
      for (j in 1:SAMPLES) {
         y[i,j] ~ dnorm(mu[i], sd = sigma.within);
      }
      mu[i] ~ dnorm(theta, sd = sigma.between);
   }
   
# Priors
   theta ~ dnorm(0.0, 1.0E-10);
   sigma.within ~ dunif(0, 100)
   sigma.between ~ dunif(0, 100)
})

# ----- Create the model -----
# the model is an R object
dyesModel <- nimbleModel(dyesCode, constants = list(BATCHES = 6, SAMPLES = 5))

# set data values
data <- matrix(c(1545, 1540, 1595, 1445, 1595, 1520, 1440, 1555, 1550, 
1440, 1630, 1455, 1440, 1490, 1605, 1595, 1515, 1450, 1520, 1560, 
1510, 1465, 1635, 1480, 1580, 1495, 1560, 1545, 1625, 1445), nrow = 6)

dyesModel$setData(list(y = data))
# dyesModel$y

# ----- Use the model -----
# set values
dyesModel$theta <- 1500
dyesModel$mu <- rnorm(6, 1500, 50)
dyesModel$sigma.within <- 20
dyesModel$sigma.between <- 20

# get values
# dyesModel$y[1,]
# dyesModel$theta

# calculate log probabilities
# dyesModel$calculate(c('theta', 'mu[1:6]', 'y[,2]'))

# simulate the model (or part of)
dyesModel$simulate(c('mu[1:3]'))
# dyesModel$mu

# query the model's relationships
# dyesModel$getDependencies(c('theta', 'mu[3]'))

# plot the model graph
# plot(dyesModel$getGraph()) # not useful when complex

# ----- Compile the model -----
# The complied model can be used in the same way, but runs much faster
compiled_dyesModel <- compileNimble(dyesModel)
compiled_dyesModel$theta <- 1450
compiled_dyesModel$calculate() ## all nodes by default
