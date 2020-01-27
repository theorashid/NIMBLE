# Theo AO Rashid -- January 2020

# ----- NIMBLE website examples -----
# Parallelisation with NIMBLE

# The key consideration is to ensure that all NIMBLE execution,
# including model building, is conducted inside the parallelized code.

library(parallel)

# ----- ----- -----
this_cluster <- makeCluster(4) # cluster with 4 cores

# Set up a situation to run 4 independent MCMC simulations
# Create a function that includes all modelling steps
# Run that function in parallel

set.seed(10120)
# Simulate some data
myData <- rgamma(1000, shape = 0.4, rate = 0.8)

# Create a function with all the needed code
run_MCMC_allcode <- function(seed, data) {
  library(nimble)
  
  # myCode could be specified outside the parallelised code
  myCode <- nimbleCode({
    a ~ dunif(0, 100)
    b ~ dnorm(0, 100)
  
    for (i in 1:length_y) {
      y[i] ~ dgamma(shape = a, rate = b)
    }
  })
  
  # Create the model from the BUGS code
  myModel <- nimbleModel(code = myCode,
                          data = list(y = data),
                          constants = list(length_y = 1000),
                          inits = list(a = 0.5, b = 0.5))
  
  # Compile the model into C++ code
  CmyModel <- compileNimble(myModel)
  
  # Run the MCMC
  myMCMC <- buildMCMC(CmyModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, niter = 10000, setSeed = seed)
  
  return(results)
}

# execute the code using parallel apply function
chain_output <- parLapply(cl = this_cluster, X = 1:4, 
                          fun = run_MCMC_allcode, 
                          data = myData)

# It's good practice to close the cluster when you're done with it.
stopCluster(this_cluster)

# Results of the 4 processes
par(mfrow = c(2,2))
for (i in 1:4) {
  this_output <- chain_output[[i]]
  plot(this_output[,"b"], type = "l", ylab = 'b')
}