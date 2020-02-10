# Theo AO Rashid -- February 2020

# ----- Linear model -----
# Simple model with only "common" terms

# Numbers of each quantity
N_sex        <- 2
N_age_groups <- 19
N_LSOA       <- 4835
N_year       <- 17 # 2001 - 2017

# Indices:
# - a -- age
# - j -- space, each LSOA
# - t -- year (time)

code <- nimbleCode({
    # PRIORS
	#dnorm is mean, precision
	alpha0 ~ dnorm(0, 0.00001)
	beta0  ~ dnorm(0, 0.00001)

    # Put all parameters together into indexed lograte
	for(a in 1:N_age_groups){
		for(j in 1:N_LSOA){
			for(t in 1:N_year){
                # remember to centre your t to improve sampling performance
                # (cuts down correlation between parameters)
                # DO NOT make it a dynamic centrering
				lograte[a, j, t] <- alpha0 + beta0 * (t - 8)
            }
        }
    }

    # LIKELIHOOD
    # N total number of cells, i.e. ages*years*areas(*sex)
    for (i in 1:N) {
        # y is number of deaths in that cell, assumed Poisson 
		# mu is predicted number of deaths in that cell
        y[i] ~ dpois(mu[i])
        log(mu[i]) <- log(n[i]) + lograte[age[i], LSOA[i], yr[i]]
    }
})

# need initial estimates for alpha0, beta0

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