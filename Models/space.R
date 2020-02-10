# Theo AO Rashid -- February 2020

# ----- Spatial model -----
# Area-specific linear trend in time
# One level of hierarchy
# No age structure

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

    # General terms
	alpha0 ~ dnorm(0, 0.00001)
	beta0  ~ dnorm(0, 0.00001)

	# Area specific random effects for intercepts and slopes
	# When building up hierarchical models, for better mixing centre the lower 
    # level parameter on the higher level one rather than adding e.g. NOT alpha0+alpha[j]
	for(j in 1:N_LSOA){
		alpha1[j] ~ dnorm(alpha0, tau_alpha1) # centred on general terms
		beta1[j]  ~ dnorm(beta0, tau_beta1)
	}
	# Need a variance for both sets of random effects 
	sigma_alpha1 ~ dunif(0,2)
	tau_alpha1 <- pow(sigma_alpha1,-2)
	sigma_beta1  ~ dunif(0,2)
	tau_beta1  <- pow(sigma_beta1,-2)

    # Put all parameters together into indexed lograte
    for(a in 1:N_age_groups){
		for(j in 1:N_LSOA){
			for(t in 1:N_year){
				lograte[a, j, t] <- alpha1[j] + beta1[j] * (t - 8)
            }
        }
    }

	# LIKELIHOOD
	for(i in 1:N){
		y[i] ~ dpois(mu[i])
		log(mu[i]) <- log(n[i]) + lograte[age[i], lsoa[i], yr[i]]		
	}
})