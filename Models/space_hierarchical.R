# Theo AO Rashid -- February 2020

# ----- Hierarchical model -----
# Model making use of the hierarchical geographies
# LSOA -> MSOA -> LAD

# Currently, the model only accounts for LSOA -> LAD

# Area-specific linear trend in time
# Two levels of hierarchy (LSOA, LAD)
# No age structure

# Numbers of each quantity
N_sex        <- 2
N_age_groups <- 19
N_LSOA       <- 4835
N_LAD        <- 33
N_year       <- 17 # 2001 - 2017

# Indices:
# - a -- age
# - j -- space, each LSOA
# - t -- year (time)

code <- nimbleCode({
	# PRIORS

    # General terms
	alpha0 ~ dnorm(0, 0.00001)
	beta0 ~ dnorm(0, 0.00001)

	# Intermediate area (LAD) random effects for intercepts and slopes
	for(j2 in 1:N_LAD){
		alpha2[j2] ~ dnorm(alpha0, tau_alpha2)
		beta2[j2]  ~ dnorm(beta0, tau_beta2)
	}
	#need a variance for both sets of random effects 
	sigma_alpha2 ~ dunif(0,2)
	tau_alpha2 <- pow(sigma_alpha2,-2)
	sigma_beta2  ~ dunif(0,2)
	tau_beta2  <- pow(sigma_beta2,-2)

	# Area (LSOA) specific random effects for intercepts and slopes
	for(j in 1:N_LSOA){
        # dist_lookup gets the LAD corresponding to the LSOA
		alpha1[j] ~ dnorm(alpha2[dist_lookup[j]], tau_alpha1) # centred on LAD term
		beta1[j]  ~ dnorm(beta2[dist_lookup[j]], tau_beta1)
	}
	#need a variance for both sets of random effects 
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