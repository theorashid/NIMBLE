# Theo AO Rashid -- February 2020

# ----- Space-Age model -----
# Introduces age structure using random walk

# Area-specific linear trend in time
# Two levels of hierarchy (LSOA, LAD)
# RW over ages for intercepts and slopes

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
	sigma_alpha1 ~ dunif(0,2)
	tau_alpha1 <- pow(sigma_alpha1,-2)
	sigma_beta1  ~ dunif(0,2)
	tau_beta1  <- pow(sigma_beta1,-2)

    # Age effects for intercepts and slopes
	alpha_age[1] <- 0 # initialise first terms for RW
	beta_age[1]  <- 0
	for(a in 2:N_age_groups){
		alpha_age[a] ~ dnorm(alpha_age[a-1], tau_alpha_age) # RW based on previous age group
		beta_age[a]  ~ dnorm(beta_age[a-1], tau_beta_age)
	}
	sigma_alpha_age ~ dunif(0,2)
	tau_alpha_age <- pow(sigma_alpha_age,-2)
	sigma_beta_age ~ dunif(0,2)
	tau_beta_age <- pow(sigma_beta_age,-2)


    # Put all parameters together into indexed lograte
	for(a in 1:N_age_groups){
		for(j in 1:N_LSOA){
			for(t in 1:N_year){
				lograte[a, j, t] <- (alpha1[j] + alpha_age[a]) + (beta1[j] + beta_age[a]) * (t - 8)
            }
        }
    }

	# LIKELIHOOD
	for(i in 1:N){
		y[i] ~ dpois(mu[i])
		log(mu[i]) <- log(n[i]) + lograte[age[i], lsoa[i], yr[i]]		
	}
})