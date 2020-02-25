# Theo AO Rashid -- February 2020

#read in the function from wherever you stored it
source("period_life_table.R")
 
##EXAMPLE
#rough estimated death rates for 18 age groups
mxE <- exp(-5.5 + c(-1,-3,-2.9,-1.8,-1.7,-1.7,-1.6,-1.4,-.8,-.2,.3,1.0,1.5,1.9,2.4,2.8,3.4,3.8))
# Pops <- c(rep(100000,10),rep(1000,8))
# rpos(n=18, mxE*Pops)
 
#function to give table
test <- PeriodLifeTable(mx=mxE,age=seq(0, 85, 5),ax=rep(NA,18),full.table=TRUE)
 # ex column and for life expectancy at birth you want the first number (e0).
 # If you wanted life expectancy age 5 (i.e. assuming that an individual
 # has survived the first age group) then you would use the
 # second number in the ex colum
 