#########################################################
#########################################################
####    BIOS 7721 Final Project : Simulation Study   ####
####                                                 ####
#### Author : Samantha Bothwell                      ####
####                                                 ####
#### Professor : Krithika Suresh                     ####
####                                                 ####
#### File Purpose : Simulate longitudinal marker     ####
####    values from a given distribution and create  ####
####    random visit times from a Poisson distibution####
####                                                 ####
#### Last Update : 2/12/2021                         ####
#########################################################
#########################################################

rm(list=ls())

## Packages 


# Simulation settings
set.seed(7721)
N <- 250 #number of simulated individuals 
K <- sample(1:10, N, replace = T) # number of planned repeated measurements per subject, per outcome
cens_horiz <- 10 #maximum follow-up time
insp.rate <- 1 #rate at which inspection times are made

  
# True parameter values ---------------------------------------------------
betas <- c(-0.2, 0.3) #betas for longitudinal model
sigma.y <-  0.5 #measurement error standard deviation
A <- matrix(c(0.5, 0, 0, 0.3), ncol = 2) #covariance matrix for random effects
alpha <- 0.6 #association parameter
h0 <- 0.1 #constant baseline hazard

  
# Longitudinal submodel ---------------------------------------------------
# 1. Generate the id variable (each individual has K measurements)
id <- rep(1:N, K)
  
# 2. Generate the times at which inidividuals have their longitudinal measurements
# Fixed measurement times: 
# times <- c(replicate(N, 0:(K-1)))
# Random measurement times: 
# times <- c(replicate(N, c(0, sort(runif(K-1, 0, t.max))))) #uniform 
# Can you do it for a poisson process? (Hint: you might need the "cumsum" function)
msmt.times <- c(replicate(N, c(0, cumsum(rpois(K-1, 1)))))
msmt.times <- lapply(K, function(i){c(0, cumsum(rpois(K-1, 1)))})

# 3. Begin creating your data set
dat <- 

# 4. Create design matrices X and Z 
# X: design matrix for fixed effects (should have as many columns as elements in "betas")
# Z: design  matrix for random effects (should have as many columns as elements in "A")
X <- 
Z <- 

# 5. Simulate random effects
# b: simulated random effects b using covariance matrix A 
b <- 

# 6. Simulate longitudinal responses
eta.y <- #Xi'(t)*beta+zi'(t)*bi
y <- #~Normal(mean=eta.y, sd=sigma.y) 

# 7. Add longitudinal responses to data set 
                                                   
                                                   
                                                   
                                                   