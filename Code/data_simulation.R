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
#### Last Update : 2/24/2021                         ####
#########################################################
#########################################################

rm(list=ls())

## Packages
library(survival)
library(nlme)
library(JM)

# Simulation settings
set.seed(7721)
N <- 250 #number of simulated individuals 
K <- 10 # sample(1:10, N, replace = T) # number of planned repeated measurements per subject, per outcome
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
id <- rep(1:N, each = K)
  
# 2. Generate the times at which inidividuals have their longitudinal measurements
# Fixed measurement times: 
# times <- c(replicate(N, 0:(K-1)))
# Random measurement times: 
# times <- c(replicate(N, c(0, sort(runif(K-1, 0, t.max))))) #uniform 
# Can you do it for a poisson process? (Hint: you might need the "cumsum" function)
msmt.times <- c(replicate(N, c(0, cumsum(rexp(K - 1, insp.rate)))))

# 3. Begin creating your data set
dat <- data.frame(ID = id, Time = msmt.times)

# 4. Create design matrices X and Z 
# X: design matrix for fixed effects (should have as many columns as elements in "betas")
# Z: design matrix for random effects (should have as many columns as elements in "A")
X <- matrix(cbind(rep(1, length(id)), msmt.times), ncol = 2)
Z <- matrix(cbind(rep(1, length(id)), msmt.times), ncol = 2)

# 5. Simulate random effects
# b: simulated random effects b using covariance matrix A 
b <- cbind(rnorm(N, 0, A[1]), rnorm(N, 0, A[4]))

# 6. Simulate longitudinal responses
index = 1
eta.y = vector()
for (i in 1:N){
  eta.y[c(index:(index+9))] <- X[c(index:(index+9)),] %*% betas + Z[c(index:(index+9)),] %*% b[i,] #Xi'(t)*beta+zi'(t)*bi
  index = index + 10
}
y <- rnorm(K*N, mean=eta.y, sd = sigma.y) #~Normal(mean=eta.y, sd=sigma.y) 

# 7. Add longitudinal responses to data set 
dat$y <- y 



# Survival submodel ---------------------------------------------------------------------------
# 1. Write a function of the hazard for each individual as a function of "s" (time) and "u" (a uniform
# random variable) 
# Integrate this function w.r.t "s" over the interval (0,t), where you want to solve for "t" 
invS <- function (t, u, i) {
  #inverse function used to obtain survival time 
  #t: time 
  #u: draw from a uniform[0,1] distribution
  #i: patient id indicator 
  h <- function (s) {
    #hazard function (based on survival and longitudinal submodels)
    XX <- cbind(1, s) %*% betas
    ZZ <- cbind(1, s) %*% b[i,]   
    m <- XX + ZZ
    haz <- h0*exp(alpha*m)
    return(haz)
  }
  integrate(h, lower = 0, upper = t)$value + log(u) # u is (1-u) but that is also uniform
}

# 2. Solve for survival times 
u <- runif(N) #generate random U[0,1] random variable 
trueTimes <- numeric(N) #initiate vector of true survival times
for (i in 1:N) {
  # Solve for true survival time (Hint: use "uniroot" function)
  # Wrap the uniroot function in "try" since for some individuals you will not return a valid 
  # survival time (we will consider those patients "cured" and assign them survival time of
  # infinity)
  surv <- function(t){ invS(t, u = u[i], i = i)}
  Root <- try(uniroot(surv, interval = c(0, 10))$root)
  # "Cure" patients have event time set to Infinity 
  trueTimes[i] <- ifelse(inherits(Root, "try-error"),Inf,Root)
}

# 3. Simulate censoring times 
Ctimes <- runif(N, min = 0, max = 10)
  
# 4. Compute observed survival time and event indicator (don't forget that there 
# is a max follow-up time)
Time <- ifelse(Ctimes < trueTimes, Ctimes, trueTimes)
event <- ifelse(Ctimes < trueTimes, 0,  1)

  
# 5. Add in survival time and indicator to data set 
dat$Survtime = rep(Time, each = K)
dat$event = rep(event, each = K)

# 6. Drop measurement times that occur after follow-up time 
# This should then be your longitudinal data set
dat <- dat[dat$Time <= dat$Survtime,]

# 7. Create a data set with one row for each individual (for survival submodel)
dat.id <- dat %>% group_by(ID) %>% slice(1)

# Plots of simulated data set -----------------------------------------------------------------
# Does your simulated data look correct? 
## Consider: KM curves 
# plot(Surv(Survtime,event)~1, data=dat)
# ## Consider: plot of the longitudinal biomarker 
# table(table(dat$ID))
# sum(dat.id$event)/N


# Fit models --------------------------------------------------------------
# 1. Fit joint model
# use method="weibull-PH-aGH" to decrease computation time 
# a. longitudinal submodel 
start_joint <- Sys.time()
lme.Fit <- lme(y ~ Time, random = ~ Time | ID, data = dat)

# b. survival submodel
surv.Fit <- coxph(Surv(Survtime, event) ~ 1, data = dat.id, x = TRUE)

# c. Joint model with piecewise-constant baseline hazard with 
# d. Adaptive GH numeric integration
joint.Fit <- jointModel(lme.Fit, surv.Fit, timeVar = "Time", method = "weibull-PH-aGH")
end_joint <- Sys.time()
time_joint <- as.numeric(end_joint - start_joint)


# 2. Fit two-stage model 
# merge data to get start and stop times
start_stop <- tmerge(dat.id, dat.id, ID, death = event(Survtime, event))

# merge with original data to get longitudinal start and stop times at each visit for an individual
start_stop <- tmerge(start_stop, dat, ID, y = tdc(Time, y))

# lme fit with start and stop times 
start_2stage <- Sys.time()
lme.Fit2 <- lme(y ~  tstart, random = ~ tstart | ID, data = start_stop)

# add predictions to the data 
start_stop$pred.y <- c(predict(lme.Fit2))

# Two-Stage model with time-varying sqrt aortic gradient
cox.Fit <- coxph(Surv(tstart, tstop, death) ~ pred.y, data = start_stop)
end_2stage <- Sys.time()
time_2stage <- as.numeric(end_2stage - start_2stage)

# Save appropriate coefficient estimates
# Check that these are close to the true value before you run for 500 runs! 

# Joint model 
joint.Fit$coefficients$betas
joint.Fit$coefficients$alpha #Association
joint.Fit$coefficients$sigma
sqrt(diag(joint.Fit$coefficients$D)) #varcov

# Two-Stage model 
coef(summary(lme.Fit2))[,1]
coef(cox.Fit)
summary(lme.Fit2)$sigma
VarCorr(lme.Fit2)[c(1,2),2]

# Compute summary metrics  ------------------------------------------------

# Compute the appropriate metrics 

# Create tables/figures to appropriately summarize the results 




  
                                                   
                                                   