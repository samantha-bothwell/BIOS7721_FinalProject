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


# Initialize dataframes 
sim = 500
dat.joint <- data.frame(matrix(NA, nrow = sim, ncol = 10))
colnames(dat.joint) <- c("Assoc", "Beta0", "Beta1", "Sigma", "A11", "A22", "SE0", "SE1", "SEy", "Time")
dat.2stage <- data.frame(matrix(NA, nrow = sim, ncol = 9))
colnames(dat.2stage) <- c("Assoc", "Beta0", "Beta1", "Sigma", "A11", "A22", "SE0", "SE1", "SEy", "Time")


# Simulation 
for (j in 1:sim){
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
    Root <- try(uniroot(surv, interval = c(0, 10), extendInt = "upX")$root)
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
  dat.joint[j, 1] <- joint.Fit$coefficients$alpha #Association
  dat.joint[j, 2] <- joint.Fit$coefficients$betas[1]
  dat.joint[j, 3] <- joint.Fit$coefficients$betas[2]
  dat.joint[j, 4] <- joint.Fit$coefficients$sigma
  dat.joint[j, 5] <- sqrt(diag(joint.Fit$coefficients$D))[1] #varcov
  dat.joint[j, 6] <- sqrt(diag(joint.Fit$coefficients$D))[2]
  
  dat.joint[j, 7] <- summary(joint.Fit)$`CoefTable-Long`[1,2] #asymptotic standard errors for long
  dat.joint[j, 8] <- summary(joint.Fit)$`CoefTable-Long`[2,2] #asymptotic standard errors for long
  dat.joint[j, 9] <- summary(joint.Fit)$`CoefTable-Event`[2,2] #asymptotic standard error for assoc
  
  dat.joint[j, 10] <- time_joint
  
  # Two-Stage model 
  dat.2stage[j, 1] <- coef(cox.Fit)
  dat.2stage[j, 2] <- coef(summary(lme.Fit2))[1,1]
  dat.2stage[j, 3] <- coef(summary(lme.Fit2))[2,1]
  dat.2stage[j, 4] <- summary(lme.Fit2)$sigma
  dat.2stage[j, 5] <- as.numeric(VarCorr(lme.Fit2)[1,2])
  dat.2stage[j, 6] <- as.numeric(VarCorr(lme.Fit2)[2,2])
  
  dat.2stage[j, 7] <- coef(summary(lme.Fit2))[1,2] #asymptotic standard errors for long
  dat.2stage[j, 8] <- coef(summary(lme.Fit2))[2,2] #asymptotic standard errors for long
  dat.2stage[j, 9] <- coef(summary(cox.Fit))[3] #asymptotic standard error for y (assoc)
  
  dat.2stage[j, 10] <- time_2stage
}
# Compute summary metrics  ------------------------------------------------

dat.2stage <- read.csv("D:/CU/Spring 2021/BIOS 7721/BIOS7721_FinalProject/Data/2stage.csv")[,-1]
dat.joint <- read.csv("D:/CU/Spring 2021/BIOS 7721/BIOS7721_FinalProject/Data/Joint.csv")[,-1]


# True parameter values ---------------------------------------------------
betas <- c(-0.2, 0.3) #betas for longitudinal model
sigma.y <-  0.5 #measurement error standard deviation
A <- matrix(c(0.5, 0, 0, 0.3), ncol = 2) #covariance matrix for random effects
alpha <- 0.6 #association parameter
h0 <- 0.1 #constant baseline hazard
sim <- 500


# Compute the appropriate metrics 
# true <- c(alpha, betas[1], betas[2], sigma.y, A[1], A[4])

## Empirical bias 
bias.2stage <- data.frame(cbind(y = mean(dat.2stage$Assoc) - alpha,
  b0 = mean(dat.2stage$Beta0) - betas[1],
  b1 = mean(dat.2stage$Beta1) - betas[2],
  sigma = mean(dat.2stage$Sigma) - sigma.y,
  A11 = mean(dat.2stage$A11) - A[1],
  A22 = mean(dat.2stage$A22) - A[4]))

bias.joint <- data.frame(cbind(y = mean(dat.joint$Assoc) - alpha,
  b0 = mean(dat.joint$Beta0) - betas[1],
  b1 = mean(dat.joint$Beta1) - betas[2],
  sigma = mean(dat.joint$Sigma) - sigma.y,
  A11 = mean(dat.joint$A11) - A[1],
  A22 = mean(dat.joint$A22) - A[4]))

## Asymptotic SE
asy.se.2stage <- data.frame(cbind(b0 = mean(dat.2stage$SE0),
  b1 = mean(dat.2stage$SE1),
  y = mean(dat.2stage$SEy)))

asy.se.joint <- data.frame(cbind(b0 = mean(dat.joint$SE0),
  b1 = mean(dat.joint$SE1),
  y = mean(dat.joint$SEy)))

## Empirical SE
emp.se.2stage <- data.frame(cbind(b0 = sd(dat.2stage$Beta0),
  b1 = sd(dat.2stage$Beta1),
  y = sd(dat.2stage$Assoc)), 
  sigma = sd(dat.2stage$Sigma), 
  A11 = sd(dat.2stage$A11), 
  A22 = sd(dat.2stage$A22))

emp.se.joint <- data.frame(cbind(b0 = sd(dat.joint$Beta0),
  b1 = sd(dat.joint$Beta1),
  y = sd(dat.joint$Assoc)), 
  sigma = sd(dat.joint$Sigma), 
  A11 = sd(dat.joint$A11), 
  A22 = sd(dat.joint$A22))


## Mean Square Error
mse.2stage <- data.frame(cbind(y = sum((dat.2stage$Assoc - alpha)^2)/sim,
  b0 = sum((dat.2stage$Beta0 - betas[1])^2)/sim,
  b1 = sum((dat.2stage$Beta1 - betas[2])^2)/sim,
  sigma = sum((dat.2stage$Sigma - sigma.y)^2)/sim,
  A11 = sum((dat.2stage$A11 - A[1])^2)/sim,
  A22 = sum((dat.2stage$A22 - A[2])^2)/sim))

mse.joint <- data.frame(cbind(y = sum((dat.joint$Assoc - alpha)^2)/sim,
  b0 = sum((dat.joint$Beta0 - betas[1])^2)/sim,
  b1 = sum((dat.joint$Beta1 - betas[2])^2)/sim,
  sigma = sum((dat.joint$Sigma - sigma.y)^2)/sim,
  A11 = sum((dat.joint$A11 - A[1])^2)/sim,
  A22 = sum((dat.joint$A22 - A[2])^2)/sim))

## Coverage rate 
conf.2stage <- data.frame(cbind(b0.low = dat.2stage$Beta0 - 1.96*dat.2stage$SE0,
  b0.up = dat.2stage$Beta0 + 1.96*dat.2stage$SE0,
  b1.low = dat.2stage$Beta1 - 1.96*dat.2stage$SE1,
  b1.up = dat.2stage$Beta1 + 1.96*dat.2stage$SE1,
  y.low = dat.2stage$Assoc - 1.96*dat.2stage$SEy,
  y.up = dat.2stage$Assoc + 1.96*dat.2stage$SEy))

conf.2stage$b0.cov = ifelse(betas[1] < conf.2stage$b0.low | betas[1] > conf.2stage$b0.up, 0, 1)
conf.2stage$b1.cov = ifelse(betas[2] < conf.2stage$b1.low | betas[2] > conf.2stage$b1.up, 0, 1)
conf.2stage$y.cov = ifelse(alpha < conf.2stage$y.low | alpha > conf.2stage$y.up, 0, 1)

conf.joint <- data.frame(b0.low = dat.joint$Beta0 - 1.96*dat.joint$SE0,
  b0.up = dat.joint$Beta0 + 1.96*dat.joint$SE0,
  b1.low = dat.joint$Beta1 - 1.96*dat.joint$SE1,
  b1.up = dat.joint$Beta1 + 1.96*dat.joint$SE1,
  y.low = dat.joint$Assoc - 1.96*dat.joint$SEy,
  y.up = dat.joint$Assoc + 1.96*dat.joint$SEy)

conf.joint$b0.cov = ifelse(betas[1] < conf.joint$b0.low | betas[1] > conf.joint$b0.up, 0, 1)
conf.joint$b1.cov = ifelse(betas[2] < conf.joint$b1.low | betas[2] > conf.joint$b1.up, 0, 1)
conf.joint$y.cov = ifelse(alpha < conf.joint$y.low | alpha > conf.joint$y.up, 0, 1)

cov.2stage <- data.frame(cbind(b0 = sum(conf.2stage$b0.cov)/sim, 
  b1 = sum(conf.2stage$b1.cov)/sim, 
  y = sum(conf.2stage$y.cov)/sim))

cov.joint <- data.frame(cbind(b0 = sum(conf.joint$b0.cov)/sim, 
  b1 = sum(conf.joint$b1.cov)/sim, 
  y = sum(conf.joint$y.cov)/sim))


## Average comp time 
time.2stage = mean(dat.2stage$Time)
time.joint = mean(dat.joint$Time)


### combine data 
sum.2stage = data.frame(cbind(y = c(bias.2stage$y, asy.se.2stage$y, emp.se.2stage$y, mse.2stage$y, cov.2stage$y),
  b0 = c(bias.2stage$b0, asy.se.2stage$b0, emp.se.2stage$b0, mse.2stage$b0, cov.2stage$b0), 
  b1 = c(bias.2stage$b1, asy.se.2stage$b1, emp.se.2stage$b1, mse.2stage$b1, cov.2stage$b1), 
  sigma = c(bias.2stage$sigma, NA, emp.se.2stage$sigma, mse.2stage$sigma, NA), 
  A11 = c(bias.2stage$A11, NA, emp.se.2stage$A11, mse.2stage$A11, NA), 
  A22 = c(bias.2stage$A22, NA, emp.se.2stage$A22, mse.2stage$A22, NA)))

sum.joint = data.frame(cbind(y = c(bias.joint$y, asy.se.joint$y, emp.se.joint$y, mse.joint$y, cov.joint$y),
  b0 = c(bias.joint$b0, asy.se.joint$b0, emp.se.joint$b0, mse.joint$b0, cov.joint$b0), 
  b1 = c(bias.joint$b1, asy.se.joint$b1, emp.se.joint$b1, mse.joint$b1, cov.joint$b1), 
  sigma = c(bias.joint$sigma, NA, emp.se.joint$sigma, mse.joint$sigma, NA), 
  A11 = c(bias.joint$A11, NA, emp.se.joint$A11, mse.joint$A11, NA), 
  A22 = c(bias.joint$A22, NA, emp.se.joint$A22, mse.joint$A22, NA)))

time = data.frame(cbind(`TwoStage` = time.2stage, `Joint` = time.joint))


# Create tables/figures to appropriately summarize the results 

write.csv(time, "D:/CU/Spring 2021/BIOS 7721/BIOS7721_FinalProject/Data/Time Summary.csv")
write.csv(sum.2stage, "D:/CU/Spring 2021/BIOS 7721/BIOS7721_FinalProject/Data/2stage Summary.csv")
write.csv(sum.joint, "D:/CU/Spring 2021/BIOS 7721/BIOS7721_FinalProject/Data/Joint Summary.csv")


  
                                                   
                                                   