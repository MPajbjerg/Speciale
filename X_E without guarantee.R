#Packages
library(readxl)
#packages for clustering cores and parallel calculations
library(future.apply)
library(foreach)
library(doParallel)
library(doFuture)
#Sourcing functions
source("Functions.R")

#Importing return and interest rate
{
  returnInterestRateData <- read_excel("Economic scenario speciale.xlsx")
  
  #Linear interpolation for interest rate  
  interestRate <- approxfun(returnInterestRateData$Time,
                            returnInterestRateData$Interest_rate,rule = 2)
  
  #Fitting a spline to the zero coupon bond price, which allows taking the
  #derivative to find the short rate.
  {
    zeroCouponPrices <- matrix(,nrow = length(returnInterestRateData$Time),
                               ncol = 2)
    zeroCouponPrices[,1] <- returnInterestRateData$Time
    
    for (i in 1:length(returnInterestRateData$Time)){
      zeroCouponPrices[i,2] <- 1/(1+returnInterestRateData$Interest_rate[i])^zeroCouponPrices[i,1]
    }
    
    shortRate <- splinefun(zeroCouponPrices[,1],-log(zeroCouponPrices[,2]),
                           method = "natural")
    
    shortRateFunc <- function(t){
      shortRate(t,deriv = 1)
    }
  }
  #Linear interpolation for return on investements
  returnInvestment <- shortRateFunc
  
  #Plotting the economic scenario
  {
    plot(seq(0,100,0.1),sapply(seq(0,100,0.1),shortRateFunc),type = "l",
         xlab = "Time in years", ylab = "Short rate", col = "black", lwd = 2,
         ylim = c(0.02,0.055))
    grid()
    
    #lines(seq(0,100,0.1),sapply(seq(0,100,0.1),returnInvestment),
    #      type = "l", col = "red",lwd = 2)
    
    #legend("topright",legend = c("Interest rate", "Return on assets"),
    #       col = c("black", "red"),lwd = 2,lty = 1, bty = "n",
    #       cex = 0.55)
  }
}

#Test mortality for sanity check of passive calculator
{
  testMuData <- read_excel("Intensitet til test.xlsx")
  
  #Interpolating
  
  testMuFunc <- approxfun(testMuData$Tid,testMuData$Intensitet, rule = 2)
}

#testing with gompertz-makeham from Dahl & Møller 2006.
#Section with simulation of mortality 
{
  startAge <- 47
  alpha <- 0.000134
  #Beta from Dahl & Møller 0.0000353
  beta <- 0.0000353
  c <- 1.1020
  
  gompertzMu <- function(t){
    alpha + beta*c^(startAge+t-10)
  }
  
  gompertzMuDerivative <- function(t){
    beta*log(c)*c^(startAge+t-10)
  }
  
  #Defining parameters
  deltaTilde <- 0.2 
  gammaTilde <- 0.008
  sigmaTilde <- 0.02
  
  delta <- deltaTilde
  gamma <- function(t){
    deltaTilde*exp(-gammaTilde*t)
  }
  
  sigma <- sigmaTilde
  
  deltaMu <- function(t){
    delta - (gompertzMuDerivative(t)/gompertzMu(t)) 
  }
  
  gammaMu <- function(t){
    gamma(t)*gompertzMu(t)
  }
  
  sigmaMu <- function(t){
    sigma*sqrt(gompertzMu(t))
  }
  
  drift <- function(mu,t){
    gammaMu(t)-deltaMu(t)*mu
  }
  
  diffusion <- function(mu,t){
    sigmaMu(t)*sqrt(mu)
  }
  
  
  initial <- gompertzMu(0)
  
  #Defining contract mu for the longevity guarantee.
  #Just the other mortality where sigma=0
  #We need a function, which always outputs 0
  zero <- function(mu,t){
    0
  }
  
  muStar <- eulerMaruyama(0.01,drift,zero,initial,10000,0)
  
  #Initializing time to make a matrix of simulations
  time <- seq(0,100,0.01)
  
  #Making simulations using every core in processor.
  #Only works for Windows
  {
    set.seed(47)
    cat("Setting up multicore use", format(Sys.time(),"%H:%M:%S"),"\n")
    
    plan(multisession)
    
    cat("Multicore ready", format(Sys.time(),"%H:%M:%S"),"\n")
    
    cat("Starting simulation at", format(Sys.time(),"%H:%M:%S"),"\n")
    
    simulatedIntensity <- future_replicate(2000,{eulerMaruyama(0.01,drift,diffusion,initial,10000,0)[,2]},
                                           simplify = "matrix")
    
    cat("Simulation ended at", format(Sys.time(),"%H:%M:%S"),"\n")
    
    cat("Ending multicore", format(Sys.time(),"%H:%M:%S"),"\n")
    
    future::plan(sequential)
    
    cat("Multicore ended", format(Sys.time(),"%H:%M:%S"),"\n")
    
    #clearing meamory of unused stuff
    gc()
  }
  
  #Making a matrix consisting of all the differences to 
  #our contractual intensity
  logDifferenceMatrix <- log(simulatedIntensity) -log(muStar[,2])
  
  #Only keeping points, where mortality rate is larger than
  #contractual mortality. #It took too much memory to store all
  #"mellemregninger". THus the not so nice code below.
  
  highMuIndex <- which.max(colSums(ifelse(logDifferenceMatrix >0,logDifferenceMatrix,0)))
  
  #Defining the mortality with the highest log difference.
  
  highMu <- simulatedIntensity[,highMuIndex]
  
  #Also finding the mortality with the highest negative log difference
  
  lowMuIndex <- which.min(colSums(ifelse(logDifferenceMatrix <0,logDifferenceMatrix,0)))
  
  lowMu <- simulatedIntensity[,lowMuIndex]
  
  #Plotting the highest and the lowest mortality rate together with the base
  {
    plot(time + startAge,log(muStar[,2]),type = "l",
         xlab = "Age in years", ylab = expression(log(mu)), col = "black",lwd=2)
    grid()
    
    lines(time+startAge,log(highMu),type = "l", col = "red", lwd = 2, lty = 2)
    
    lines(time + startAge, log(lowMu),type = "l", col = "blue", lwd = 2, lty = 3)
    
    legend("bottomright",
           legend = c(expression(paste(mu,"*")), expression(mu[H]), expression(mu[L])),
           col = c("black", "red", "blue"),
           lty = c(1, 2, 3),
           lwd = 2,
           bty = "n")
    
  }
  
  #Plotting without taking logarithms
  {
    plot(time + startAge,muStar[,2],type = "l",
         xlab = "Age in years", ylab = "mu", col = "black",lwd=2)
    grid()
    
    lines(time+startAge,highMu,type = "l", col = "red", lwd = 2, lty = 2)
    
    lines(time + startAge, lowMu,type = "l", col = "blue", lwd = 2, lty = 3)
  }
  
  #Making mortality functions by linear interpolation
  muStarfunc <- approxfun(time, muStar[,2], rule = 2)
  muLowfunc <- approxfun(time,lowMu, rule = 2)
  muHighfunc <- approxfun(time, highMu, rule = 2)
}

#Calculating survival times and passive
{
  #Calculating passive one time for good, since the passive relies on the
  #contractual mortality, which is the same for all simulations. Plotting
  #all the expected remaining lifetimes.
  #creating an interest rate 0 function, since kappa needs function input.
  
  zeror <- function(t){
    0
  }
  
  #calculating for test intensity for sanity check
  
  lifetimeMuTest <- kappa(0,testMuFunc,zeror,TRUE)
  #Looks about right for mortalities from 2023.
  
  #calculating expected remaining lifetime
  
  lifetimeMuStar <- kappa(0,muStarfunc,zeror,TRUE)
  
  lifetimeMuLow <- kappa(0,muLowfunc,zeror,TRUE)
  
  lifetimeMuHigh <- kappa(0,muHighfunc,zeror,TRUE)
  
  #Calculating passive, as this is the same in all scenario.
  #The passive is calculated for future times as well
  
  #Making a wrapper for the Sapply
  
  wrapperPassive <- function(t){
    kappa(t,muStarfunc,shortRateFunc,FALSE)
  }
  
  #Takes some time, so uses more cores
  {
    plan(multisession)
    passive <- future_sapply(seq(0,100,0.01),wrapperPassive)
    plan(sequential)
    gc()
  }
  #Making the passive into a function, since some time steps are missing,
  #due to double trapez integration.
  
  
  passiveFunc <- approxfun(seq(0,100,0.01),passive,rule = 2)
  
  #Low mort scenario
  wrapperPassiveLow <- function(t){
    kappa(t,muLowfunc,shortRateFunc,FALSE)
  }
  
  #Takes some time, so uses more cores
  {
    plan(multisession)
    passiveLow <- future_sapply(seq(0,100,0.01),wrapperPassiveLow)
    plan(sequential)
    gc()
  }
  #Making the passive into a function, since some time steps are missing,
  #due to double trapez integration.
  
  
  passiveFuncLow <- approxfun(seq(0,100,0.01),passiveLow,rule = 2)
 
  
  #High mort scenario
  wrapperPassiveHigh <- function(t){
    kappa(t,muHighfunc,shortRateFunc,FALSE)
  }
  
  #Takes some time, so uses more cores
  {
    plan(multisession)
    passiveHigh <- future_sapply(seq(0,100,0.01),wrapperPassiveHigh)
    plan(sequential)
    gc()
  }
  #Making the passive into a function, since some time steps are missing,
  #due to double trapez integration.
  
  
  passiveFuncHigh <- approxfun(seq(0,100,0.01),passiveHigh,rule = 2) 
  
}


#Defining payment functions and sum at risk. 
#We run into 'dividing by zero' issues, when the passive becomes very small.
#Thus we will close the account, when the account is too small
{
  b_ad <- function(t,x){
    ifelse(x>= 10^(-25) & t<=67-startAge ,
           x,
           0)
  }
  
  #also added if the payment is larger than the account payout the account.
  b_a <- function(t, x) {
    value <- ifelse(x >= 10^(-25),
                    ifelse(t <= 67 - startAge, -72000 * 1.02^t,
                           ifelse(passiveFunc(t)>10^(-10),x / passiveFunc(t),0)),
                    0)
    return(value)
  }
  
  b_aLow <- function(t, x) {
    value <- ifelse(x >= 10^(-25),
                    ifelse(t <= 67 - startAge, -72000 * 1.02^t,
                           ifelse(passiveFuncLow(t)>10^(-10),x / passiveFuncLow(t),0)),
                    0)
    return(value)
  }
  
  b_aHigh <- function(t, x) {
    value <- ifelse(x >= 10^(-25),
                    ifelse(t <= 67 - startAge, -72000 * 1.02^t,
                           ifelse(passiveFuncHigh(t)>10^(-10),x / passiveFuncHigh(t),0)),
                    0)
    return(value)
  }
  
  
  rho <- function(t,x){
    ifelse (x>=10^(-25),
            b_ad(t,x)-x
            ,
            0)
  }
}
#
{
  #We run into problems, when numerically dividing by zero.
  #We thus define the CA as practically zero, when this happens
  dynamicsX_a <- function(t,x){
    ifelse (x>= 10^(-25),
            returnInvestment(t)*x-b_a(t,x)-rho(t,x)*muStarfunc(t)
            ,
            0)
  }
  
  X_a_initial <- 1000000
  
  X_a <- cbind(time,rungeKutta(0.01,dynamicsX_a,X_a_initial,10000,0))
  
  
  X_a_Func <- approxfun(X_a,rule = 2)
}
{
  #We run into problems, when numerically dividing by zero.
  #We thus define the CA as practically zero, when this happens
  dynamicsX_aLow <- function(t,x){
    ifelse (x>= 10^(-25),
            returnInvestment(t)*x-b_aLow(t,x)-rho(t,x)*muLowfunc(t)
            ,
            0)
  }
  
  X_a_initial <- 1000000
  
  X_aLow <- cbind(time,rungeKutta(0.01,dynamicsX_aLow,X_a_initial,10000,0))
  
  
  X_a_FuncLow <- approxfun(X_aLow,rule = 2)
}
{
  #We run into problems, when numerically dividing by zero.
  #We thus define the CA as practically zero, when this happens
  dynamicsX_aHigh <- function(t,x){
    ifelse (x>= 10^(-25),
            returnInvestment(t)*x-b_aHigh(t,x)-rho(t,x)*muHighfunc(t)
            ,
            0)
  }
  
  X_a_initial <- 1000000
  
  X_aHigh <- cbind(time,rungeKutta(0.01,dynamicsX_aHigh,X_a_initial,10000,0))
  
  X_a_FuncHigh <- approxfun(X_aHigh,rule = 2)
}


#We now calculate the portfolio-wide mean of a portfolio
#consisting of male individuals at age 47.
#We make this in a setup, where a true mortality is specified each time
#Then we can replicate this setup, when calculating from every realization of

{
  #Initializing matrix for portfolio-wide means.
  portMeans <- matrix(NA, nrow= length(time), ncol = 4)
  
  colnames(portMeans) <- c("Time","MuLow","MuStar","MuHigh")
  
  portMeans[,1] <- time
  
  #Defining dynamics of portfolio-wide mean, which only works, when
  #true Mu and true survivalprobability have been defined
  dynamicsX_E <- function(t,x){
    
    p_0_t <- trueSurvivalProb(t)
    X_a_t <- X_a_Func(t)
    
    dX_e <- returnInvestment(t)*p_0_t*X_a_t - p_0_t*b_a(t,X_a_t) -
      b_ad(t,X_a_t)*p_0_t*trueMu(t)
    
    return(dX_e)
  }
  
  dynamicsX_ELow <- function(t,x){
    
    p_0_t <- trueSurvivalProb(t)
    X_a_t <- X_a_FuncLow(t)
    
    dX_e <- returnInvestment(t)*p_0_t*X_a_t - p_0_t*b_aLow(t,X_a_t) -
      b_ad(t,X_a_t)*p_0_t*trueMu(t)
    
    return(dX_e)
  }
  
  dynamicsX_EHigh <- function(t,x){
    
    p_0_t <- trueSurvivalProb(t)
    X_a_t <- X_a_FuncHigh(t)
    
    dX_e <- returnInvestment(t)*p_0_t*X_a_t - p_0_t*b_aHigh(t,X_a_t) -
      b_ad(t,X_a_t)*p_0_t*trueMu(t)
    
    return(dX_e)
  }
  
  #Filling the portfolio wide mean based on low mortality rate.
  {
    #Fill in the true mortality rate in this simulation
    trueMu <- muLowfunc
    
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
    
    #Calculating the portfolio-wide mean at different times.
    portMeans[,2] <- rungeKutta(0.01,dynamicsX_ELow,X_a_initial,10000,0)
    
  }
  #Setting the true mu equal to the contractual mu
  {
    #Fill in the true mortality rate in this simulation
    trueMu <- muStarfunc
    
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
    
    #Calculating the portfolio-wide mean at different times.
    portMeans[,3] <- rungeKutta(0.01,dynamicsX_E,X_a_initial,10000,0)
    
  }
  
  #Setting the true mu equal to the high mortality rate.
  {
    #Fill in the true mortality rate in this simulation
    trueMu <- muHighfunc
    
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
    
    #Calculating the portfolio-wide mean at different times.
    portMeans[,4] <- rungeKutta(0.01,dynamicsX_EHigh,X_a_initial,10000,0)
    
  }
  
  #Plotting all of the portfolio-wide means together.
  {
    plot(time + startAge,portMeans[,3],type = "l",
         xlab = "Age in years", ylab = expression(X[E]* " Without guarantee"), col = "black",lwd=2)
    grid()
    
    lines(time+startAge,portMeans[,4],type = "l", col = "red", lwd = 2, lty = 2)
    
    lines(time + startAge, portMeans[,2],type = "l", col = "blue", lwd = 2, lty = 3)
    
    legend("topright",
           legend = c(expression(paste(mu,"*")), expression(mu[H]), expression(mu[L])),
           col = c("black", "red", "blue"),
           lty = c(1, 2, 3),
           lwd = 2,
           bty = "n")
    
  }
}

