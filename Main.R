#Packages
library(readxl)
library(future.apply)
#Sourcing functions
source("Functions.R")

#Importing return and interest rate
{
returnInterestRateData <- read_excel("Economic scenario speciale.xlsx")

#Linear interpolation for interest rate  
interestRate <- approxfun(returnInterestRateData$Time,
                          returnInterestRateData$Interest_rate,rule = 2)

#Linear interpolation for return on investements
returnInvestment <- approxfun(returnInterestRateData$Time,
                              returnInterestRateData$Return,rule = 2)

#Plotting the economic scenario
{
plot(seq(0,100,0.1),sapply(seq(0,100,0.1),interestRate),type = "l",
     xlab = "Time in years", ylab = "Rate", col = "black", lwd = 2,
     ylim = c(0.02,0.055))
grid()

lines(seq(0,100,0.1),sapply(seq(0,100,0.1),returnInvestment),
      type = "l", col = "red",lwd = 2)

legend("topright",legend = c("Interest rate", "Return on assets"),
       col = c("black", "red"),lwd = 2,lty = 1, bty = "n",
       cex = 0.55)
}
}

#Test mortality for sanity check of passive calculator
{
  testMuData <- returnInterestRateData <- read_excel("Intensitet til test.xlsx")
  
  #Interpolating
  
  testMuFunc <- approxfun(testMuData$Tid,testMuData$Intensitet, rule = 2)
}

#testing with gompertz-makeham from Dahl & Møller 2006.
#Section with simulation of mortality 
{
  startAge <- 47
  alpha <- 0.000134
  beta <- 0.0000353
  c <- 1.1020
  
  gompertzMu <- function(t){
    alpha + beta*c^(startAge+t)
  }
  
  gompertzMuDerivative <- function(t){
    beta*log(c)*c^(startAge+t)
  }
  
  #Defining parameters
  deltaTilde <- 0.2 
  gammaTilde <- 0.008
  sigmaTilde <- 0.02

  delta <- deltaTilde
  gamma <- function(t){
    g <- deltaTilde*exp(-gammaTilde*t)
    return(g)
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
    
    simulatedIntensity <- future_replicate(10000,{eulerMaruyama(0.01,drift,diffusion,initial,10000,0)[,2]},
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
  
  passive <- kappa(0,muStarfunc,interestRate,TRUE)
  
  #Making the passive into a function, since som time steps are missing,
  #due to double trapez integration.
  
  passiveFunc <- approxfun(seq(0.01,99.99,0.01),passive,rule = 2)
  
}
#Defining payment functions
{
  b_ad <- function(t,x){
    if(t<=67){
      value <- x
    }
    else {value <- 0}
    return(value)
  }
  
  b_a <- function(t,x){
    if(t<=67){
      -72000*1.02^(t-startAge)
    }
    else {}
  }
  
  rho <- function(t,x){
    b_ad(t,x)-x
  }
}