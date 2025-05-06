#Packages
library(readxl)

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
#testing with gompertz-makeham from Dahl & Møller 2006.

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
  
  set.seed(1)
  #Testing how stochastic mortality develops.
  testStocMu <- eulerMaruyama(0.01,drift,diffusion,initial,10000,0)
  
  plot(testStocMu[,1]+startAge,testStocMu[,2],type = "l",
       xlab = "age in years", ylab = "intensity", col = "black", lwd = 2)
  grid()
  
  testStocMu2 <- eulerMaruyama(0.01,drift,diffusion,initial,10000,0)
  
  lines(testStocMu2[,1]+startAge,testStocMu2[,2],type = "l", col = "red",
        lwd = 2)
  
  testStocMu3 <- eulerMaruyama(0.01,drift,diffusion,initial,10000,0)
  
  lines(testStocMu3[,1]+startAge,testStocMu3[,2],type = "l", col = "blue",
        lwd = 2)
  
  #Testing, if we can calculate survival probabilities with interpolated functions
  
  testStocMufunc <- approxfun(testStocMu,rule = 2)
  
  survivalProb(0,1,testStocMufunc)
}

#Defining payment functions and passive
{
  kappa <- function(t,mu){
    #creating wrapper function, which can be integrated
    exponential <- function(u){survivalProb(t,u,mu) + interestRate(u)}
     max(numIntegrate(t,300,exponential,0.01),1)
  }
  
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