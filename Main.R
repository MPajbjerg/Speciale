#Sourcing functions
source("Functions.R")

#testing with gompertz-makeham from Dahl & Møller 2006.

{
  startAge <- 10
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
  
  
  #Testing how stochastic mortality develops.
  testStocMu <- eulerMaruyama(0.1,drift,diffusion,initial,1000,0)
  
  plot(testStocMu[,1]+startAge,testStocMu[,2],type = "l",
       xlab = "age in years", ylab = "intensity", col = "black", lwd = 2)
  grid()
  
  testStocMu2 <- eulerMaruyama(0.1,drift,diffusion,initial,1000,0)
  
  lines(testStocMu2[,1]+startAge,testStocMu2[,2],type = "l", col = "red",
        lwd = 2)
  
  testStocMu3 <- eulerMaruyama(0.1,drift,diffusion,initial,1000,0)
  
  lines(testStocMu3[,1]+startAge,testStocMu3[,2],type = "l", col = "blue",
        lwd = 2)
}
