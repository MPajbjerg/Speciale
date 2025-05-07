#Implementing an Euler-Maruyama scheme
eulerMaruyama <- function(stepSize,drift,diffusion,initialValue,numberOfSteps,time){
  #Initializing vector to store values
  output <- matrix(nrow=numberOfSteps+1, ncol= 2)
  output[1,1] <- time
  output[1,2] <- initialValue
  #Simulating brownian motion increments
  deltaW <- rnorm(numberOfSteps,0,sqrt(stepSize))
  for (i in 1:(numberOfSteps)){
    #calculating value
    output[i+1,2] <- output[i,2]+drift(output[i,2],output[i,1])*stepSize +
      diffusion(output[i,2],output[i,1])*deltaW[i]
    #increasing time for next step
    output[i+1,1] <- output[i,1]+stepSize
  }
  return(output)
}

#Implementing numerical integration by the trapezoidal rule
#Taken from LIV2 Week_5 file "ODEandInt"
numIntegrate <- function(lowerlim, upperlim, integrand, stepSize,middle){
  yval <- sapply(seq(lowerlim,upperlim,stepSize),integrand)
  middlevalues <- (yval[-1]+yval[-(length(yval))])*stepSize/2
  if (middle==FALSE){
  return(sum(middlevalues))}
  else return(middlevalues)
}

#Defining survival probability. Always calculated with stepsize 1/100 year
survivalProb <- function(t,s,mu,middle){
  exp(-numIntegrate(t,s,mu,0.01,middle))
}

rateAdjustedSurvivalProb <- function(t,s,mu,rate,middle){
  r_mu <- function(u){
    mu(u) + rate(u)
  }
  exp(-numIntegrate(t,s,r_mu,0.01,middle))
}
