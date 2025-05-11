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

rungeKutta <- function(stepSize,f,initialValue,numberOfSteps,time){
  y <- numeric(numberOfSteps + 1)
  y[1] <- initialValue
  t <- seq(time,numberOfSteps*stepSize,stepSize)
  for (i in 1:numberOfSteps) {
    k1 <- f(t[i], y[i])
    k2 <- f(t[i] + stepSize / 2, y[i] + stepSize / 2 * k1)
    #Stopping the algortihm close to zero.
    ifelse( y[i]>10^(-10),
    y[i + 1] <- y[i] + stepSize * k2,y[i+1] <- 0)
  }
  
  return(y)
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
  p_t_s <- exp(-numIntegrate(t,s,mu,0.01,middle))
  #The probability should be 1 at time of evaluation
  return(c(1,p_t_s))
}

#Note here that I have imported the zero coupon spot rate, which is NOT
#the forward rate. We work in time intervals of 1/100 pr. year.
#Thus we will have to make a rate with the same step size
#While taking into account compounding rates.
#Thus we have to take the 100th root
rateAdjustedSurvivalProb <- function(t,s,mu,rate,middle){
  r_mu <- function(u){
    mu(u) + (1+rate(u))^(1/100)-1
  }
  pr_t_s <- exp(-numIntegrate(t,s,r_mu,0.01,middle))
  return(c(1,pr_t_s))
}


kappa <- function(t,mu,rate,forAllt){
  #The adjustet survival probability function outputs the integral of
  #each increment. thus to get each point, which needs to be integrated
  #we need to multiply each step. Then we can use the trapezoidal
  #rule on this 
  
  middleValues1 <- rateAdjustedSurvivalProb(t,100,mu,rate,TRUE)
  
  middleValues2 <- cumprod(middleValues1)
  
  passive <- (middleValues2[-1]+middleValues2[-(length(middleValues2))])*0.01/2
  
  passiveForAllt <- sum(passive) - c(0,cumsum(passive[-length(passive)]))
  
  if (forAllt==FALSE){
    return(sum(passive))}
  else return(passiveForAllt)
}





