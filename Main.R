#Packages
library(readxl)
#packages for clustering cores
library(future.apply)
library(foreach)
library(doParallel)
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
  
  passiveFunc <- approxfun(seq(0.01,100,0.01),passive,rule = 2)
  
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
  
  b_a <- function(t, x) {
    value <- ifelse(x >= 10^(-25),
                    ifelse(t <= 67 - startAge, -72000 * 1.02^t,
                           ifelse(passiveFunc(t)>10^(-10),x / passiveFunc(t),0)),
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
#Defining the customer account and calculating and plotting.Note that this
#is the same in all scenario.
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
  
  plot(time + startAge,X_a[,2],type = "l",
       xlab = "Age in years", ylab = expression(X[a]), col = "black",lwd=2)
  
  X_a_Func <- approxfun(X_a,rule = 2)
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
    
    dX_e <- returnInvestment(t)*p_0_t*X_a_t - rho(t,X_a_t)*p_0_t*
        (muStarfunc(t) - trueMu(t)) - p_0_t*b_a(t,X_a_t) -
        b_ad(t,X_a_t)*p_0_t*trueMu(t)
    
    return(dX_e)
  }
  
  #Filling the portfolio wide mean based on low mortality rate.
  {
    #Fill in the true mortality rate in this simulation
  trueMu <- muLowfunc

  trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)

#Calculating the portfolio-wide mean at different times.
  portMeans[,2] <- rungeKutta(0.01,dynamicsX_E,X_a_initial,10000,0)

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
    portMeans[,4] <- rungeKutta(0.01,dynamicsX_E,X_a_initial,10000,0)
    
  }
  
  #Plotting all of the portfolio-wide means together.
  {
    plot(time + startAge,portMeans[,3],type = "l",
         xlab = "Age in years", ylab = expression(X[E]), col = "black",lwd=2)
    grid()
    
    lines(time+startAge,portMeans[,4],type = "l", col = "red", lwd = 2, lty = 2)
    
    lines(time + startAge, portMeans[,2],type = "l", col = "blue", lwd = 2, lty = 3)
    
    legend("bottomright",
           legend = c(expression(paste(mu,"*")), expression(mu[H]), expression(mu[L])),
           col = c("black", "red", "blue"),
           lty = c(1, 2, 3),
           lwd = 2,
           bty = "n")
    
  }
}

#In this section, we calculate the value of the guarantee.
#We do it by defining the dynamics of the guarantee and then
#using a runge-kutta scheme to calculate the accumulated value of the guarantee.
#We will not use the same runge-Kutta solver, as we will not stop, when the
#value becomes zero
{
  #We make a matrix, where w can store the values
  valueGuarantee <- matrix(NA, nrow= length(time), ncol = 4)
  
  valueGuarantee[,1] <- time
  
  colnames(valueGuarantee) <- c("Time","MuLow","MuStar","MuHigh")

  #Defining the dynamics of the guarantee, this is the dynamics of the
  #Profit process, when alpha=beta and Tilde(X)=X. Thus, we call the process P.
  dynamicsP_E <- function(t,x){
    p_0_t <- trueSurvivalProb(t)
    X_a_t <- X_a_Func(t)
    return(rho(t,X_a_t)*p_0_t*
             (muStarfunc(t) - trueMu(t)))}
  
  # We start off by using tho low mu
  {
    trueMu <- muLowfunc
  
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
  
    valueGuarantee[,2] <- rungeKuttaProfit(0.01,dynamicsP_E,0,10000,0)
    
  }
  
  #We do the same for the contractual mortality, where we expect zero profits.
  
  {
    trueMu <- muStarfunc
    
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
    
    valueGuarantee[,3] <- rungeKuttaProfit(0.01,dynamicsP_E,0,10000,0)
  }
  
  #And for the highMu
  {
    trueMu <- muHighfunc
    
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
    
    valueGuarantee[,4] <- rungeKuttaProfit(0.01,dynamicsP_E,0,10000,0)
    
  }
  
  #Plotting them together
  {
    plot(time + startAge,valueGuarantee[,3],type = "l",
         xlab = "Age in years", ylab = expression(P[E]), col = "black",lwd=2,
         ylim = c(-110000,110000))
    grid()
    
    lines(time+startAge,valueGuarantee[,4],type = "l", col = "red", lwd = 2, lty = 2)
    
    lines(time + startAge, valueGuarantee[,2],type = "l", col = "blue", lwd = 2, lty = 3)
    
    legend("topleft",
           legend = c(expression(paste(mu,"*")), expression(mu[H]), expression(mu[L])),
           col = c("black", "red", "blue"),
           lty = c(1, 2, 3),
           lwd = 2,
           bty = "n")
    
  }
  
}

#Now we make a setup, where we can calculate the reserve.
{
  #Before calculating the reserve, we need to calculate our X-modified
  #transition probabilities. We will reuse the above setup
  {
    a <- function(t){
      ifelse(t>67-startAge,ifelse(passiveFunc(t)>=10^(-10),1/passiveFunc(t),0),0)
    }
    
    b <- function(t){
      ifelse(t<= 67-startAge,-72000*1.02^t,0)
    }
    
    c_func <- function(t){
      ifelse(t<=67-startAge,1,0)
    }
    
    d <- function(t){
      0
    }
    
    #Now we are ready to calculate the modified transition probabilities.
    {
      #We define the differential equation
      dp_aa_X <- function(s,p){
        return(p*(interestRate(s)-a(s)+(c_func(s)-1)*trueMu(s)) -
          trueSurvivalProb(s)*(b(s)+d(s)*trueMu(s)))
      }
      
      #We then calculate all the modified prob. using a Runge-Kutta scheme.
      
      p_aa_X <- rungeKutta(0.01,dp_aa_X,X_a_initial,10000,0)
      
      p_aa_X_Func <- approxfun(time,p_aa_X,rule = 2)
      
    }
    
    #We can then calculate the reserve - we will do it as a
    #vector calculation to optimize time use.
    {
      V_mu_integrand <- rateAdjustedSurvivalProb(0,100,zeror,interestRate,TRUE)*
        (trueSurvivalProb(time)*(b(time)+d(time)*trueMu(time))+p_aa_X*(a(time)+c_func(time)*trueMu(time)))
      
      #Numerical integration
      
      V_mu <- sum((V_mu_integrand[-1]+V_mu_integrand[-(length(V_mu_integrand))])*0.01/2)
    }
    
  }
  
  
  #We will now replicate the setup in a pure unit-link setup to see,
  #if the reserve simplifies to the value of the customer account.
  {
    trueMuPure <- muStarfunc
    
    trueSurvivalProbPure <- approxfun(time,cumprod(survivalProb(0,100,trueMuPure,TRUE)),rule = 2)
    
    a_pure <- function(t){
      ifelse(t>67-startAge,ifelse(passiveFunc(t)>=10^(-10),1/passiveFunc(t),0),0)
    }
    
    b_pure <- function(t){
      ifelse(t<= 67-startAge,-72000*1.02^t,0)
    }
    
    c_pure <- function(t){
      1
    }
    
    d_pure <- function(t){
      0
    }
    
    
    #Now we are ready to calculate the modified transition probabilities.
    {
      #We define the differential equation
      dp_aa_X_pure <- function(s,p){
        return(p*(interestRate(s)-a(s)+(c_pure(s)-1)*trueMuPure(s)) -
                 trueSurvivalProbPure(s)*(b_pure(s)+d_pure(s)*trueMuPure(s)))
      }
      
      #We then calculate all the modified prob. using a Runge-Kutta scheme.
      
      p_aa_X_pure <- rungeKutta(0.01,dp_aa_X_pure,X_a_initial,10000,0)
      
      p_aa_X_Func_pure <- approxfun(time,p_aa_X_pure,rule = 2)
      
    }
    
    #We can then calculate the reserve - we will do it as a
    #vector calculation to optimize time use.
    {
      V_pure_integrand <- rateAdjustedSurvivalProb(0,100,zeror,interestRate,TRUE)*
        (trueSurvivalProbPure(time)*(b_pure(time)+d_pure(time)*trueMuPure(time))+p_aa_X_pure*
           (a_pure(time)+c_pure(time)*trueMuPure(time)))
      
      #Numerical integration
      
      V_pure <- sum((V_pure_integrand[-1]+V_pure_integrand[-(length(V_pure_integrand))])*0.01/2)
    }
    
    
  }
  
  
}



#We are now ready to calcualte the profit from the guarantee, the expected lifetimes
#and the reserves for all of the realizations of the stochastic mortality

{
  #First off, we make a matrix, where we can store all the values. We let each 
  #coloumn be a realization of mortality and each row correspond to something,
  #we want to calculate.
  
  dataAllSimulations <- matrix(,nrow = 3, ncol = dim(simulatedIntensity)[2])
  
  rownames(dataAllSimulations) <- c("Expected lifetime","Profit guarantee","Reserve")
  
  #We need to loop over all the different simulations. Start looping
  
  for (i in 1:100){
       #dim(simulatedIntensity)[2]){
    
    #Define sizes, which will be used for all mortality rates
    
    trueMu <- approxfun(time,simulatedIntensity[,i],rule = 2)
    
    trueSurvivalProb <- approxfun(time,cumprod(survivalProb(0,100,trueMu,TRUE)),rule = 2)
    
    #calculating the first output, which is the expected survival time.
    {
    dataAllSimulations[1,i] <- kappa(0,trueMu,zeror,TRUE)[1]
    }
    
    #Then we calculate the price of the guarantee. Note that the dynamics of
    #X_a and P are allready pre defined. Note that X_a are the same for all
    #mortality rates. This is the second output. Dynamics P_E is based on
    #the true mortality rate and survival probability defined earlier.
    
    {
      dataAllSimulations[2,i] <- rungeKuttaProfit(0.01,dynamicsP_E,0,
                                                  10000,0)[length(rungeKuttaProfit(0.01,dynamicsP_E,0,
                                                                                   10000,0))]
    }
    
    #For the final output the reserve is calculated.
    {
      #Firstly the modified transition probabilities are calculated
      {
        p_aa_X <- rungeKutta(0.01,dp_aa_X,X_a_initial,10000,0)
      
        p_aa_X_Func <- approxfun(time,p_aa_X,rule = 2)
      }
    
    #We can then calculate the reserve - we will do it as a
    #vector calculation to optimize time use.
      {
        V_mu_integrand <- rateAdjustedSurvivalProb(0,100,zeror,interestRate,TRUE)*
          (trueSurvivalProb(time)*(b(time)+d(time)*trueMu(time))+p_aa_X*(a(time)+c_func(time)*trueMu(time)))
      
        #Numerical integration
      
        dataAllSimulations[3,i] <- sum((V_mu_integrand[-1]+V_mu_integrand[-(length(V_mu_integrand))])*0.01/2)
      }
    }
    
    
    #End of loop
  }
  
  
  
  
  
}
