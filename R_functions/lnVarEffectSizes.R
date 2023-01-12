
#### Created by A M Senior @ The University of Sydney, Australia, 18/05/2014

### Below are functions for caluculating log variance estimates of variance and it's associated measurement error.

### The first function requires the estimate of standard deviation and a sample size, the second just a sample size


Calc.lnSD<-function(SD, N){
  
  lnSD <- log(SD) + (1 / (2 * (N - 1)))
  
  return(lnSD)
}

# for variance
Calc.lnVAR<-function(VAR, N){
  
  lnVAR <- log(VAR) + (1 / (N - 1))
  
  return(lnVAR)
}

Calc.lVR<-function(CSD, CN, PSD, PN){
  
  lnVR <- log(PSD) - log(CSD) + 1 / (2 * (PN - 1)) - 1 / (2 * (CN - 1))
  
  return(lnVR)
}

# for variance
Calc.lVARR<-function(CVAR, CN, PVAR, PN){
  
  lVARR <- log(PVAR) - log(CVAR) + 1 / (PN - 1) - 1 / (CN - 1)
  
  return(lVARR)
}


Calc.var.lnSD<-function(N){
  
  var.lnSD <- (1 / (2 * (N - 1)))
  
  return(var.lnSD)
  
}

# for variance
Calc.var.lnVAR<-function(N){
  
  var.lnVAR <- (1 / (N - 1))
  
  return(var.lnVAR)
  
}

Calc.var.lVR<-function(CN, PN){
  
  var.lnVR <- 1 / (2 * (PN - 1)) + 1 / (2 * (CN - 1))
  
  return(var.lnVR)
}

# for variance
Calc.var.lVARR<-function(CN, PN){
  
  var.lVARR <- 1 / (PN - 1) + 1 /  (CN - 1)
  
  return(var.lVARR)
}

relinear <- function(ES){
  linear <- exp(ES) - 1
  return(linear)
}