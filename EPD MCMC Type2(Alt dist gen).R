# Clearing Old Variables from Environment
rm(list = ls(all = TRUE))

library(mcmc)

# Generate random nos from Exponential power distribution under PT-II censuring with binomial removals

Gen_R <- function(n, m, p) {
  
  r <- rep(0,m)
  
  for (i in 1:m - 1) {
    r[i] <- rbinom(1, n - m - sum(r), p)
  }
  
  r[m] <- n - m - sum(r)
  return (r)
}

Gen_EPD <- function(n, m, r, alpha, lambda) {  
  
  X <- rep(0,m)
  U <- rep(0,n)
  Z <- rep(0,n)
  W <- rep(0,n)
  
  U <- runif(n, min = 0, max = 1)
  Z <- (1 / lambda) * ((log(1 - log(1 - U))) ^ (1 / alpha))
  
  W <- sort(Z)
  
  for(i in 1:m) {
    X[i] = min(W)
    W = W[which(W!=min(W))]
    if(length(W) == 1){
      W = W
    }
    else {
      r2 = 0;
      for(k in 1:i) {
        r2 = r2 + r[k]
      }
      W = sample(W, n-i-r2, replace = F)
    }
  }
  return (X)
}

# Main function for estimating parameters

BayesAl <- function(N, n, m, p, alpha, lambda) {
  
  # Adjusting Starting values for MCMC to shorten converging times 
  alpha_ins <- alpha - 0.1
  lambda_ins <- lambda - 0.1
  burnIn <- 5000
  
  # We chose parameters of prior distribution large in such a way the mean of the prior distribution
  # is equal to the parameter of the EPD
  b1 <- alpha / 0.01
  a1 <- alpha * b1
  b2 <- lambda / 0.01
  a2 <- lambda * b2
  
  # Mean value 
  Mean_A <- rep(0, N)
  Mean_L <- rep(0, N)
  Mean_P <- rep(0, N)
  
  # Standard Errors
  SE_A <- rep(0, N)
  SE_L <- rep(0, N)
  SE_P <- rep(0, N)
  
  # Equal Density and the Coverage Probability
  EPD_A <- array(dim = c(N, 3))
  EPD_L <- array(dim = c(N, 3))
  EPD_P <- array(dim = c(N, 3))
  LED_A <- rep(0, N)
  RED_A <- rep(0, N)
  CPED_A <- rep(0, N)
  LED_L <- rep(0, N)
  RED_L <- rep(0, N)
  CPED_L <- rep(0, N)
  LED_P <- rep(0,N)
  RED_P <- rep(0,N)
  CPED_P <- rep(0,N)
    
  # Highest Posterior Density and the Coverage Probability
  HPD_A <- array(dim = c(N, 3))
  HPD_L <- array(dim = c(N, 3))
  HPD_P <- array(dim = c(N, 3))
  LHD_A <- rep(0, N)
  RHD_A <- rep(0, N)
  CPHD_A <- rep(0, N)
  LHD_L <- rep(0, N)
  RHD_L <- rep(0, N)
  CPHD_L <- rep(0, N)
  LHD_P <- rep(0,N)
  RHD_P <- rep(0,N)
  CPHD_P <- rep(0,N)
  
  for(i in 1:N) {
    
    # Generate EPD values with Binomial Removals
    R <- Gen_R(n = n, m = m, p = p)
    X <- Gen_EPD(n = n, m = m, r = R, alpha, lambda)
    
    # Function for mcmc Package
    lupost <- function(param) {
      
      stopifnot(is.numeric(param))
      stopifnot(is.finite(param))
      stopifnot(length(param) == 2)
      
      # Absolute values to prevent Nan as log of negative numbers not possible    
      a <- abs(param[1])
      l <- abs(param[2])
      
      # Log-Likelihood 
      logL <- sum ( log(a) +
                      (a * log(l)) +
                      ((a - 1) * log(X)) +
                      ((l * X) ^ a) +
                      (1) - 
                      (exp((l * X) ^ a)) )
      
      # Prior Distribution - Gamma
      aprior <- (a1 * log(b1)) +
        (a * (-b1)) +
        ((a1 - 1) * log(a)) -
        lgamma(a1)
      lprior <- (a2 * log(b2)) +
        (l * (-b2)) +
        ((a2 - 1) * log(l)) -
        (lgamma(a2))
      # Posterior Distribution
      return (logL + aprior + lprior)
    }
    
    # Metropolis Algorithm 
    out = metrop(lupost, initial = c(alpha_ins, lambda_ins), nbatch = 5000, blen = 1)
    
    #set scale so that acceptance rate near about 20%
    out <- metrop(out, scale = 0.21)
    out$accept
    out <- metrop(out, scale = 0.095)
    out$accept
    
    out <- metrop(out, nbatch = 20000, blen = 1)
    
    # Removing BurnIn Observations
    phi <- abs(out$batch[-(1:burnIn), ])
    phi_sort <- apply(phi, 2, sort)
    
    # Mean Estimates
    Mean_A[i] <- mean(phi[, 1])
    Mean_L[i] <- mean(phi[, 2])
    
    # Standard Error Estimates
    SE <- apply(phi, 2, sd)
    SE_A[i] <- SE[1]
    SE_L[i] <- SE[2]
    
    # Equal-Tailed Credible Interval
    ETCI <- function(phi_sort) {
      CI <- (quantile(phi_sort, probs = c(0.025, 0.975)))
      CI_width = CI[2] - CI[1]
      CPET <- (CI[1] < alpha & CI[2] > alpha)
      return (c(CI[1], CI[2], CPET))
    }
    
    EPD_A[i,] <- ETCI(phi_sort[,1])
    EPD_L[i,] <- ETCI(phi_sort[,2])
    
    # Highest Posterior Credible Interval
    HPDCI <- function(phi_sort , credMass = 0.95) {
      ciIdxInc <- floor(credMass * length(phi_sort))
      nCIs <- length(phi_sort) - ciIdxInc
      ci_width <- rep(0 , nCIs)
      for (i in 1:nCIs) {
        ci_width[i] <- phi_sort[i + ciIdxInc] - phi_sort[i]
      }
      HPDmin <- phi_sort[which.min(ci_width)]
      HPDmax <- phi_sort[which.min(ci_width) + ciIdxInc]
      CPHPD <- (HPDmin < alpha & HPDmax > alpha)
      return (c(HPDmin, HPDmax, CPHPD))
    }
    
    HPD_A[i,] <- HPDCI(phi_sort[,1])
    HPD_L[i,] <- HPDCI(phi_sort[,2])
    
    # Estimating 'p' using standard generation method (follows beta first kind distribution)
    # Choosing Jeffrey's prior as non-informative prior
    a <- 1/2
    b <- 1/2
    
    P <- rbeta(10000, a + sum(R[1:(m - 1)]), b + (m - 1) * (n - m) -  sum((R[1:(m - 1)]) * (m - (1:(m - 1)))))
    
    Mean_P[i] <- mean(P, na.rm = T)
    SE_P[i] <- sd(P)
    
    EP <- sort(P)
    EPD_P <- ETCI(P)
    HPD_P <- HPDCI(P)
    
  }
  
  Mean_A1 <- mean(Mean_A, na.rm = T)
  Mean_L1 <- mean(Mean_L, na.rm = T)
  Mean_P1 <- mean(Mean_P, na.rm = T)
  
  
  Bias_A1 <- mean(abs(Mean_A[i] - alpha))
  Bias_L1 <- mean(abs(Mean_L[i] - lambda))
  Bias_P1 <- mean(abs(Mean_P[i] - p))
  
  MSE_A1 <- mean((Mean_A[i] - alpha) ^ 2)
  MSE_L1 <- mean((Mean_L[i] - lambda) ^ 2)
  MSE_P1 <- mean((Mean_P[i] - p) ^ 2)
  
  SE_A1 <- mean(SE_A, na.rm = T)
  SE_L1 <- mean(SE_L, na.rm = T)
  SE_P1 <- mean(SE_P, na.rm = T)
  
  LED_A1 <- mean(EPD_A[, 1], na.rm = T)
  RED_A1 <- mean(EPD_A[, 2], na.rm = T)
  WED_A1 <- RED_A1 - LED_A1
  CPED_A1 <- mean(EPD_A[, 3], na.rm = T)
  LED_L1 <- mean(EPD_L[, 1], na.rm = T)
  RED_L1 <- mean(EPD_L[, 2], na.rm = T)
  WED_L1 <- RED_L1 - LED_L1
  CPED_L1 <- mean(EPD_L[, 3], na.rm = T)
  LED_P1 <- mean(LED_P, na.rm = T)
  RED_P1 <- mean(RED_P, na.rm = T)
  CPED_P1 <- mean(CPED_P, na.rm = T)
  WED_P1 <- RED_P1 - LED_P1
  
  LHD_A1 <- mean(HPD_A[, 1], na.rm = T)
  RHD_A1 <- mean(HPD_A[, 2], na.rm = T)
  WHD_A1 <- RHD_A1 - LHD_A1
  CPHD_A1 <- mean(HPD_A[, 3], na.rm = T)
  LHD_L1 <- mean(HPD_L[, 1], na.rm = T)
  RHD_L1 <- mean(HPD_L[, 2], na.rm = T)
  WHD_L1 <- RHD_L1 - LHD_L1
  CPHD_L1 <- mean(HPD_L[, 3], na.rm = T)
  LHD_P1 <- mean(LHD_P, na.rm = T)
  RHD_P1 <- mean(RHD_P, na.rm = T)
  CPHD_P1 <- mean(CPHD_P, na.rm = T)
  WHD_P1 <- RHD_P1 - LHD_P1
  
  A <- data.frame(n, m, alpha, Mean_A1, Bias_A1, MSE_A1, SE_A1, LED_A1, RED_A1, WED_A1, CPED_A1,
                  LHD_A1, RHD_A1, WHD_A1, CPHD_A1,
                  lambda, Mean_L1, Bias_L1, MSE_L1, SE_L1, LED_L1, RED_L1, WED_L1, CPED_L1, 
                  LHD_L1, RHD_L1, WHD_L1, CPHD_L1,
                  p, Mean_P1, Bias_P1, MSE_P1, SE_P1, LED_P1, RED_P1, WED_P1, CPED_P1, 
                  LHD_P1, RHD_P1, WHD_P1, CPHD_P1)
  return (A)
  
}

R1 <- BayesAl(N = 20, n = 50, m = 10, p = 0.1, alpha = 0.5, lambda = 0.2)
R2 <- BayesAl(N = 20, n = 50, m = 10, p = 0.3, alpha = 0.5, lambda = 0.2 )
R3 <- BayesAl(N = 20, n = 28, m = 35, p = 0.3, alpha = 0.5, lambda = 0.2 )
R3 <- BayesAl(N = 20, n = 42, m = 35, p = 0.3, alpha = 0.5, lambda = 0.2 )

DD1 <- rbind(R1, R2, R3, R4)
   
write.csv(x = DD1, file = "Exponential Power Distribution MCMC Estimates.csv") 