




#' @title Mean incidence from the renewal equation.
#'
#' @param t Numerical. Time.
#' @param R0 Numerical. Basic reproduction number.
#' @param B Numerical vector. Multiplicative factor changing the contact rate.
#' (if no change, set \code{B[t]=1} for all t)
#' @param S Numerical for the number of susceptible.
#' @param N Numerical. Total population size
#' @param alpha Numerical. Parameter representing mixing heterogeneity.
#' \code{alpha = 0} means homogeneous mixing.
#' @param g Numerical vector representing the intrinsic generation interval distribution.
#' @param I Numerical vector. Incidence at each time step.
#'
#' @return
#'
#' @examples
#' 
mean_inc <- function(t, R0, B, S, N, alpha, g, I) {
  
  #stopifnot(t>1)
  
  # work with existing data
  n = sum(!is.na(I))
  #stopifnot(n>0)
  
  tmp1 = R0 * B[t] * (S[t-1]/N)^(exp(alpha))  
  
  revI = rev(I[1:n])
  n2   = min(n,length(g))
  tmp2 = g[1:n2] * revI[1:n2] 
  
  m = tmp1 * sum(tmp2)
  return(m)
}

calc_I <- function(lambda, deterministic){
  res = lambda
  if(!deterministic) res = rpois(n=1, lambda = lambda)
  return(res)
}

calc_Y <- function(lambdaY,n, deterministic){
  res = lambdaY
  if(!deterministic) res = rpois(n=n, lambda = lambdaY)
  return(res)
}


#' @title Simulate an epidemic based on the renewal equation.
#' 
#' @description
#' Simulate an epidemic base on a renewal equation that incorporates
#' both the clinical and wastewater data (ie fecal shedding is 
#' explicitly taken into account).
#'
#' @param prm List of model parameters.
#' @param deterministic Logical. Is the simulation deterministic. 
reem_simulate <- function(prm, deterministic) {
  
  # Unpack parameters
  R0      = prm$R0
  B       = prm$B
  N       = prm$N
  alpha   = prm$alpha
  I.init  = prm$I.init
  horizon = prm$horizon
  rho     = prm$rho
  lag     = prm$lag
  g       = prm$g
  fec     = prm$fec
  kappa   = prm$kappa
  psi     = prm$psi
  t.obs.ww = prm$t.obs.ww
  
  ni = length(I.init)
  
  m = numeric(horizon)
  I = rep(NA, horizon)
  S = rep(NA, horizon)
  A = rep(NA, horizon)
  Y = rep(NA, horizon)
  Wd = rep(NA, horizon)
  
  # Initial period when incidence is known:
  m[1:ni] = I.init
  I[1:ni] = I.init
  S[1:ni] = N - cumsum(I.init)
  
  for(t in 1:ni){
    tlag = max(1, t-lag)
    A[t] = sum(I[t:tlag])
  }
  
  if(!deterministic) rho = rnorm(n=horizon, mean = rho, sd = 0.05)
  rho[rho <= 0] <- 0.001
  
  for(t in (ni+1):horizon){ # t = ni+1
    
    m[t] = mean_inc(t, R0, B, S, N, alpha, g, I)
    I[t] = calc_I(lambda = m[t], deterministic)
    S[t] = max(0, S[t-1] - I[t])
    
    tlag = max(1, t-lag)
    A[t] = sum(I[t:tlag])
  }
  
  # Observed aggregated incidence (Y):
  lambdaY = rho*A
  lambdaY[lambdaY==0] <- 1e-3
  lambdaY[is.na(lambdaY)] <- 1e-3
  Y = calc_Y(lambdaY, n=length(A), deterministic)
  
  # --- Wastewater ---
  
  # -- deposited
  
  nf = length(fec)
  w.m = numeric(horizon)
  for(t in 1:horizon){
    idx = 1:min(t-1,nf)
    w.m[t] = sum(fec[idx]*I[t-idx])
  }
  Wd = w.m
  if(!deterministic) Wd = rnorm(n = horizon, mean = w.m, sd = w.m * 0.1)
  
  # -- present at sampling site
  
  n.psi= length(psi)
  wp.m = numeric(horizon)
  
  for(t in 1:horizon){
    idx = 1:min(t-1,n.psi)
    wp.m[t] = sum(psi[idx] * Wd[t-idx] * exp(-kappa * idx))
  }
  Wp = wp.m
  if(!deterministic) Wp = rnorm(n=horizon, mean = wp.m, sd = wp.m * 0.1)
  
  # -- observed ("reported") at sampling site
  
  t.obs.ww = t.obs.ww[t.obs.ww < horizon]
  n.ww = length(t.obs.ww)
  wr.m = numeric(n.ww)
  
  for(i in seq_along(t.obs.ww)) 
    wr.m[i] = Wp[t.obs.ww[i]]
  
  Wr = wr.m
  if(!deterministic) Wr = rnorm(n=n.ww, mean = wr.m, sd = wr.m * 0.2)
  
  tmp = data.frame(t = t.obs.ww, Wr = Wr)
  
  df = data.frame(
    t = 1:horizon, 
    m = m, 
    I = I,
    S = S,
    A = A,
    Y = Y,
    Wd = Wd, 
    Wp = Wp)
  
  df = left_join(df, tmp, by='t')
  
  return(df)    
}



