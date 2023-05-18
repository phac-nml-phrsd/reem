
#' Class to represent a Renewal Equation Epidemic Model.
#'
#' @field name character. 
#' @field prms list. 
#' @field is.fitted logical. 
#' 
#' @return
#' @export
#'
#' @examples
#' 
setRefClass(
  Class = 'reem', 
  fields =  list(
    name      = "character",
    prms      = "list",
    is.fitted = "logical"
  ),
  
  # - - - - - - - - - - - - - - - - - -
  # - - - `reem` CLASS METHODS 
  # - - - - - - - - - - - - - - - - - -
  
  methods = 
    list(
      
      print_prms = function(){
        cat('\n--- Parameters for REEM `',name,'`\n')
        for(i in seq_along(prms)){
          cat(names(prms)[i],' = ', prms[[i]],'\n')
        }
        cat(' --------------------------\n')
      },
      
      
      simulate = function(deterministic){
        
        # Unpack parameters
        R0      = prms$R0
        B       = prms$B
        N       = prms$N
        alpha   = prms$alpha
        I.init  = prms$I.init
        horizon = prms$horizon
        rho     = prms$rho
        lag     = prms$lag
        g       = prms$g
        fec     = prms$fec
        kappa   = prms$kappa
        psi     = prms$psi
        t.obs.ww = prms$t.obs.ww
        
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
        
        df = dplyr::left_join(df, tmp, by='t')
        
        return(df)  
      }
      
    ))




