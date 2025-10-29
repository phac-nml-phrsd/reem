#include <Rcpp.h>
#include <numeric>
#include <cmath>



using namespace Rcpp;


double mean_inc(int t, double R0, 
                NumericVector B, 
                NumericVector S, 
                int N, double alpha, 
                NumericVector g, 
                NumericVector I){
  
  // STOPPED HERE
  // Check the idex limits, probably off by 1 or something like that... 
  // Compare with R code... 
  
  int n = t;
  
  double tmp1 = R0 * B[t] * std::pow(S[t-1]/N, exp(alpha)) ;
  
  std::vector<double> tmp2(n, 0.0);
  
  int kmax = n-1; 
  if(g.size()-1 < n-1) kmax = g.size()-1;
  
  for(int k = 0; k < kmax; k++){
    tmp2[k] = g[k] * I[n-k];
  }
  double sumtmp2 = std::accumulate(tmp2.begin(), tmp2.end(), 0.0);
  
  double res = tmp1 * sumtmp2;
  //   
  //   revI = I[n:1]  # faster than using `rev()`
  // n2   = min(n,length(g))
  //   tmp2 = g[1:n2] * revI[1:n2] 
  // 
  // m = tmp1 * sum(tmp2)
   
  return res;
}


// [[Rcpp::export]]
DataFrame simul_C(List prms) {

  // unpack all parameters 
  int horizon = prms["horizon"];
  double R0 = prms["R0"];
  
  //NumericVector B = prms["B"]["mult"];
  NumericVector B(horizon, 1.0);
  
  int N = prms["N"];
  double alpha = prms["alpha"];
  NumericVector I_init = prms["I.init"];
  double rho = prms["rho"];
  int lag = prms["lag"];
  NumericVector g = prms["g"];
  NumericVector fec = prms["fec"];
  double h_prop = prms["h.prop"];
  NumericVector h_lags = prms["h.lags"];
  double kappa = prms["kappa"];
  NumericVector psi = prms["psi"];
  double shed_mult = prms["shed.mult"];
  
  int ni = I_init.size();
  
  // Rcout << "DEBUG 1" << std::endl;
  
  
  NumericVector cum_I_init(ni);
  double tmp = 0.0;
  for (int i = 0; i < ni; ++i) {
    tmp += I_init[i];
    cum_I_init[i] = tmp;
  }
  
  
  // // create epi vectors
  NumericVector m(horizon);
  NumericVector I(horizon);
  NumericVector S(horizon);
  NumericVector A(horizon);

 // Initial period when incidence is known:
  for(int i = 0; i < ni; i++) {
    m[i] = I_init[i];
    I[i] = I_init[i];
    S[i] = N - cum_I_init[i];
  }
  A[0] = I[0]; 
  for(int t = 1; t < ni; t++){
    int tlag = std::max(0, t-lag);
    // TODO: try to get std::accumulate to work...
    double s = 0;
    for(int k = tlag; k <= t; k++){
      s += I[k];
    }
    A[t] = s;
  }
  
  for(int t = ni; t < horizon; t++){
    m[t] = mean_inc(t, R0, B, S, N, alpha, g, I);
    I[t] = m[t];
    S[t] = S[t-1] - I[t];
    if(S[t]<0) S[t] = 0;
    
    int tlag = std::max(0, t-lag);
    double s = 0;
    for(int k = tlag; k <= t; k++){
      s += I[k];
    }
    A[t] = s;
  }
  // Observed aggregated incidence (Y): 
  std::vector<double> Y(horizon, 0.0);
  
  for(int i=0; i<horizon; i++ ){
    tmp = rho * A[i];
    if(tmp < 1e-9) tmp = 1e-3;
    Y[i] = tmp;
  }
   
  // Hospital admissions
  int nh = h_lags.size();
  double sumh = std::accumulate(h_lags.begin(), h_lags.end(),0.0);
  std::vector<double> hnorm(nh);
  for(int i = 0; i< nh; i++){
    hnorm[i] = h_lags[i] / sumh * h_prop;
  } 
  
  std::vector<double> H(horizon, 0);
  for(int t =2; t < horizon; t++){
    H[t] = 0;
    int upperidx = std::min(nh, t-1);
    for(int k=0; k< upperidx; k++) {
      H[t] +=  hnorm[k] * I[t-k];
    }
    H[t] = std::round(H[t]);
  } 
  
  std::vector<double> Hpercapita(horizon, 0);
  for(int i=0; i < horizon; i++) 
    Hpercapita[i] = H[i] / N;
    
    
  // Wastewater
  
  int nf = fec.size();
  int imax;
  std::vector<double> Wd(horizon, 0);
  for(int t = 0; t < horizon; t++){
    imax = std::min(t-1, nf);
    std::vector<double> z(nf,0);
    for(int i=0; i < imax; i++) 
      z[i] = fec[i] * I[t - i];
    double s = 0;
    for(int i = 0; i < z.size(); i++) 
      s+= z[i];
    Wd[t] = shed_mult * s;
  }
  
  int npsi = psi.size();
  std::vector<double> Wp(horizon, 0);
  for(int t = 1; t<horizon; t++){
    imax = std::min(t-1, npsi);
    std::vector<double> y(imax, 0);
    for(int i=0; i<imax; i++){
      y[i] = psi[i] * Wd[t-i] * exp(-kappa * i);
    }
    double s = 0;
    for(int i=0; i<imax; i++) s += y[i];
    Wp[t] = s;
  }
  std::vector<double> tvec(horizon, 0);
  for(int i=0; i<horizon; i++) tvec[i] = i+1;
  
  return DataFrame::create(
    Named("t") = tvec,
    Named("m") = m,
    Named("I") = I,
    Named("S") = S, 
    Named("A") = A,
    Named("Y") = Y,
    Named("H") = H,
    Named("Hpercapita") = Hpercapita,
    Named("Wd") = Wd,
    Named("Wp") = Wp,
    Named("Wr") = Wp
  );
  
  
}
