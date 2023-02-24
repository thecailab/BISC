data {
  int N;    //number of genes
  int K;    //number of cells
  int  y[N,K]; //read count data matrix
  real l[K]; //library size
  real tau;
 // BCV as input data
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
   real x0;
  real bcv[N,K];
 matrix<lower=0>[N,K] lambda;  
 real<lower=0> kon[N];    
 real<lower=0> koff[N];  //kon and koff are beta(kon,koff) parameters 
 matrix<lower=0,upper=1>[N,K] pi;  //bursting 
 vector<lower=0>[N] s;  //cell size parameters. gene specific expression rate, should be cell specific?
}
transformed parameters {
  matrix<lower=0>[N,K] lambda1;  // read count poisson paramters 
  matrix<lower=0,upper=1>[N,K] p; //zero-inflated ber(p)
  for(n in 1:N){
    for(k in 1:K)
    lambda1[n,k]=pi[n,k]*s[n]*l[k];
  }
  for(n in 1:N){
    for(k in 1:K)
    p[n,k]=1/(1+exp(tau*(log(lambda[n,k])-x0)));
  }
}
model {
  tau~normal(0,50);
  x0~normal(0,50);
    s~gamma(100,3);
    kon~gamma(2,0.1);
    koff~gamma(2,0.1);
    for(n in 1:N){
   pi[n,:] ~beta(kon[n],koff[n]);
  }
    for(n in 1:N){
    for(k in 1:K)
    lambda[n,k]~gamma(1/(bcv[n,k]*bcv[n,k]),1/(lambda1[n,k]*bcv[n,k]*bcv[n,k]));
  }
   // read count data generation process
    for(n in 1:N){
    for(k in 1:K)
      if(y[n,k]==0)
      target += log_sum_exp(bernoulli_lpmf(1 | p[n,k]),bernoulli_lpmf(0 | p[n,k]) + poisson_log_lpmf(y[n,k] | lambda[n,k]));
      else
      target += bernoulli_lpmf(0 | p[n,k])+ poisson_log_lpmf(y[n,k] | lambda[n,k]);
  }
}
