data {
  int N;    //number of genes
  int K;    //number of cells
  int  y[N,K]; //read count data matrix
  real l[K]; 
  real B[N];
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
 real<lower=0> kon[N]; 
 real<lower=0> phi[N];
 real<lower=0> koff[N];
   //kon and koff are beta(kon,koff) parameters 
 matrix<lower=0,upper=1>[N,K] pi;  //bursting 
 vector<lower=0>[N] s;  //cell size parameters. gene specific expression rate, should be cell specific?
}
transformed parameters {
  matrix<lower=0>[N,K] lambda;  // read count poisson paramters 
  real<lower=0> prekoff[N];
  for(n in 1:N){
  prekoff[n]=(kon[n]*kon[n]*B[n]+kon[n]*B[n])/(1-kon[n]*B[n]);
  }
  for(n in 1:N){
    for(k in 1:K)
    lambda[n,k]=pi[n,k]*s[n]*l[k];
  }
  
}

model {
    s~gamma(100,3);
    kon~gamma(10,10);
    phi~gamma(0.0001,0.0001);
    for(n in 1:N){
   koff[n] ~ gamma(1/(phi[n]*phi[n]),prekoff[n]*phi[n]*phi[n]);
  }
    for(n in 1:N){
   pi[n,:] ~beta(kon[n],koff[n]);
  }
  
   // read count data generation process
    for(n in 1:N){
    for(k in 1:K)
      target += poisson_lpmf(y[n,k] | lambda[n,k]);
  }
}
