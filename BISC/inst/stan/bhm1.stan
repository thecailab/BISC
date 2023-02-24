data {
  int N;    //number of genes
  int K;    //number of cells
  int  y[N,K]; //read count data matrix
  real bcv[N,K];// BCV as input data
  real l[K]; 
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
 matrix<lower=0>[N,K] lambda;
 real<lower=0> kon[N];    
 real<lower=0> koff[N];  //kon and koff are beta(kon,koff) parameters 
 matrix<lower=0,upper=1>[N,K] pi;  //bursting 
 vector<lower=0>[N] s;  //cell size parameters. gene specific expression rate, should be cell specific?
}
transformed parameters {
  matrix<lower=0>[N,K] lambda1;  // read count poisson paramters 
  for(n in 1:N){
    for(k in 1:K)
    //lambda1[n,k]=pi[n,k]*s[n]*l[k];
	lambda1[n,k]=pi[n,k]*s[n];
  }
  
}

model {
    s~gamma(100,3);
    kon~gamma(10,3);
    koff~gamma(10,3);
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
      target += poisson_lpmf(y[n,k] | lambda[n,k]);
  }
}
