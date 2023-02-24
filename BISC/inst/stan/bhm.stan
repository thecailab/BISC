
data {
  int N;    //number of genes
  int K;    //number of cells
  int  y[N,K]; //read count data matrix
  real l[K]; 
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
 real<lower=0> kon[N];    
 real<lower=0> koff[N];  //kon and koff are beta(kon,koff) parameters 
 matrix<lower=0,upper=1>[N,K] pi;  //bursting 
 vector<lower=0>[N] s;  //cell size parameters. gene specific expression rate, should be cell specific?
}
transformed parameters {
  matrix<lower=0>[N,K] lambda;  // read count poisson paramters 
  for(n in 1:N){
    for(k in 1:K)
    //lambda[n,k]=pi[n,k]*s[n]*l[k];
    lambda[n,k]=pi[n,k]*s[n];
  }
}

model {
    s~gamma(10,0.03);
    kon~gamma(1,0.01);
    koff~gamma(1,0.01);
    for(n in 1:N){
   pi[n,:] ~beta(kon[n],koff[n]);
  }
   // read count data generation process
    for(n in 1:N){
    for(k in 1:K)
      target += poisson_lpmf(y[n,k] | lambda[n,k]);
  }
}

