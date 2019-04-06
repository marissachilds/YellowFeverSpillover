data {
    int<lower=0> n; //Number of observations
    int<lower=0> Y[n]; // number of mosquitoes caught at each site
    real<lower=0> r[n]; // distance from sites to release location
    real<lower=0> t[n]; // time spent at different release locations 
}
    
parameters {
    real<lower=0> beta; 
    real alpha;
    real<lower=0> k;
}
  
model {
  // Define the priors
    alpha ~ uniform(-100, 100);
    beta  ~ uniform(0, 100);
    k     ~ uniform(0, 100);
    // Define the likelihood
    for(i in 1:n){
      Y[i] ~ neg_binomial_2(t[i]*exp(alpha)*exp(-r[i]*beta), k);
    }
}

