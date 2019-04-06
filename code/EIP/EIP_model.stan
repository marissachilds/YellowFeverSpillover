data {
    int<lower=0> k; //Number of observations
    real temp[k];
    real<lower=0> X1[k];
    int<lower=0> N1[k];
    real<lower=0> X2[k];
    int<lower=0> N2[k];
}
    
parameters {
    real<lower=0, upper = 1> M_max_inf; 
    real M_T0;
    real M_Tm;
    real mu0;
    real muT;
    real sigma0;
    real sigmaT;
}
  
transformed parameters {
    real c; // scaling parameter in M,
    c = -4*M_max_inf*(M_Tm-M_T0)^(-2); 
}
    
model {
    real M[k]; // max proportion infectious
    real mu[k]; // log of EIP50 = log of median of EIP distribution
    real sigma[k]; // the variance of EIP as a function of temperature
    real Fx1[k]; // probability of infection by time X
    real Fx2[k]; // probability of infection by time Y

    // Priors
    M_max_inf ~ uniform(0, 1);
    M_T0 ~ normal(10, 10); 
    M_Tm ~ normal(40, 10); 
    mu0 ~ normal(2, 0.5); 
    muT ~ normal(0, 0.05);
    sigma0 ~ normal(-1,0.5);
    sigmaT ~ normal(0,0.1);
       
    for(i in 1:k){
      // Convert temperature into parameters
      M[i] = ((M_T0 < temp[i] && temp[i] < (M_Tm)) ? 
        c*(temp[i]-M_T0).*(temp[i]-(M_Tm)) : 0);
      mu[i] = temp[i]*muT + mu0;
      sigma[i] = exp(sigma0 + sigmaT*temp[i]);
      // Log normal F
      if(X1[i] == 0) Fx1[i] = 0;
      else Fx1[i] = 0.5 + 0.5*erf( (log(X1[i]) - mu[i])/(sqrt(2)*sigma[i]));
      Fx2[i] = 0.5 + 0.5*erf( (log(X2[i]) - mu[i])/(sqrt(2)*sigma[i]) );
      }
    
    // likelihood
    for(i in 1:k){
      if(N2[i] == 0){ // if its a zero, then there was never transmission
        target += N1[i]*log(1-M[i]*Fx1[i]);
      }
      else{ // if its not zero, there's a min and max
        target += N1[i]*log(1-M[i]*Fx1[i]) + log(1-(1-M[i]*Fx2[i])^N2[i]*(1-M[i]*Fx1[i])^(N1[i]-N2[i]));
      }
    }
}

