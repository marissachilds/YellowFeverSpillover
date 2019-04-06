data {
    int<lower=0> k; //Number of observations
    real L[k]; // Mean lifespan
    real temp[k]; // temperature 
    real loc[k]; // location (field vs lab)
}
    
parameters {
    real<lower = 0> max_L_lab; // Maximum lifespan in the lab
    real<lower = 0> max_L_field; // Maximum lifespan in the field
    real T0; // Critical thermal minimum, assumed same in lab and field
    real Tm; // Critical thermal maximum, assumed same in lab and field
    real<lower=0> sigma; // Standard deviation
}
  
transformed parameters {
    real c_lab;  
    real c_field; 
    c_lab = -4*max_L_lab*(Tm-T0)^(-2); // convert lab max lifespan to the function coefficient
    c_field = -4*max_L_field*(Tm-T0)^(-2); // convert field max lifespan to the function coefficient
}
    
model {
    // uninformative priors for lifespan
    max_L_lab ~ uniform(0, 200); 
    max_L_field ~ uniform(0, 200);
    // Slightly more limited priors for thermal limits 
    T0 ~ normal(10, 5);
    Tm ~ normal(40, 5);
    sigma ~ exponential(0.1);

    // Likelihood based on location of observations
    for(i in 1:k){
      if(loc[i] == 1){
        L[i] ~ normal(c_lab*(temp[i]-T0)*(temp[i]-Tm), sigma);
      }
      else{
        L[i] ~ normal(c_field*(temp[i]-T0)*(temp[i]-Tm), sigma);
      }
    }
}

