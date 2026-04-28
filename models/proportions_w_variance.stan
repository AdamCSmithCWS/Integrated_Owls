
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n_regions;
  array[n_regions] int<lower=0,upper=n_regions> regions;
  
  vector[n_regions] mu_obs;
  vector[n_regions] sigma_obs;
  
  vector[n_regions] L;
  
  vector[n_regions] U;
  
}

// mu is the true mean abundance given the approximate variance and the upper and lower limits
parameters {
  vector[n_regions] mu;
}



//Truncated normal distribution of variation 
model {
  for(i in regions){
  mu_obs[i] ~ normal(mu[i], sigma_obs[i]) T[L[i], U[i]];
  }
}

generated quantities {

  vector[n_regions] p;

  for(i in regions){

    p[i] = mu[i]/(sum(mu[1:n_regions])-mu[i]); 
  }
  
}

