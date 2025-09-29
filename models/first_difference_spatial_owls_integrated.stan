// This is a Stan implementation of the first difference model that shares information among strata on the annual differences
// This is an elaboration of the model used in Link and Sauer
// 

// iCAR function, from Morris et al. 2019
// Morris, M., K. Wheeler-Martin, D. Simpson, S. J. Mooney, A. Gelman, and C. DiMaggio (2019). 
// Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. 
// Spatial and Spatio-temporal Epidemiology 31:100301.

 functions {
   real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
     return -0.5 * dot_self(bb[n1] - bb[n2])
       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
  }
 }




data {
  
  ///BBS only data
  int<lower=1> n_sites_owl;
  int<lower=1> n_sites_bbs;
  int<lower=1> n_counts_bbs;
  int<lower=1> n_counts_owl;
  int<lower=1> n_years;
  int<lower=1> n_strata;
  int<lower=1> n_protocols;
  int<lower=1> ebird_year; // index of year for the eBird relative abundance surface
  array[ebird_year-1] int<lower=1> yrev; // reverse year vector (ebird_year-1):1 (necessary for Stan indexing)
  
  
  array[n_counts_bbs] int<lower=0> count_bbs;              // count observations
  array[n_counts_bbs] int<lower=1> strat_bbs;               // strata indicators
  array[n_counts_bbs] int<lower=1> year_bbs; // year index
  array[n_counts_bbs] int<lower=1> site_bbs; // site index


  array[n_counts_owl] int<lower=0> count_owl;              // count observations
  array[n_counts_owl] int<lower=1> strat_owl;               // strata indicators
  array[n_counts_owl] int<lower=1> year_owl; // year index
  array[n_counts_owl] int<lower=1> site_owl; // site index
  array[n_counts_owl] int<lower=1> proto; // site index
  vector[n_counts_owl] off_set; // site index

  vector[n_strata] log_mean_rel_abund; // log-transformed mean relative abund in stratum

  // a vector of zeros to fill fixed beta values for fixed_year
  vector[n_strata] zero_betas;
  //array[n_years_m1] int<lower=0> y_2020; //indicators for 2020 = 0 if 2020 and missing if fixed_year

  
  //data for spatial iCAR among strata
  int<lower=1> n_edges;
  array [n_edges] int<lower=1, upper=n_strata> node1;  // node1[i] adjacent to node2[i]
  array [n_edges] int<lower=1, upper=n_strata> node2;  // and node1[i] < node2[i]

  // Extra Poisson variance options
  int<lower=0,upper=1> heavy_tailed; //indicator if extra poisson variance should be t-distributed or normal (yes = 1, no = 0 and therefore normal)
  int<lower=0,upper=1> calc_nu; //indicator if nu should be calculated (yes = 1, no = 0)
  int<lower=0,upper=1> use_pois; //indicator if count variation should be based on over-dispersed Poisson (if ==1) or Negative binomial (if == 0)
  
  // loo or CV calculations
  int<lower=0,upper=1> calc_log_lik; //indicator if log_lik should be calculated (log_lik for all obs to support loo = 1, no log-lik = 0)
  int<lower=0,upper=1> calc_cv; //indicator if CV should be calculated (CrossValidation = 1, no CV = 0)
  // CV folds - if calc_cv == 1 then the following values define the training and test sets
  int<lower=1, upper=n_counts_bbs> n_train_bbs; //
  int<lower=1, upper=n_counts_bbs> n_test_bbs; //
  array[n_train_bbs] int<lower=1, upper=n_counts_bbs> train_bbs; // indices of counts to include in train data
  array[n_test_bbs] int<lower=1, upper=n_counts_bbs> test_bbs; // indices of counts to include in test data
 
  int<lower=1, upper=n_counts_owl> n_train_owl; //
  int<lower=1, upper=n_counts_owl> n_test_owl; //
  array[n_train_owl] int<lower=1, upper=n_counts_owl> train_owl; // indices of counts to include in train data
  array[n_test_owl] int<lower=1, upper=n_counts_owl> test_owl; // indices of counts to include in test data
  

}

transformed data {
   //These statements split the data into training and testing sets for cross-validation
   // if calc_cv == 0 (and so no CV is required), then n_train = ncount and all data are included in the training set
   // in that case, n_test = 1, but all values of xxx_te are ignored for the remainder of the model
  
  int n_years_m1 = n_years-1;
  
  array[n_train_bbs] int<lower=0> count_bbs_tr = count_bbs[train_bbs];              // count observations
  array[n_train_bbs] int<lower=1> strat_bbs_tr = strat_bbs[train_bbs];               // strata indicators
  array[n_train_bbs] int<lower=1> year_bbs_tr = year_bbs[train_bbs]; // year index
  array[n_train_bbs] int<lower=1> site_bbs_tr = site_bbs[train_bbs]; // site index

  array[n_test_bbs] int<lower=0> count_bbs_te = count_bbs[test_bbs];              // count observations
  array[n_test_bbs] int<lower=1> strat_bbs_te = strat_bbs[test_bbs];               // strata indicators
  array[n_test_bbs] int<lower=1> year_bbs_te = year_bbs[test_bbs]; // year index
  array[n_test_bbs] int<lower=1> site_bbs_te = site_bbs[test_bbs]; // site index

  array[n_train_owl] int<lower=0> count_owl_tr = count_owl[train_owl];              // count observations
  array[n_train_owl] int<lower=1> strat_owl_tr = strat_owl[train_owl];               // strata indicators
  array[n_train_owl] int<lower=1> year_owl_tr = year_owl[train_owl]; // year index
  array[n_train_owl] int<lower=1> site_owl_tr = site_owl[train_owl]; // site index
  array[n_train_owl] int<lower=1> proto_tr = proto[train_owl]; // site index
  vector[n_train_owl] off_set_tr = off_set[train_owl]; // site index

  array[n_test_owl] int<lower=0> count_owl_te = count_owl[test_owl];              // count observations
  array[n_test_owl] int<lower=1> strat_owl_te = strat_owl[test_owl];               // strata indicators
  array[n_test_owl] int<lower=1> year_owl_te = year_owl[test_owl]; // year index
  array[n_test_owl] int<lower=1> site_owl_te = site_owl[test_owl]; // site index
  array[n_test_owl] int<lower=1> proto_te = proto[test_owl]; // site index
  vector[n_test_owl] off_set_te = off_set[test_owl]; // site index

  
  
  
}


parameters {
  vector[n_train_bbs*use_pois] noise_bbs_raw;             // over-dispersion if use_pois == 1
  vector[n_train_owl*use_pois] noise_owl_raw;             // over-dispersion if use_pois == 1
 
  real OWL; // OWL intercept 
  real BBS; // BBS intercept 
 
  vector[n_sites_bbs] ste_bbs_raw;   // site (route) effects
  vector[n_sites_owl] ste_owl_raw;   // site (route) effects

  real<lower=0> sdste_bbs;    // sd of site (route) effects
  real<lower=0> sdste_owl;    // sd of site effects

  vector[n_protocols] protocol_raw;   // protocol effects for owls
  real<lower=0> sdprotocol;    // sd of overall protocol

  real<lower=0> sdnoise_bbs;    // sd of over-dispersion, if use_pois == 1
  real<lower=0> sdnoise_owl;    // sd of over-dispersion, if use_pois == 1
  real<lower=3> nu_bbs; // df of t-distribution, if calc_nu == 1, > 3 so that it has a finite mean, variance, kurtosis
  real<lower=3> nu_owl; // df of t-distribution > 3 so that it has a finite mean, variance, kurtosis

  vector<lower=0>[n_years] sdbeta;    // sd of annual changes among strata 
 // real<lower=0> sdbeta;    // sd of annual changes among strata 
  real<lower=0> sdBETA;    // sd of overall annual changes

  vector[n_years_m1] BETA_raw;//_hyperparameter of overall annual change values - "differences" between years 
  matrix[n_strata,n_years_m1] beta_raw;         // strata level parameters
  

   

}

transformed parameters { 
  vector[n_train_bbs] E_bbs;           // log_scale additive likelihood
  vector[n_train_owl] E_owl;           // log_scale additive likelihood

  matrix[n_strata,n_years] beta;         // strata-level mean differences (0-centered deviation from continental mean BETA)
  matrix[n_strata,n_years] yeareffect;  // matrix of estimated annual values of trajectory
  vector[n_years_m1] BETA; // annual estimates of continental mean differences (n_years - 1, because middle year is fixed at 0)
  vector[n_years] YearEffect;

  vector[n_protocols] protocol;

  real<lower=0> phi_bbs; //transformed sdnoise if use_pois == 0 (and therefore Negative Binomial)
  real<lower=0> phi_owl; //transformed sdnoise
 
  if(use_pois){
    phi_bbs = 0;
    phi_owl = 0;
  }else{
    phi_bbs = 1/sqrt(sdnoise_bbs); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    phi_owl = 1/sqrt(sdnoise_owl); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  }
  
  protocol = (sdprotocol*protocol_raw); //
  

  BETA = sdBETA * BETA_raw;

  beta[,ebird_year] = zero_betas; //set to zero but never used
  yeareffect[,ebird_year] = zero_betas; //fixed at zero
  YearEffect[ebird_year] = 0; //fixed at zero

// first half of time-series - runs backwards from ebird_year
  for(t in yrev){
    beta[,t] = (sdbeta[t] * beta_raw[,t]) + BETA[t];
    yeareffect[,t] = yeareffect[,t+1] - beta[,t];
    YearEffect[t] = YearEffect[t+1] - BETA[t]; // hyperparameter trajectory interesting to monitor but no direct inference
  }
// second half of time-series - runs forwards from ebird_year
if(n_years - ebird_year){ // only used if ebird_year is not equal to n_years, i.e., ebird_year != final year of the BBS data
for(t in (ebird_year+1):n_years){

    beta[,t] = (sdbeta[t] * beta_raw[,t-1]) + BETA[t-1];//t-1 indicators to match dimensionality
    yeareffect[,t] = yeareffect[,t-1] + beta[,t];
    YearEffect[t] = YearEffect[t-1] + BETA[t-1]; 
  }
}


  for(i in 1:n_train_bbs){
    real noise;
    //real obs = sdobs*obs_raw[observer_tr[i]];
    real ste = sdste_bbs*ste_bbs_raw[site_bbs_tr[i]]; // site intercepts are zero-centered for BBS
    if(use_pois){
    noise = sdnoise_bbs*noise_bbs_raw[i];
    }else{
    noise = 0;
    }
    
    E_bbs[i] =  BBS + yeareffect[strat_bbs_tr[i],year_bbs_tr[i]] + ste + noise;
  }
  
  
  
    for(i in 1:n_train_owl){
    real noise;
    //real obs = sdobs*obs_raw[observer_tr[i]];
    real ste = sdste_owl*ste_owl_raw[site_owl_tr[i]]; // site intercepts are zero-centered for BBS
    if(use_pois){
    noise = sdnoise_owl*noise_owl_raw[i];
    }else{
    noise = 0;
    }
    
    E_owl[i] =  OWL + yeareffect[strat_owl_tr[i],year_owl_tr[i]] + ste + off_set_tr[i] + protocol[proto_tr[i]] + noise;
  }
  

  
  }
  
  
  
model {
  nu_bbs ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  nu_owl ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
 
  if(use_pois){
    if(heavy_tailed){
    if(calc_nu){
      noise_bbs_raw ~ student_t(nu_bbs,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
      noise_owl_raw ~ student_t(nu_owl,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
   }else{
      noise_bbs_raw ~ student_t(3,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
      noise_owl_raw ~ student_t(3,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
    }
   }else{
    noise_bbs_raw ~ std_normal();//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
    noise_owl_raw ~ std_normal();//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
    }
  }
  if(use_pois){
  sdnoise_bbs ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  sdnoise_owl ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  }else{
  sdnoise_bbs ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  sdnoise_owl ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  }  
//  sdobs ~ normal(0,0.3); // informative prior on scale of observer effects - suggests observer variation larger than 3-4-fold differences is unlikely
  sdbeta ~ student_t(3,0,0.05); // prior on sd among strata in the yearly differences
  sdBETA ~ student_t(3,0,0.5); // prior on sd of mean hyperparameter time-series of range-wide mean
  BETA_raw ~ std_normal();// prior on fixed effect mean intercept

  // bbs
  sdste_bbs ~ student_t(3,0,1); //prior on sd of site effects
  ste_bbs_raw ~ std_normal();//site effects
  sum(ste_bbs_raw) ~ normal(0,0.001*n_sites_bbs); //
  BBS ~ std_normal();//bbs mean
  
  //owl 
  protocol_raw ~ std_normal();//site effects
  sum(protocol_raw) ~ normal(0,0.001*n_protocols); //
  OWL ~ std_normal();// prior on fixed effect mean of owl surveys
  sdste_owl ~ student_t(3,0,1); //prior on sd of site effects
  ste_owl_raw ~ std_normal();//site effects
  sum(ste_owl_raw) ~ normal(0,0.001*n_sites_owl); //constraint 





// Shared annual change values, spatially explicit and 
// fixed at zero in year of eBird abundance surface
for(t in 1:(n_years_m1)){
   //   if(ebird_year-t){ // zero/false for the ebird year
 
    beta_raw[,t] ~ icar_normal(n_strata, node1, node2);
    // }else{
    //  beta_raw[,t] ~ normal(0,0.001); // arbitrarily small because this is the ebird year and these are not included in the likelihood
    // 
    //  }
}

 
  
if(use_pois){
  count_bbs_tr ~ poisson_log(E_bbs); //vectorized count likelihood with log-transformation
  count_owl_tr ~ poisson_log(E_owl); //vectorized count likelihood with log-transformation
}else{
   count_bbs_tr ~ neg_binomial_2_log(E_bbs,phi_bbs); //vectorized count likelihood with log-transformation
   count_owl_tr ~ neg_binomial_2_log(E_owl,phi_owl); //vectorized count likelihood with log-transformation
 
}

}

 generated quantities {
// only apply to BBS

   array[n_strata,n_years] real<lower=0> n; //full annual indices
   vector<lower=0>[n_years] Hyper_N; // hyperparameter mean survey-wide population trajectory - only for the first difference model
  
   vector[n_counts_bbs*calc_log_lik] log_lik_bbs; // alternative value to track the observervation level log-likelihood
   vector[n_test_bbs*calc_cv] log_lik_cv_bbs; // alternative value to track the log-likelihood of the coutns in the test dataset

   vector[n_counts_owl*calc_log_lik] log_lik_owl; // alternative value to track the observervation level log-likelihood
   vector[n_test_owl*calc_cv] log_lik_cv_owl; // alternative value to track the log-likelihood of the coutns in the test dataset

 
  if(calc_log_lik){
  // potentially useful for estimating loo-diagnostics, such as looic
  if(use_pois){
  for(i in 1:n_train_owl){
   log_lik_owl[i] = poisson_log_lpmf(count_owl_tr[i] | E_owl[i]);
   }
     for(i in 1:n_train_bbs){
   log_lik_bbs[i] = poisson_log_lpmf(count_bbs_tr[i] | E_bbs[i]);
   }
  }else{
   for(i in 1:n_train_owl){
   log_lik_owl[i] = neg_binomial_2_log_lpmf(count_owl_tr[i] | E_owl[i] , phi_owl);
   } 
  for(i in 1:n_train_bbs){
   log_lik_bbs[i] = neg_binomial_2_log_lpmf(count_bbs_tr[i] | E_bbs[i] , phi_bbs);
   } 
  }
  }
  
  if(calc_cv){
    
    for(i in 1:n_test_bbs){
      
    real noise;
 //   real obs = sdobs*obs_raw[observer_te[i]];
    real ste = sdste_bbs*ste_bbs_raw[site_bbs_te[i]]; // site intercepts
   
   if(use_pois){
      if(heavy_tailed){
        if(calc_nu){
    noise = student_t_rng(nu_bbs,0,sdnoise_bbs);
        }else{
    noise = student_t_rng(3,0,sdnoise_bbs);
        }
      }else{
    noise = normal_rng(0,sdnoise_bbs);
      }
   
      
   log_lik_cv_bbs[i] = poisson_log_lpmf(count_bbs_te[i] | BBS + yeareffect[strat_bbs_te[i],year_bbs_te[i]] + ste + noise);
  
   }else{
     noise = 0;
  log_lik_cv_bbs[i] = neg_binomial_2_log_lpmf(count_bbs_te[i] | BBS + yeareffect[strat_bbs_te[i],year_bbs_te[i]] + ste + noise, phi_bbs);
  
   }
  
  }
  
  
  /// owls
   for(i in 1:n_test_owl){
      
    real noise;
 //   real obs = sdobs*obs_raw[observer_te[i]];
    real ste = sdste_owl*ste_owl_raw[site_owl_te[i]]; // site intercepts
   
   if(use_pois){
      if(heavy_tailed){
        if(calc_nu){
    noise = student_t_rng(nu_owl,0,sdnoise_owl);
        }else{
    noise = student_t_rng(3,0,sdnoise_owl);
        }
      }else{
    noise = normal_rng(0,sdnoise_owl);
      }
   
      
   log_lik_cv_owl[i] = poisson_log_lpmf(count_owl_te[i] | OWL + yeareffect[strat_owl_te[i],year_owl_te[i]] + ste + off_set_te[i] + proto_te[i] + noise);
  
   }else{
     noise = 0;
  log_lik_cv_owl[i] = neg_binomial_2_log_lpmf(count_owl_te[i] | OWL + yeareffect[strat_owl_te[i],year_owl_te[i]] + ste + off_set_te[i] + proto_te[i] + noise, phi_owl);
  
   }
  
  }
  
  }
  
   

// Annual indices of abundance - strata-level annual predicted counts


for(y in 1:n_years){

      for(s in 1:n_strata){

      n[s,y] = exp(yeareffect[s,y] + log_mean_rel_abund[s]);
        }

      Hyper_N[y] = exp(YearEffect[y]);
      

    }
  }





