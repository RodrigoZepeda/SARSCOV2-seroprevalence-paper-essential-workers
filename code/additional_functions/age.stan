// ANALYSIS FOR CONTINUOUS NORMALLY DISTRIBUTED VARIABLES
// --------------------------------------------------------------------------------
// Rodrigo Zepeda Tello
// email: rodrigo.zepeda®imss.gob.mx [change ® for appropriate similar symbol]
// 
// Model:
// -------------------------------------------------------------------------------
//
// w .- Work Center
// c .- Cycle of survey
// Age(w,c) ~ Normal(n(w,c), p(w,c))
// where
// sigma(w,c) ~ Cauchy(0, 2.5)
// with delta = sensitivity and gamma = specificity of the test ( = 1 for perfect test).
// The value ¶ represents the adjusted (true) prevalence.
// ¶(w,c) ~ Beta(alpha(w,c), beta(w,c))
// where: 
// alpha(w,c) ~ Normal(A(c), sigma)
// beta(w,c)  ~ Normal(B(c), sigma)
// represent the work center specific parameters and A(c), B(c) the prevalence for the cycle
// given by;
// A(c) = A(c - 1) + epsilon_A
// B(c) = B(c - 1) + epsilon_B
// with epsilon_A and epsilon_B being random (normally distributed) errors with variance sigma_hiper. 
// and A(1) ~ Normal(alpha_prior, 0.001); B(1) = ~ Normal(beta_prior, 0.001)
//
// --------------------------------------------------------------------------------

data {
  int<lower=0>  N;                           //Número de mediciones de la encuesta
  int<lower=0>  Ncycles;                     //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                    //Número de centers de trabajo analizados
  int<lower=0>  ages[N,3];                   //Edad de los individuios
  real<lower=0> pondef[Ncenters,Ncycles];    //Ponderadores asociados por cycle/centro
  real<lower=0> prior_mu;                    //A priori para alpha
  real<lower=0> prior_sigma;                 //A priori para beta
  real<lower=0> sigma_hiper;
}

transformed data {
  int<lower=0> age_col    = 1; //Age column
  int<lower=0> center_col = 2; //Center column
  int<lower=0> cycle_col  = 3;  //Cycle column
}

#La solución usar un vector de diferente longitud

parameters {
  
  //Variance of alpha and beta
  real<lower=0> sigmasq;
  real<lower=0> sigma_1;
  
  //Hierarchical distribution of p
  real<lower=0> mu_centro_cycle[Ncenters,Ncycles];
  real<lower=0> sigma_centro_cycle[Ncenters,Ncycles];
  real<lower=0> mu_cycle[Ncycles];
  real<lower=0> sigma_cycle[Ncycles];
  
}

transformed parameters {
  
  //Proportion of positive as well as prevalence proportions
  vector[Ncycles] mu_all_centers;
  
  
  for (cycle in 1:Ncycles){
    
    //Start with 0 and sum the proportions
    mu_all_centers[cycle] = 0;
    
    //Weighted prevalence
    for (center in 1:Ncenters){
      mu_all_centers[cycle] += mu_centro_cycle[center, cycle]*pondef[center,cycle];
    }
  }
}

model {
  //Prior variance
  sigmasq ~ normal(sigma_hiper, 0.0001);
  sigma_1 ~ normal(sigma_hiper, 0.0001);

  for (n in 1:N){
    
    //Cases are binomial
    ages[n, age_col] ~ normal(mu_centro_cycle[ages[n, center_col], ages[n, cycle_col]], sigma_centro_cycle[ages[n, center_col], ages[n, cycle_col]]);
      
  }
  
  for (center in 1:Ncenters){
    for (cycle in 1:Ncycles){
      
      //The hyperparameters
      mu_centro_cycle[center, cycle]     ~ normal(mu_cycle[cycle], sigmasq);
      sigma_centro_cycle[center, cycle]  ~ normal(sigma_cycle[cycle], sigmasq);
      
      if (cycle > 1){
        mu_cycle[cycle]     ~ normal(mu_cycle[cycle - 1], sigma_1);
        sigma_cycle[cycle]  ~ normal(sigma_cycle[cycle - 1], sigma_1);
      } else {
        mu_cycle[cycle]     ~ normal(prior_mu, sigma_1);
        sigma_cycle[cycle]  ~ normal(prior_sigma, sigma_1);
      }
    }
  }
}
