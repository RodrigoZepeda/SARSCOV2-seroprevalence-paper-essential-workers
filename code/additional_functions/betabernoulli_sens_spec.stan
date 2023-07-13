// ANALYSIS FOR BINARY VALUES THAT MIGHT BE ADJUSTED FOR SENSITIVITY/SPECIFICITY
// Checked 30 ago 2022
// --------------------------------------------------------------------------------
// Rodrigo Zepeda Tello
// email: rodrigo.zepeda®imss.gob.mx [change ® for appropriate similar symbol]
// 
// Model:
// -------------------------------------------------------------------------------
//
// w .- Work Center
// c .- Cycle of survey
// Cases(w,c) ~ Binomial(n(w,c), p(w,c))
// where
// p(w,c) = delta*¶(w,c) + (1 - gamma)*¶(w,c)
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
  int<lower=0>  Ncycles;                     //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                    //Número de centers de trabajo analizados
  int<lower=0>  cases[Ncenters,Ncycles];     //Número de individuos positivos por centro y cycle
  int<lower=0>  ns[Ncenters,Ncycles];        //Número de individuos muestreados por centro y cycle
  real<lower=0> pondef[Ncenters,Ncycles];    //Ponderadores asociados por cycle/centro
  real<lower=0,upper=1> sens;                //Sensibilidad de la prueba
  real<lower=0,upper=1> spec;                //Especificidad de la prueba
  real<lower=0> prior_alpha;                 //A priori para alpha
  real<lower=0> prior_beta;                  //A priori para beta
  real<lower=0> sigma_hiper;
}

parameters {
  
  //Variance of alpha and beta
  real<lower=0> sigmasq;
  real<lower=0> sigma_1;
  
  //Probability adjusted and unadjusted
  real<lower=0,upper=1> prevalence_positive[Ncenters,Ncycles];
  
  //Hierarchical distribution of p
  real<lower=0> alpha_centro_cycle[Ncenters,Ncycles];
  real<lower=0> beta_centro_cycle[Ncenters,Ncycles];
  real<lower=0> alpha_cycle[Ncycles];
  real<lower=0> beta_cycle[Ncycles];
  
}

transformed parameters {
  
  //Proportion of positive as well as prevalence proportions
  matrix[Ncenters,Ncycles] proportion_positive;
  vector[Ncycles] proportion_positive_cycle;
  vector[Ncycles] prevalence_positive_cycle;
  
  //Adjusted by sensititivy and specificity
  for (cycle in 1:Ncycles){
    for (center in 1:Ncenters){
      proportion_positive[center, cycle] = sens*prevalence_positive[center,cycle] + 
              (1 - spec)*(1 - prevalence_positive[center,cycle]);
    }
  }
      
  for (cycle in 1:Ncycles){
    
    //Start with 0 and sum the proportions
    proportion_positive_cycle[cycle] = 0;
    prevalence_positive_cycle[cycle] = 0;
    
    //Weighted prevalence
    for (center in 1:Ncenters){
      proportion_positive_cycle[cycle] += proportion_positive[center, cycle]*pondef[center,cycle];
      prevalence_positive_cycle[cycle] += prevalence_positive[center, cycle]*pondef[center,cycle];
    }
  }
}

model {
  //Prior variance
  sigmasq ~ normal(sigma_hiper, 0.0001);
  sigma_1 ~ normal(sigma_hiper, 0.0001);

  
  for (center in 1:Ncenters){
    for (cycle in 1:Ncycles){
      
      //Cases are binomial
      cases[center, cycle] ~ binomial(ns[center, cycle], proportion_positive[center, cycle]);
      
      //The adjusted has a beta prior
      prevalence_positive[center,cycle] ~ beta(alpha_centro_cycle[center, cycle], 
                          beta_centro_cycle[center, cycle]);
      
      //The hyperparameters
      alpha_centro_cycle[center, cycle] ~ normal(alpha_cycle[cycle], sigmasq);
      beta_centro_cycle[center, cycle]  ~ normal(beta_cycle[cycle], sigmasq);
      
      if (cycle > 1){
        alpha_cycle[cycle] ~ normal(alpha_cycle[cycle - 1], sigma_1);
        beta_cycle[cycle]  ~ normal(beta_cycle[cycle - 1], sigma_1);
      } else {
        alpha_cycle[cycle] ~ normal(prior_alpha, sigma_1);
        beta_cycle[cycle]  ~ normal(prior_beta, sigma_1);
      }
    }
  }
}
