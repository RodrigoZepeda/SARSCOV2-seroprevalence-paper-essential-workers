// ANALYSIS FOR BINARY VALUES THAT MIGHT BE ADJUSTED FOR SENSITIVITY/SPECIFICITY
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
  int<lower=0>  Ncycles;                         //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                        //Número de centers de trabajo analizados
  int<lower=0>  Nweeks;                          //Número de semanas epidemiológicas
  int<lower=0>  cases[Nweeks, Ncenters,Ncycles]; //Número de individuos positivos por centro y cycle
  int<lower=0>  ns[Ncenters,Ncycles];    //Número de individuos muestreados por centro y cycle
  real<lower=0> pondef[Ncenters,Ncycles];        //Ponderadores asociados por cycle/centro
  real<lower=0> prior_alpha;                     //A priori para alpha
  real<lower=0> prior_beta;                      //A priori para beta
  real<lower=0> sigma_hiper;
}

parameters {
  
  //Variance of alpha and beta
  //real<lower=0> sigmasq;
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
  
  //Probability adjusted and unadjusted
  real<lower=0,upper=1> proportion_positive[Nweeks, Ncenters,Ncycles];
  
  //Hierarchical distribution of p
  real<lower=0> alpha_week_centro_cycle[Nweeks, Ncenters,Ncycles];
  real<lower=0> beta_week_centro_cycle[Nweeks, Ncenters,Ncycles];
  real<lower=0> alpha_week_center[Nweeks, Ncycles];
  real<lower=0> beta_week_center[Nweeks, Ncycles];
  //real<lower=0> alpha_week[Nweeks];
  //real<lower=0> beta_week[Nweeks];
  
}

model {
  //Prior variance
  //sigmasq ~ normal(sigma_hiper, 0.0001);
  sigma_1 ~ normal(sigma_hiper, 0.0001);
  sigma_2 ~ normal(sigma_hiper, 0.0001);
  
  for (week in 1:Nweeks){
    for (cycle in 1:Ncycles){
      for (center in 1:Ncenters){
      
        //Cases are binomial
        cases[week, center, cycle] ~ binomial(ns[center, cycle], 
            proportion_positive[week, center, cycle]);
      
        //The adjusted has a beta prior
        proportion_positive[week, center,cycle] ~ beta(
          alpha_week_centro_cycle[week, center, cycle], 
          beta_week_centro_cycle[week, center, cycle]
        );
        
        alpha_week_centro_cycle[week, center, cycle] ~ normal(alpha_week_center[week, cycle], sigma_1);
        beta_week_centro_cycle[week, center, cycle]  ~ normal(beta_week_center[week, cycle], sigma_1);
        
      }
      
      if (cycle > 1){
        alpha_week_center[week, cycle] ~ normal(alpha_week_center[week, cycle - 1], sigma_2);
        beta_week_center[week, cycle]  ~ normal(alpha_week_center[week, cycle - 1], sigma_2);
      } else {
        alpha_week_center[week, cycle] ~ normal(prior_alpha, sigma_2);
        beta_week_center[week, cycle]  ~ normal(prior_beta, sigma_2);
      }
    }
  }
}

generated quantities {
  
  real<lower=0,upper=1> proportion_week_cycle[Nweeks, Ncycles];
  real<lower=0,upper=1> proportion_week[Nweeks];
  
  for (week in 1:Nweeks){
    proportion_week[week] = 0;
    
    for (cycle in 1:Ncycles){
    
     proportion_week_cycle[week,cycle] = 0;
    
     for (center in 1:Ncenters){
      proportion_week_cycle[week,cycle] += 
        proportion_positive[week, center, cycle]*pondef[center,cycle];
      }
      
      proportion_week[week] += proportion_week_cycle[week,cycle]/Ncycles;
    }
    
  }
}
