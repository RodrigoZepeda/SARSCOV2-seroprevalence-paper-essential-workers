// ANALYSIS FOR BINARY VALUES THAT MIGHT BE ADJUSTED FOR SENSITIVITY/SPECIFICITY
// CONDITIONAL ON A SECOND VARIABLE
// --------------------------------------------------------------------------------
// Rodrigo Zepeda Tello
// email: rodrigo.zepeda®imss.gob.mx [change ® for appropriate similar symbol]
// 
// Model:
// -------------------------------------------------------------------------------
// w .- Work Center
// c .- Cycle of survey
// PrimaryCases(w,c) ~ Binomial(n(w,c), p(w,c))
// and then
// SecondaryCases+(w,c) ~ Binomial(PrimaryCases(w,c), q+(w,c))
// SecondaryCases-(w,c) ~ Binomial(n(w,c) - PrimaryCases(w,c), q-(w,c))
// where
// q+(w,c) = delta*¶+(w,c) + (1 - gamma)*¶+(w,c)
// q-(w,c) = delta*¶-(w,c) + (1 - gamma)*¶-(w,c)
// with delta = sensitivity and gamma = specificity of the test ( = 1 for perfect test).
// The value ¶ represents the adjusted (true) prevalence.
// ¶+(w,c) ~ Beta(alpha+(w,c), beta+(w,c))
// ¶-(w,c) ~ Beta(alpha-(w,c), beta-(w,c))
// and p(w,c) ~ Beta(alphap(w,c), betap(w,c))
// where: 
// alpha+(w,c), alpha-(w,c) ~ Normal(Aq(c), sigma)
// beta+(w,c), beta-(w,c)   ~ Normal(Bq(c), sigma)
// and
// alphap ~ Normal(Ap(c), sigma); betap ~ Normal(Bp(c), sigma)
// finally
// Aq(c) = Aq(c - 1) + epsilon_Aq
// Bq(c) = Bq(c - 1) + epsilon_Bq
// Ap(c) = Ap(c - 1) + epsilon_Ap
// Bp(c) = Bp(c - 1) + epsilon_Bp
// with epsilon_A and epsilon_B being random (normally distributed) errors. 
// and A(1) ~ Normal(alpha_prior, 100); B(1) = ~ Normal(beta_prior, 100)
//
// --------------------------------------------------------------------------------

data {
  int<lower=0>  Ncycles;                                    //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                                   //Número de centers de trabajo analizados
  int<lower=0>  primary_cases[Ncenters,Ncycles];            //Número de individuos positivos (primary, ITT) por centro y cycle
  int<lower=0>  secondary_positive_cases[Ncenters,Ncycles]; //Número de individuos positivos (secondary, IgG) con ITT por centro y cycle
  int<lower=0>  secondary_negative_cases[Ncenters,Ncycles]; //Número de individuos positivos (secondary, IgG) sin ITT por centro y cycle
  int<lower=0>  ns[Ncenters,Ncycles];                       //Número de individuos muestreados por centro y cycle para el primary
  real<lower=0> pondef[Ncenters,Ncycles];                   //Ponderadores asociados por cycle/centro
  real<lower=0,upper=1> sens;                               //Sensibilidad de la prueba
  real<lower=0,upper=1> spec;                               //Especificidad de la prueba
  real<lower=0> prior_alpha;                                //A priori para alpha
  real<lower=0> prior_beta;                                 //A priori para beta
  real<lower=0> sigma_hiper;
}

parameters {
  
  //Variance of alpha and beta
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> sigma3;
  real<lower=0> sigma4;
  real<lower=0> sigma5;
  
  //Probability adjusted and unadjusted
  real<lower=0,upper=1> prevalence_primary[Ncenters,Ncycles];
  real<lower=0,upper=1> prevalence_secondary_positive[Ncenters,Ncycles];
  real<lower=0,upper=1> prevalence_secondary_negative[Ncenters,Ncycles];
  
  //Hierarchical distribution of p
  real<lower=0> alpha_primary_centro_cycle[Ncenters,Ncycles];
  real<lower=0> beta_primary_centro_cycle[Ncenters,Ncycles];
  real<lower=0> alpha_secondary_positive_centro_cycle[Ncenters,Ncycles];
  real<lower=0> beta_secondary_positive_centro_cycle[Ncenters,Ncycles];
  real<lower=0> alpha_secondary_negative_centro_cycle[Ncenters,Ncycles];
  real<lower=0> beta_secondary_negative_centro_cycle[Ncenters,Ncycles];
  real<lower=0> alpha_secondary_cycle[Ncycles];
  real<lower=0> beta_secondary_cycle[Ncycles];
  real<lower=0> alpha_primary_cycle[Ncycles];
  real<lower=0> beta_primary_cycle[Ncycles];
  
}

transformed parameters {
  
  //Proportion of positive as well as prevalence proportions
  matrix[Ncenters,Ncycles] proportion_secondary_positive;
  matrix[Ncenters,Ncycles] proportion_secondary_negative;
  vector[Ncycles] proportion_primary_cycle;
  vector[Ncycles] proportion_secondary_cycle;
  vector[Ncycles] prevalence_secondary_cycle;
  vector[Ncycles] proportion_secondary_cycle_positive;
  vector[Ncycles] prevalence_secondary_cycle_positive;
  vector[Ncycles] proportion_secondary_cycle_negative;
  vector[Ncycles] prevalence_secondary_cycle_negative;
  vector[Ncycles] conditional_proportion_primary_secondary_positive;
  vector[Ncycles] conditional_proportion_primary_secondary_negative;
  vector[Ncycles] conditional_prevalence_primary_secondary_positive;
  vector[Ncycles] conditional_prevalence_primary_secondary_negative;
  
  //Adjusted by sensititivy and specificity
  for (cycle in 1:Ncycles){
    for (center in 1:Ncenters){
      
      proportion_secondary_positive[center, cycle] = sens*prevalence_secondary_positive[center,cycle] + 
              (1 - spec)*(1 - prevalence_secondary_positive[center,cycle]);
              
      proportion_secondary_negative[center, cycle] = sens*prevalence_secondary_negative[center,cycle] + 
              (1 - spec)*(1 - prevalence_secondary_negative[center,cycle]);
    }
  }
      
  for (cycle in 1:Ncycles){
    
    //Start with 0 and sum the proportions
    proportion_primary_cycle[cycle] = 0;
    
    proportion_secondary_cycle_positive[cycle] = 0;
    prevalence_secondary_cycle_positive[cycle] = 0;
    
    proportion_secondary_cycle_negative[cycle] = 0;
    prevalence_secondary_cycle_negative[cycle] = 0;
  
    //Weighted prevalence
    for (center in 1:Ncenters){
      
      proportion_primary_cycle[cycle] += prevalence_primary[center, cycle]*pondef[center,cycle];
      
      proportion_secondary_cycle_positive[cycle] += proportion_secondary_positive[center, cycle]*pondef[center,cycle];
      prevalence_secondary_cycle_positive[cycle] += prevalence_secondary_positive[center, cycle]*pondef[center,cycle];
      
      proportion_secondary_cycle_negative[cycle] += proportion_secondary_negative[center, cycle]*pondef[center,cycle];
      prevalence_secondary_cycle_negative[cycle] += prevalence_secondary_negative[center, cycle]*pondef[center,cycle];
    }
    
    //Max out at 1
    proportion_primary_cycle[cycle] = fmin(proportion_primary_cycle[cycle], 1.0);
    
    proportion_secondary_cycle_positive[cycle] = fmin(proportion_secondary_cycle_positive[cycle], 1.0);
    prevalence_secondary_cycle_positive[cycle] = fmin(prevalence_secondary_cycle_positive[cycle], 1.0);
    
    proportion_secondary_cycle_negative[cycle] = fmin(proportion_secondary_cycle_negative[cycle], 1.0);
    prevalence_secondary_cycle_negative[cycle] = fmin(prevalence_secondary_cycle_negative[cycle], 1.0);
    
    //P(Second) adjusted and unadjusted
    prevalence_secondary_cycle[cycle] = fmin(prevalence_secondary_cycle_positive[cycle]*proportion_primary_cycle[cycle] +
      prevalence_secondary_cycle_negative[cycle]*(1.0 - proportion_primary_cycle[cycle]), 1.0);
      
    proportion_secondary_cycle[cycle] = fmin(proportion_secondary_cycle_positive[cycle]*proportion_primary_cycle[cycle] +
      proportion_secondary_cycle_negative[cycle]*(1.0 - proportion_primary_cycle[cycle]), 1.0);
      
    //P(First|Second) adjusted and unadjusted  
    conditional_proportion_primary_secondary_positive[cycle] = 
      fmin(proportion_secondary_cycle_positive[cycle]*proportion_primary_cycle[cycle] / proportion_secondary_cycle[cycle], 1.0);
      
    conditional_proportion_primary_secondary_negative[cycle] = 
      fmin(proportion_secondary_cycle_negative[cycle]*proportion_primary_cycle[cycle] / proportion_secondary_cycle[cycle], 1.0); 
      
    conditional_prevalence_primary_secondary_positive[cycle] = 
      fmin(proportion_secondary_cycle_positive[cycle]*proportion_primary_cycle[cycle] / proportion_secondary_cycle[cycle], 1.0);
      
    conditional_prevalence_primary_secondary_negative[cycle] = 
      fmin(proportion_secondary_cycle_negative[cycle]*proportion_primary_cycle[cycle] / proportion_secondary_cycle[cycle], 1.0); 
  }
}

model {
  
  sigma1 ~ normal(sigma_hiper, 0.0001);
  sigma2 ~ normal(sigma_hiper, 0.0001);
  sigma3 ~ normal(sigma_hiper, 0.0001);
  sigma4 ~ normal(sigma_hiper, 0.0001);
  sigma5 ~ normal(sigma_hiper, 0.0001);
  
  for (center in 1:Ncenters){
    for (cycle in 1:Ncycles){
      
      //Cases are binomial
      primary_cases[center, cycle]            ~ binomial(ns[center, cycle], prevalence_primary[center, cycle]);
      secondary_positive_cases[center, cycle] ~ binomial(primary_cases[center, cycle], proportion_secondary_positive[center,cycle]);
      secondary_negative_cases[center, cycle] ~ binomial(ns[center, cycle] - primary_cases[center, cycle], proportion_secondary_negative[center,cycle]);
      
      //The adjusted has a beta prior
      prevalence_primary[center,cycle]            ~ beta(alpha_primary_centro_cycle[center, cycle], beta_primary_centro_cycle[center, cycle]);
      prevalence_secondary_positive[center,cycle] ~ beta(alpha_secondary_positive_centro_cycle[center, cycle], beta_secondary_positive_centro_cycle[center, cycle]);
      prevalence_secondary_negative[center,cycle] ~ beta(alpha_secondary_negative_centro_cycle[center, cycle], beta_secondary_negative_centro_cycle[center, cycle]);
      
      //The hyperparameters
      alpha_primary_centro_cycle[center, cycle] ~ normal(alpha_primary_cycle[cycle], sigma1);
      beta_primary_centro_cycle[center, cycle]  ~ normal(beta_primary_cycle[cycle], sigma1);
      
      alpha_secondary_positive_centro_cycle[center, cycle] ~ normal(alpha_secondary_cycle[cycle], sigma2);
      beta_secondary_positive_centro_cycle[center, cycle]  ~ normal(beta_secondary_cycle[cycle], sigma2);
      
      alpha_secondary_negative_centro_cycle[center, cycle] ~ normal(alpha_secondary_cycle[cycle], sigma3);
      beta_secondary_negative_centro_cycle[center, cycle]  ~ normal(beta_secondary_cycle[cycle], sigma3);
      
      if (cycle > 1){
        alpha_primary_cycle[cycle]    ~ normal(alpha_primary_cycle[cycle - 1], sigma4);
        beta_primary_cycle[cycle]     ~ normal(beta_primary_cycle[cycle - 1], sigma4);
        
        alpha_secondary_cycle[cycle]  ~ normal(alpha_secondary_cycle[cycle - 1], sigma5);
        beta_secondary_cycle[cycle]   ~ normal(beta_secondary_cycle[cycle - 1], sigma5);
        
      } else {
        alpha_primary_cycle[cycle]    ~ normal(prior_alpha, sigma4);
        beta_primary_cycle[cycle]     ~ normal(prior_beta, sigma4);
        
        alpha_secondary_cycle[cycle]  ~ normal(prior_alpha, sigma5);
        beta_secondary_cycle[cycle]   ~ normal(prior_beta, sigma5);
      }
    }
  }
}
