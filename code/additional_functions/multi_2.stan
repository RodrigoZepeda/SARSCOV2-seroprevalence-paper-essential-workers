// ANALYSIS FOR BINARY VALUES THAT MIGHT BE ADJUSTED FOR SENSITIVITY/SPECIFICITY
// --------------------------------------------------------------------------------
// Rodrigo Zepeda Tello
// email: rodrigo.zepeda®imss.gob.mx [change ® for appropriate similar symbol]
// 
// Model:
// -------------------------------------------------------------------------------
// TODO edit here
// TODO multiply variance by timespan
// --------------------------------------------------------------------------------

data {
  int<lower=0>  Ncycles;                        //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                       //Número de centers de trabajo analizados
  int<lower=0>  Groups;                         //Número de grupos en la multinomial
  int<lower=0>  cases[Ncycles,Ncenters,Groups]; //Número de individuos por grupo por center y cycle
  int<lower=0>  cases_in_group[Ncycles,Ncenters,Groups]; //Número de + por gpo center y cycle
  real<lower=0> pondef[Ncenters,Ncycles];       //Ponderadores asociados por cycle/center
  real<lower=0,upper=1> sens;                   //Sensibilidad de la prueba
  real<lower=0,upper=1> spec;                   //Especificidad de la prueba
  real<lower=0> prior_alpha;                    //A priori para alpha
  real<lower=0> prior_beta;                     //A priori para beta
  real<lower=0> sigma_hiper;
  real<lower=0> kappa_prior;
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
}

parameters {
  
  //Probability adjusted and unadjusted
  array[Ncycles, Ncenters] simplex[Groups] theta;
  array[Ncycles, Ncenters, Groups] real<lower=0,upper=1> p_adjusted;
  
  //Hierarchical distribution of ø ~ Dirichlet(k_c*phi(cycle))
  //real error_c[Ncycles];
  //array[Ncycles] simplex[Groups] phi_c;
  //real<lower=0> kappa;
  real<lower=0> kappa_c[Ncycles];
  
  //Variance of kappa
  real<lower=0,upper=1> p;
  real<lower=0> sigmasq;
  real<lower=0> sigmasq2;
  real<lower=0> sigmasq3;
  real<lower=0> sigmasq4;

  //Hierarchical distribution of p
  real<lower=0> alpha_center_cycle_group[Ncycles, Ncenters, Groups];
  real<lower=0> beta_center_cycle_group[Ncycles, Ncenters, Groups];
  real<lower=0> alpha_cycle_center[Ncycles, Ncenters];
  real<lower=0> beta_cycle_center[Ncycles, Ncenters];
  real<lower=0> alpha_cycle[Ncycles];
  real<lower=0> beta_cycle[Ncycles];
  
}

transformed parameters {
  
  vector[Groups] alpha;
  array[Ncycles] vector[Groups] alpha_c;
  array[Ncycles] vector[Groups] proportions_per_cycle;
  array[Ncycles, Ncenters, Groups] real<lower=0,upper=1> p_unadjusted;
  matrix[Ncycles, Groups] p_unadjusted_cycle;
  matrix[Ncycles, Groups] p_adjusted_cycle;
  
  //From: https://mc-stan.org/docs/2_29/stan-users-guide/reparameterizations.html
  p_unadjusted_cycle           = rep_matrix(0.0, Ncycles, Groups);
  p_adjusted_cycle             = rep_matrix(0.0, Ncycles, Groups);
    
  //Prior values for alpha  
  alpha = rep_vector(1.0, Groups);
  
  for (cycle in 1:Ncycles){
    
    if (cycle > 1){
      alpha_c[cycle] = kappa_c[cycle] * alpha_c[cycle - 1];
    } else {
      alpha_c[1] = alpha;
    }
    
    proportions_per_cycle[cycle] = rep_vector(0.0, Groups);
    for (center in 1:Ncenters){
      proportions_per_cycle[cycle] += pondef[center,cycle]*theta[cycle, center];
    }
  }
  
  for (group in 1:Groups){
    for (cycle in 1:Ncycles){
      for (center in 1:Ncenters){
      
        p_unadjusted[cycle,center,group] = fmin(sens*p_adjusted[cycle,center,group] + 
          (1 - spec)*(1 - p_adjusted[cycle,center,group]), 1.0);
          
        p_adjusted_cycle[cycle,group]   += pondef[center,cycle]*p_adjusted[cycle, center, group];
        p_unadjusted_cycle[cycle,group] += pondef[center,cycle]*p_unadjusted[cycle, center, group];
        
      }
    }
  }
}

model {
  sigmasq   ~ cauchy(sigma_1, sigma_2);
  sigmasq2  ~ normal(sigma_hiper, 0.0001);
  sigmasq3  ~ normal(sigma_hiper, 0.0001);
  sigmasq4  ~ normal(sigma_hiper, 0.0001);
  
  for (cycle in 1:Ncycles){
    
   if (cycle > 1){
      kappa_c[cycle] ~ normal(kappa_c[cycle - 1], sigmasq);
    } else {
      kappa_c[cycle] ~ normal(kappa_prior, sigmasq);
    }
    
    if (cycle > 1){
      
        alpha_cycle[cycle] ~ normal(alpha_cycle[cycle - 1], sigmasq3);
        beta_cycle[cycle]  ~ normal(beta_cycle[cycle - 1], sigmasq3);

      } else {
        
        alpha_cycle[cycle] ~ normal(prior_alpha, sigmasq3);
        beta_cycle[cycle]  ~ normal(prior_beta, sigmasq3);
    }
      
    for (center in 1:Ncenters){
      
      alpha_cycle_center[cycle, center] ~ normal(alpha_cycle[cycle], sigmasq4);
      beta_cycle_center[cycle, center]  ~ normal(beta_cycle[cycle], sigmasq4);
      
      //Cases are multinomial
      theta[cycle, center]  ~ dirichlet(alpha_c[cycle]);
      cases[cycle, center,] ~ multinomial(theta[cycle, center]);
      
      for (group in 1:Groups){
        
        //The hyperparameters
        alpha_center_cycle_group[cycle, center, group] ~ normal(alpha_cycle_center[cycle, center], sigmasq2);
        beta_center_cycle_group[cycle, center, group]  ~ normal(beta_cycle_center[cycle, center], sigmasq2);
        
        p_adjusted[cycle, center, group]   ~ beta(alpha_center_cycle_group[cycle, center, group], 
                                                  beta_center_cycle_group[cycle, center, group]);
        cases_in_group[cycle,center,group] ~ binomial(cases[cycle, center, group],  
                                                    p_unadjusted[cycle, center, group]);
        
      }
    }
  }
}
