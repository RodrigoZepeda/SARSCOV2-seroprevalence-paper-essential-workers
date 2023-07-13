// ANALYSIS FOR BINARY VALUES THAT MIGHT BE ADJUSTED FOR SENSITIVITY/SPECIFICITY
// --------------------------------------------------------------------------------
// Rodrigo Zepeda Tello
// email: rodrigo.zepeda®imss.gob.mx [change ® for appropriate similar symbol]
// 
// Model:
// -------------------------------------------------------------------------------
// TODO edit here
//
// --------------------------------------------------------------------------------

data {
  int<lower=0>  Ncycles;                        //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                       //Número de centers de trabajo analizados
  int<lower=0>  Groups;                         //Número de grupos en la multinomial
  int<lower=0>  cases[Ncycles,Ncenters,Groups]; //Número de individuos por cada grupo por centro y cycle
  real<lower=0> pondef[Ncenters,Ncycles];       //Ponderadores asociados por cycle/centro
  real<lower=0> kappa_prior;
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
}

parameters {
  
  //Probability adjusted and unadjusted
  array[Ncycles, Ncenters] simplex[Groups] theta;
  
  //Hierarchical distribution of ø ~ Dirichlet(k_c*phi(cycle))
  //array[Ncycles, Ncenters] error_c;
  
  real<lower=0> kappa_c[Ncycles];

  //simplex[Groups] phi;
  
  //Variance of kappa
  real<lower=0> sigmasq;
  
}

transformed parameters {
  
  vector[Groups] alpha;
  array[Ncycles] vector[Groups] alpha_c;
  array[Ncycles] vector[Groups] proportions_per_cycle;
  
    //Prior values for alpha  
  alpha = rep_vector(1.0, Groups);


  //From: https://mc-stan.org/docs/2_29/stan-users-guide/reparameterizations.html
  for (cycle in 1:Ncycles){
    
    if (cycle > 1){
      alpha_c[cycle] = kappa_c[cycle] * alpha_c[cycle - 1];
    } else {
      alpha_c[1]     = alpha;
    }
    
    proportions_per_cycle[cycle] = rep_vector(0.0, Groups);
    for (center in 1:Ncenters){
      proportions_per_cycle[cycle] += pondef[center,cycle]*theta[cycle, center];
    }
  }
}

model {
  
  sigmasq  ~ cauchy(sigma_1, sigma_2);
  
  for (cycle in 1:Ncycles){
    
    if (cycle > 1){
      kappa_c[cycle] ~ normal(kappa_c[cycle - 1], sigmasq);
    } else {
      kappa_c[cycle] ~ normal(kappa_prior, sigmasq);
    }
      
    for (center in 1:Ncenters){
      //Cases are multinomial
      theta[cycle, center]  ~ dirichlet(alpha_c[cycle]);
      cases[cycle, center,] ~ multinomial(theta[cycle, center]);
      
    }
  }
}
