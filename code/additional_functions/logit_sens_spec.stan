// ANALYSIS FOR ADJUSTED LOGIT OF IGG 
// --------------------------------------------------------------------------------
// Rodrigo Zepeda Tello
// email: rodrigo.zepeda®imss.gob.mx [change ® for appropriate similar symbol]
// 
// Model:
// -------------------------------------------------------------------------------
//
// w .- Work Center
// c .- Cycle of survey
//
// --------------------------------------------------------------------------------
// NOTE:
// There is a Bernoulli-logit function in Stan; however there is no way to adjust
// for the test's sensitivity and specificity. So we use 
// a specification of the form: bernoulli(inv_logit(alpha + beta * x[n]));

data {
  int<lower=0>  Ncycles;                     //Número de cycles de la encuesta
  int<lower=0>  Ncenters;                    //Número de centers de trabajo analizados
  int<lower=0>  Ncovariates;                 //Número de covariables de la encuesta
  int<lower=0>  N_obs_cycle[Ncenters, Ncycles];     //Número de individuos por ciclo
  int<lower=0>  max_obs_cycle;                      //Máximo número de obs por ciclo (truco)
  real<lower=0> pondef[Ncenters,Ncycles];           //Ponderadores asociados por cycle/centro
  
  //IGG and covariates per cycle
  int<lower=0, upper=1>  igg[Ncenters, Ncycles, max_obs_cycle];  //Respuesta de anticuerpos (1 = sí / 0 = no)
  int<lower=0, upper=1>  cov[Ncenters, Ncycles, Ncovariates, max_obs_cycle];  //Covariables edad, sexo, etc
  
  real<lower=0,upper=1> sens;  //Sensibilidad de la prueba
  real<lower=0,upper=1> spec;  //Especificidad de la prueba
  
  real<lower = 0> sigma_prior;
}

transformed data {
  int<lower=0, upper=Ncovariates> sexo_col    = 1; //Columna con la variable de sexo
  int<lower=0, upper=Ncovariates> edad_1_col  = 2; //Columna con la variable de edad < 35
  int<lower=0, upper=Ncovariates> edad_2_col  = 3; //Columna con la variable de edad 35-44
  int<lower=0, upper=Ncovariates> edad_3_col  = 4; //Columna con la variable de edad 45-60
  int<lower=0, upper=Ncovariates> edad_4_col  = 5; //Columna con la variable de edad 60+
  int<lower=0, upper=Ncovariates> sintoma_col = 6; //Columna con la variable de síntoma
  int<lower=0, upper=Ncovariates> itt_col     = 7; //Columna con la variable de itt
}

parameters {
  real<lower=0> sigma_sigma_hyper;
  real<lower=0> sigma_hyper[Ncovariates];
  real<lower=0> sigma[Ncovariates];
  real<lower=0> sigma_super;
  
  //Centralized betas for regression Normal(0, 1)
  real beta_baseline_center[Ncenters, Ncycles];
  real beta_sexo_center[Ncenters, Ncycles];
  real beta_edad_35_44_center[Ncenters, Ncycles];
  real beta_edad_45_60_center[Ncenters, Ncycles];
  real beta_edad_60_plus_center[Ncenters, Ncycles];
  real beta_sintoma_center[Ncenters, Ncycles];
  real beta_itt_center[Ncenters, Ncycles];
  
  //Hyperparameters per cycle
  real beta_baseline_cycle[Ncycles];
  real beta_sexo_cycle[Ncycles];
  real beta_edad_35_44_cycle[Ncycles];
  real beta_edad_45_60_cycle[Ncycles];
  real beta_edad_60_plus_cycle[Ncycles];
  real beta_sintoma_cycle[Ncycles];
  real beta_itt_cycle[Ncycles];

}

transformed parameters {
  
  //Hyperparameters for report
  real beta_baseline_mean_cycle[Ncycles];
  real beta_sexo_mean_cycle[Ncycles];
  real beta_edad_35_44_mean_cycle[Ncycles];
  real beta_edad_45_60_mean_cycle[Ncycles];
  real beta_edad_60_plus_mean_cycle[Ncycles];
  real beta_sintoma_mean_cycle[Ncycles];
  real beta_itt_mean_cycle[Ncycles];
  
  //Hyperparameters for report
  real average_OR_baseline[Ncycles];
  real average_OR_sexo[Ncycles];
  real average_OR_35[Ncycles];
  real average_OR_45[Ncycles];
  real average_OR_60[Ncycles];
  real average_OR_sintoma[Ncycles];
  real average_OR_itt[Ncycles];
  
  //Betas for regression
  real beta_baseline[Ncenters, Ncycles];
  real beta_sexo[Ncenters, Ncycles];
  real beta_edad_35_44[Ncenters, Ncycles];
  real beta_edad_45_60[Ncenters, Ncycles];
  real beta_edad_60_plus[Ncenters, Ncycles];
  real beta_sintoma[Ncenters, Ncycles];
  real beta_itt[Ncenters, Ncycles];
  
  //Mean of logit
  real mu[Ncenters, Ncycles, max_obs_cycle];
  real prevalence[Ncenters, Ncycles, max_obs_cycle];
  real positive[Ncenters, Ncycles, max_obs_cycle];
  
  for (cycle in 1:Ncycles){
    
      //Start at zero
      beta_baseline_mean_cycle[cycle]     = 0;
      beta_sexo_mean_cycle[cycle]         = 0;
      beta_edad_35_44_mean_cycle[cycle]   = 0;
      beta_edad_45_60_mean_cycle[cycle]   = 0;
      beta_edad_60_plus_mean_cycle[cycle] = 0;
      beta_sintoma_mean_cycle[cycle]      = 0;
      beta_itt_mean_cycle[cycle]          = 0;
      
      average_OR_baseline[cycle]          = 0;
      average_OR_sexo[cycle]              = 0;
      average_OR_35[cycle]                = 0;
      average_OR_45[cycle]                = 0;
      average_OR_60[cycle]                = 0;
      average_OR_sintoma[cycle]           = 0;
      average_OR_itt[cycle]               = 0;
      
    for (center in 1:Ncenters){
      
      beta_baseline[center, cycle]     = beta_baseline_cycle[cycle] + sigma[1]*beta_baseline_center[center, cycle];
      beta_sexo[center, cycle]         = beta_sexo_cycle[cycle] + sigma[2]*beta_sexo_center[center, cycle];;
      beta_edad_35_44[center, cycle]   = beta_edad_35_44_cycle[cycle] + sigma[3]*beta_edad_35_44_center[center, cycle];;
      beta_edad_45_60[center, cycle]   = beta_edad_45_60_cycle[cycle] + sigma[4]*beta_edad_45_60_center[center, cycle];;
      beta_edad_60_plus[center, cycle] = beta_edad_60_plus_cycle[cycle] + sigma[5]*beta_edad_60_plus_center[center, cycle];
      beta_sintoma[center, cycle]      = beta_sintoma_cycle[cycle] + sigma[6]*beta_sintoma_center[center, cycle];
      beta_itt[center, cycle]          = beta_itt_cycle[cycle] + sigma[7]*beta_itt_center[center, cycle];
      
      //Average effect (beta)
      beta_baseline_mean_cycle[cycle]     += beta_baseline[center, cycle]*pondef[center,cycle];
      beta_sexo_mean_cycle[cycle]         += beta_sexo[center, cycle]*pondef[center,cycle];
      beta_edad_35_44_mean_cycle[cycle]   += beta_edad_35_44[center, cycle]*pondef[center,cycle];
      beta_edad_45_60_mean_cycle[cycle]   += beta_edad_45_60[center, cycle]*pondef[center,cycle];
      beta_edad_60_plus_mean_cycle[cycle] += beta_edad_60_plus[center, cycle]*pondef[center,cycle];
      beta_sintoma_mean_cycle[cycle]      += beta_sintoma[center, cycle]*pondef[center,cycle];
      beta_itt_mean_cycle[cycle]          += beta_itt[center, cycle]*pondef[center,cycle];
      
      //Average OR (exp beta)
      average_OR_baseline[cycle]     += exp(beta_baseline[center, cycle])*pondef[center,cycle];
      average_OR_sexo[cycle]         += exp(beta_sexo[center, cycle])*pondef[center,cycle];
      average_OR_35[cycle]           += exp(beta_edad_35_44[center, cycle])*pondef[center,cycle];
      average_OR_45[cycle]           += exp(beta_edad_45_60[center, cycle])*pondef[center,cycle];
      average_OR_60[cycle]           += exp(beta_edad_60_plus[center, cycle])*pondef[center,cycle];
      average_OR_sintoma[cycle]      += exp(beta_sintoma[center, cycle])*pondef[center,cycle];
      average_OR_itt[cycle]          += exp(beta_itt[center, cycle])*pondef[center,cycle];
      
      for (n in 1:N_obs_cycle[center, cycle]){
        mu[center, cycle, n] = beta_baseline[center, cycle] + 
          beta_sexo[center, cycle]*cov[center, cycle, sexo_col, n] +
          beta_edad_35_44[center, cycle]*cov[center, cycle, edad_2_col, n] +
          beta_edad_45_60[center, cycle]*cov[center, cycle, edad_3_col, n] +
          beta_edad_60_plus[center, cycle]*cov[center, cycle, edad_4_col, n] +
          beta_sintoma[center, cycle]*cov[center, cycle, sintoma_col, n] +
          beta_itt[center, cycle]*cov[center, cycle, itt_col, n];
          
        prevalence[center, cycle, n] = inv_logit(mu[center, cycle, n]);
        positive[center, cycle, n]   = fmax(fmin(sens*prevalence[center, cycle, n] + 
              (1 - spec)*(1 - prevalence[center, cycle, n]), 1.0), 0.0);
      } 
    }
  }
  
}

model {
  
  sigma_sigma_hyper ~ normal(0, 2.5);
  sigma_hyper ~ normal(sigma_sigma_hyper, sigma_prior);
  sigma       ~ normal(sigma_super, sigma_prior);
  sigma_super ~ cauchy(0, 2.5);
  
  for (cycle in 1:Ncycles){
    //Dynamic temporal structure
    if (cycle > 1){
      beta_baseline_cycle[cycle]     ~ normal(beta_baseline_cycle[cycle - 1], sigma_hyper[1]); 
      beta_sexo_cycle[cycle]         ~ normal(beta_sexo_cycle[cycle - 1], sigma_hyper[2]); 
      beta_edad_35_44_cycle[cycle]   ~ normal(beta_edad_35_44_cycle[cycle - 1], sigma_hyper[3]); 
      beta_edad_45_60_cycle[cycle]   ~ normal(beta_edad_45_60_cycle[cycle - 1], sigma_hyper[4]); 
      beta_edad_60_plus_cycle[cycle] ~ normal(beta_edad_60_plus_cycle[cycle - 1], sigma_hyper[5]); 
      beta_sintoma_cycle[cycle]      ~ normal(beta_sintoma_cycle[cycle - 1], sigma_hyper[6]); 
      beta_itt_cycle[cycle]          ~ normal(beta_itt_cycle[cycle - 1], sigma_hyper[7]); 
    } else {
      beta_baseline_cycle[cycle]     ~ normal(0, sigma_hyper[1]); 
      beta_sexo_cycle[cycle]         ~ normal(0, sigma_hyper[2]);
      beta_edad_35_44_cycle[cycle]   ~ normal(0, sigma_hyper[3]);
      beta_edad_45_60_cycle[cycle]   ~ normal(0, sigma_hyper[4]);
      beta_edad_60_plus_cycle[cycle] ~ normal(0, sigma_hyper[5]);
      beta_sintoma_cycle[cycle]      ~ normal(0, sigma_hyper[6]);
      beta_itt_cycle[cycle]          ~ normal(0, sigma_hyper[7]);
    }
    
    for (center in 1:Ncenters){
      
      beta_baseline_center[center, cycle]     ~ std_normal();
      beta_sexo_center[center, cycle]         ~ std_normal();
      beta_edad_35_44_center[center, cycle]   ~ std_normal();
      beta_edad_45_60_center[center, cycle]   ~ std_normal();
      beta_edad_60_plus_center[center, cycle] ~ std_normal();
      beta_sintoma_center[center, cycle]      ~ std_normal();
      beta_itt_center[center, cycle]          ~ std_normal();
  
      for (n in 1:N_obs_cycle[center, cycle]){
        igg[center, cycle, n] ~ bernoulli(positive[center, cycle, n]);
      }
    }
  }
}
