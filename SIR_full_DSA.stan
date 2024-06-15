
functions{
  real[] SVRmodel(real t, real[] y, real[] parms, real[] x_r, int[] x_i){
    real alpha = parms[1];
    real beta = parms[2];
    real gamma = parms[3];
    real rho = parms[4];
    
    real dydt[8];
    dydt[1] = -beta*y[1]*y[3]; // S
    dydt[2] = beta*y[1]*y[3] - alpha*y[2]; // E
    dydt[3] = alpha*y[2] - gamma*y[3]; // I
    dydt[4] = gamma*y[3]; // R
    dydt[5] = beta*y[1]*y[3] - gamma*y[3]; // V
    
    dydt[6] = -beta*y[3]*y[6]*y[7]/y[5]; // \tilde{S}
    dydt[7] = beta*y[3]*y[6]*y[7]/y[5] - gamma*y[3]*y[7]/y[5]; // \tilde{I}
    dydt[8] = gamma*y[3]*y[7]/y[5];
    return dydt;
  }
  real[] recovery_density(real t, real[] y, real[] parms, real[] x_r, int[] x_i){
    real alpha = parms[1];
    real beta = parms[2];
    real gamma = parms[3];
    real rho = parms[4];
    
    real dydt[14];
    dydt[1] = -beta*y[1]*y[3]; // S
    dydt[2] = beta*y[1]*y[3] - alpha*y[2]; // E
    dydt[3] = alpha*y[2] - gamma*y[3]; // I
    dydt[4] = gamma*y[3]; // R
    dydt[5] = beta*y[1]*y[3] - gamma*y[3]; // V
    
    dydt[6] = -beta*y[3]*y[6]*y[7]/y[5]; // \tilde{S}
    dydt[7] = beta*y[3]*y[6]*y[7]/y[5] - gamma*y[3]*y[7]/y[5]; // \tilde{I}
    dydt[8] = gamma*y[3]*y[7]/y[5]; //\tilde{R}
    
    dydt[9] = gamma*(dydt[3]*y[5] - y[3]*dydt[5])/(y[5]*y[5]); //gamma_t
    dydt[10] = y[9]; // integral of gamma_t
    dydt[11] = dydt[9]*exp(-y[10]) - y[9]*y[11]; //r_1
    dydt[12] = -dydt[6]*exp(y[10]); //r_2
    dydt[13] = dydt[11]*y[12] + y[11]*dydt[12]; //g
    dydt[14] = y[13]; //G
    return dydt;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; 
  real<lower=0.0> infection_times[N];
  int<lower=0> K;
  real<lower=0.0> I2R_times[K];
  //real<lower=0.0> alpha_0;
  //real<lower=0.0> beta_0;
  //real<lower=0.0> gamma_0;
  //real<lower=0.0> rho_0;
}

transformed data {
    real x_r[0];
    int x_i[0];
}

parameters {
    real<lower=0.0> alpha;
    real<lower=0.0, upper=10.0> beta;
    //real<lower=0.0> R0;
    real<lower=0.0, upper=beta> gamma;
    real<lower=0.0, upper=1.0> rho;
}

transformed parameters{
  real<lower=0.0> R0 = beta/gamma;
  //real<lower=0.0> beta = R0*gamma;
}

model {
  real parms[4];
  real ic[8];
  real s[N,8];
  real r_ic[14];
  real r[K, 14];
  
  real t0;
  real Smax;
  real factor;

  parms[1] = alpha;
  parms[2] = beta;
  parms[3] = gamma;
  parms[4] = rho;
  
  ic[1] = 1.0; 
  ic[2] = 0.0;
  ic[3] = rho;
  ic[4] = 0.0;
  ic[5] = rho; 
  ic[6] = 1.0;
  ic[7] = rho;
  ic[8] = 0.0;
  
  t0 = 0.0; 
  
  s = integrate_ode_rk45(SVRmodel,ic,t0,infection_times,parms,x_r,x_i);
  Smax = s[N,6];
  factor = 1.0 - Smax;
  
  for (i in 1:N){
        target += log((beta*s[i,6]*s[i,3]*s[i,7]/s[i,5])/factor);
    }
    
  r_ic[1] = 1.0; 
  r_ic[2] = 0.0;
  r_ic[3] = rho;
  r_ic[4] = 0.0;
  r_ic[5] = rho; 
  r_ic[6] = 1.0;
  r_ic[7] = rho;
  r_ic[8] = 0.0;
  r_ic[9] = gamma;
  r_ic[10] = 0.0;
  r_ic[11] = gamma;
  r_ic[12] = beta*rho;
  r_ic[13] = beta*gamma*rho;
  r_ic[14] = 0.0;
  
  r = integrate_ode_rk45(recovery_density,r_ic, t0, I2R_times, parms, x_r,x_i);
  Smax = r[K, 14];
  
  for (i in 1:K){
    target += log(r[i, 13]/Smax);
  }
  
  //target += exponential_lpdf(alpha | alpha_0);
  //target += exponential_lpdf(beta | beta_0) ;
  //target += exponential_lpdf(gamma | gamma_0);
  //target += exponential_lpdf(rho | rho_0);
  
  
  
  
}

