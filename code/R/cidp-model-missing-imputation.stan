functions{

  /**********************************************************************************************/
  // function to return the derivatives for one time step
  vector cidp_derivs(vector states, real[] p)
  {
    // parameters
    real a = p[1];
    real b = p[2];
    real e_ = p[3];
    real f = p[4];
    real g = p[5];
    real k = p[6];
    real h = p[7];
    real w = p[8];
    real m = p[9];
    real q = p[10];
    real r = p[11];
    real s = p[12];

    real C = states[1];
    real I = states[2];
    real D = states[3];
    real P = states[4];

    vector[4] derivs;

    derivs[1] = a * C * (P/b - 1) - e_ * C;
    derivs[2] = f*g*C - w*I - I/(D*s + I)*D + k*(h - f*g*C);
    derivs[3] = m*(h*q/P - D);
    derivs[4] = r*P*( s*D/I - 1);

    return derivs[];

  }//derivs

  /**********************************************************************************************/
  // Runge-Kutta 4th order numerical integration scheme
  vector integrate_rk4(vector states, real[] p, real dt){

    vector[dims(states)[1]] K[4];
    real int_constant = 1.0/6.0;

    K[1,] = cidp_derivs(states, p) * dt;
    K[2,] = cidp_derivs(states + K[1,]*0.5, p) * dt;
    K[3,] = cidp_derivs(states + K[2,]*0.5, p) * dt;
    K[4,] = cidp_derivs(states + K[3,], p) * dt;

    return (states + (K[1,] + (K[2,] + K[3,]) * 2 + K[4,]) * int_constant);
  }
  /**********************************************************************************************/
  // numerically-intergrate the capital-inventory-demand-price ODE model
  vector[] cidp_ode(int n_sims, int n_states, real[] sv_inits, real[] p, real dt){

    vector[n_states] states[n_sims];

    // initial state variable values
    states[1, 1] = sv_inits[1];
    states[1, 2] = sv_inits[2];
    states[1, 3] = sv_inits[3];
    states[1, 4] = sv_inits[4];

    // loop through the simulation
    for(i in 2:n_sims) states[i,:] = integrate_rk4(states[i-1,:], p, dt);

    // need to return every 1/dt^th value in states

    return states[];
  }//ode

  // use the in-built integrate_ode_rk45 solver in Stan instead
  real[] cidp_ode_rk45(real t, real[] y, real[] theta,  real[] x_r, int[] x_i)
  {
    // parameters
    real a = theta[1];
    real b = theta[2];
    real e_ = theta[3];
    real f = theta[4];
    real g = theta[5];
    real k = theta[6];
    real h = theta[7];
    real w = theta[8];
    real m = theta[9];
    real q = theta[10];
    real r = theta[11];
    real s = theta[12];

    real C = y[1];
    real I = y[2];
    real D = y[3];
    real P = y[4];

    real derivs[4];

    derivs[1] = a * C * (P/b - 1) - e_ * C;
    derivs[2] = f*g*C - w*I - I/(D*s + I)*D + k*(h - f*g*C);
    derivs[3] = m*(h*q/P - D);
    derivs[4] = r*P*( s*D/I - 1);

    return derivs[];
  }

}

data{
  int<lower=1> n_t;
  real ts[n_t];
  int n_sims;
  int sims_seq[n_t];
  int<lower=1> n_states;
  int<lower=0> n_parameters;
  real<upper=1> dt;
  //real<lower=0> initial_values[n_states];
  real y[n_t, n_states];
  real<lower=0> domestic_production[n_t];
  real<lower=0> consumption[n_t];
  real b_fix;
  real e_fix;
  real f_fix;
  real g_fix;
  real h_fix;
  int<lower=0> n_missing_C;
  int<lower=0> n_missing_I;
  int<lower=0> n_missing_D;
  int<lower=0> n_missing_P;
  int<lower=0> missing_idx_C[n_missing_C];
  int<lower=0> missing_idx_I[n_missing_I];
  int<lower=0> missing_idx_D[n_missing_D];
  int<lower=0> missing_idx_P[n_missing_P];
  int<lower=0,upper=1> prior_only;
  int<lower=0,upper=1> rk4;
  //vector[n_parameters] p;

}

transformed data{
  // these are fixed parameters -- parameters that do not need to be estimated
  real x_r[0];
  int x_i[0];
}

parameters{
  // parameters on the transformed scale, easier to estimate
  real<lower=0> a;
  //real<lower=0> b;
  //real<lower=0> e_;
  //real<lower=0> f;
  //real<lower=0> g;
  real<lower=0,upper=1> k;
  //real<lower=0> h;
  real<lower=0> w;
  real<lower=0> m;
  real<lower=0> q;
  real<lower=0> r;
  real<lower=0> s;
  real<lower=0> states_init[4];
  vector<lower=0>[4] sigma;
  real<lower=0> sigma_domestic_production;
  real<lower=0> sigma_consumption;
  // missing data parameters
  real<lower=0> C_impute[n_missing_C];
  real<lower=0> I_impute[n_missing_I];
  real<lower=0> D_impute[n_missing_D];
  real<lower=0> P_impute[n_missing_P];
}

transformed parameters{
  //vector[n_states] states[n_sims];
  real states[n_t, n_states];
  real p[n_parameters];
  real initial_values[n_states];

  // model parameters -- need to re-scale from the prior values to aid the HMC sampler
  p[1] = a / 100;
  // divide by 100
  p[2] = b_fix; //b * 100;
  p[3] = e_fix; //e_ / 100;
  p[4] = f_fix; // f / 100;
  // divided by 1000
  p[5] = g_fix; //g * 1000;
  p[6] = k;
  // divided by 10e6
  p[7] = h_fix; //h * 10e6;
  p[8] = w;
  p[9] = m;
  // divided by 100
  p[10] = q * 1000;
  p[11] = r / 10;
  p[12] = s;

  initial_values[1] = states_init[1] * 1e6;
  initial_values[2] = states_init[2] * 1e8;
  initial_values[3] = states_init[3] * 1e8;
  initial_values[4] = states_init[4] * 1e3;

  // intergate
  //if(rk4) states = cidp_ode(n_sims, n_states, initial_values, p, dt);
  if(!rk4) states = integrate_ode_rk45(cidp_ode_rk45, initial_values, 0, ts, p, x_r, x_i);
}

model{
  real y_merged[n_t, n_states]; // merged missing and observed y variables

  // build y_merged
  y_merged = y;
  for(j in 1:n_missing_C) y_merged[ missing_idx_C[j], 1] = C_impute[j];
  for(j in 1:n_missing_I) y_merged[ missing_idx_I[j], 2] = I_impute[j];
  for(j in 1:n_missing_D) y_merged[ missing_idx_D[j], 3] = D_impute[j];
  for(j in 1:n_missing_P) y_merged[ missing_idx_P[j], 4] = P_impute[j];


  // priors, all on a reasonable scale
  a ~ normal(0, 1);
  //b ~ normal(0, 1);
  //e_ ~ normal(0, 1);
  //f ~ normal(0, 1);
  //g ~ normal(1.98, 0.1);
  k ~ beta(2, 2);
  //h ~ normal(3, 1);
  w ~ normal(0, 1);
  m ~ normal(0, 1);
  q ~ normal(0, 1);
  r ~ normal(0, 1);
  s ~ normal(0, 1);

  states_init[1] ~ lognormal( log(y[1,1]/1e6), 1);
  states_init[2] ~ lognormal( log(y[1,2]/1e8), 1);
  states_init[3] ~ lognormal( log(h_fix/1e8), 1);
  states_init[4] ~ lognormal( log(y[1,4]/1e3), 1);

  sigma ~ lognormal(-1, 1);
  sigma_domestic_production ~ lognormal(-1, 1);
  sigma_consumption ~ lognormal(-1, 1);

  // likelihood
  if(!prior_only){

    // rk4 likelihood
    if(rk4){
      domestic_production ~ lognormal( log(to_vector(states[sims_seq, 1])*f_fix*g_fix), sigma_domestic_production);
      for(t in 1:n_t) consumption[t] ~ lognormal( log(states[sims_seq[t], 3] * states[sims_seq[t], 2]/(states[sims_seq[t], 3]*s + states[sims_seq[t], 2])), sigma_consumption);
      for(n in 1:n_states) y_merged[,n] ~ lognormal(log(states[sims_seq,n]), sigma[n]);
    }
    else{
      // rk45 likelihood
      domestic_production ~ lognormal( log(to_vector(states[, 1])*f_fix*g_fix), sigma_domestic_production);
      for(t in 1:n_t) consumption[t] ~ lognormal( log(states[t, 3] * states[t, 2]/(states[t, 3]*s + states[t, 2])), sigma_consumption);
      for(n in 1:n_states) y_merged[,n] ~ lognormal(log(states[,n]), sigma[n]);
    }
  }
}

generated quantities{
  // transform to dimensionless quantities
  real alpha = p[10]/b_fix;
  real omega = p[8]/p[1];
  real gamma = 1/(p[1] * p[12]);
  real kappa = p[6];
  real beta = e_fix/p[1];
  // critical ratio
  real critical_ratio = alpha*(omega + gamma/2) / (kappa*gamma*(1 + beta));
}
