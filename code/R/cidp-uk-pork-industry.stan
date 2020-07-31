functions{

  // use the in-built integrate_ode_rk45 solver in Stan instead
  real[] cidp_ode(real t, real[] y, real[] theta,  real[] x_r, int[] x_i)
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
  int<lower=1> n_states;
  int<lower=0> n_parameters;
  real y[n_t, n_states];
  real slaughtered[n_t];
  real consumption[n_t];
  int<lower=0> n_missing_C;
  int<lower=0> n_missing_I;
  int<lower=0> n_missing_D;
  int<lower=0> n_missing_P;
  int<lower=0> n_missing_slaughtered;
  int<lower=0> n_missing_consumption;
  int<lower=0> missing_idx_C[n_missing_C];
  int<lower=0> missing_idx_I[n_missing_I];
  int<lower=0> missing_idx_D[n_missing_D];
  int<lower=0> missing_idx_P[n_missing_P];
  int<lower=0> missing_idx_slaughtered;
  int<lower=0> missing_idx_consumption[n_missing_consumption];
  int<lower=0,upper=1> prior_only;
}

transformed data{
  // these are fixed parameters -- parameters that do not need to be estimated
  real b_fix = 140;
  real e_fix = 1/2.5;
  real f_fix = 24;
  real g_fix = 110/1000.0 * 0.75;
  real h_fix = 800;
  real w_fix = 0.3;
  real s_fix = 1;
  real x_r[0];
  int x_i[0];
}

parameters{
  // parameters on the transformed scale, easier to estimate
  real<lower=0> a;
  real<lower=0,upper=1> k;
  //real<lower=0> w;
  real<lower=0> m;
  real<lower=0> q;
  real<lower=0> r;
  //real<lower=0> s;
  real<lower=0> states_init[4];
  vector<lower=0>[4] sigma;
  real<lower=0> sigma_domestic_production;
  real<lower=0> sigma_consumption;
  // missing data parameters
  real<lower=0> C_impute[n_missing_C];
  real<lower=0> I_impute[n_missing_I];
  real<lower=0> D_impute[n_missing_D];
  real<lower=0> P_impute[n_missing_P];
  real<lower=0> slaughtered_impute;
  real<lower=0> consumption_impute[n_missing_consumption];
}

transformed parameters{
  real states[n_t, n_states];
  real p[n_parameters];
  real initial_values[n_states];

  // model parameters -- need to re-scale from the prior values to aid the HMC sampler
  p[1] = a;
  p[2] = b_fix;
  p[3] = e_fix;
  p[4] = f_fix;
  p[5] = g_fix;
  p[6] = k;
  p[7] = h_fix;
  p[8] = w_fix;
  p[9] = m;
  // q = reference price per kg / 100
  p[10] = q * 100;
  p[11] = r;
  p[12] = s_fix;

  initial_values[1] = states_init[1] * 100;
  initial_values[2] = states_init[2] * 1000;
  initial_values[3] = states_init[3] * 1000;
  initial_values[4] = states_init[4] * 100;

  // intergate
  states = integrate_ode_rk45(cidp_ode, initial_values, 0, ts, p, x_r, x_i, 1e-7, 1e-7, 1e5);
}

model{
  real y_merged[n_t, n_states] = y; // merged missing and observed state variables
  real slaughtered_merged[n_t] = slaughtered;
  real consumption_merged[n_t] = consumption;

  // build merged vectors
  for(j in 1:n_missing_C) y_merged[ missing_idx_C[j], 1] = C_impute[j];
  for(j in 1:n_missing_I) y_merged[ missing_idx_I[j], 2] = I_impute[j];
  for(j in 1:n_missing_D) y_merged[ missing_idx_D[j], 3] = D_impute[j];
  for(j in 1:n_missing_P) y_merged[ missing_idx_P[j], 4] = P_impute[j];
  slaughtered_merged[ missing_idx_slaughtered ] = slaughtered_impute;
  for(j in 1:n_missing_consumption) consumption_merged[ missing_idx_consumption[j] ] = consumption_impute[j];

  // priors, all on a reasonable scale
  a ~ normal(0, 1);
  k ~ beta(2, 2);
  //w ~ normal(0, 1);
  m ~ normal(0, 1);
  q ~ normal(0, 1);
  r ~ normal(0, 1);
  //s ~ normal(0, 1);

  states_init[1] ~ lognormal( log(800/100.0), 1);
  states_init[2] ~ lognormal( log(1500/1000.0), 1);
  states_init[3] ~ lognormal( log(1500/1000.0), 1);
  states_init[4] ~ lognormal( log(150/100.0), 1);

  sigma ~ lognormal(-1, 1);
  sigma_domestic_production ~ lognormal(-1, 1);
  sigma_consumption ~ lognormal(-1, 1);

  // likelihood
  if(!prior_only){
      slaughtered_merged ~ lognormal( log( to_vector(states[,1]) * f_fix), sigma_domestic_production );
      for(t in 1:n_t) consumption_merged[t] ~ lognormal( log(states[t, 3] * states[t, 2]/(states[t, 3]*s_fix + states[t, 2])), sigma_consumption);
      for(n in 1:n_states) y_merged[,n] ~ lognormal(log(states[,n]), sigma[n]);
  }
}

generated quantities{
  // transform to dimensionless quantities
  real alpha = p[10]/p[2];
  real omega = p[8]/p[1];
  real gamma = 1/(p[1] * p[12]);
  real kappa = p[6];
  real beta = p[3]/p[1];
  // critical ratio
  real critical_ratio = alpha*(omega + gamma/2) / (kappa*gamma*(1 + beta));
}
