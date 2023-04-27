functions{

  // ODE model definition to be passed to integrate_ode_rk45 solver in Stan
  vector cidp_ode(real t, vector y, real[] theta)
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

    vector[4] derivs;

    derivs[1] = a * C * (P/b - 1) - e_ * C;
    derivs[2] = f*g*C - w*I - I/(D*s + I)*D + k*(h - f*g*C);
    derivs[3] = m*(h*q/P - D);
    derivs[4] = r*P*( s*D/I - 1);

    return derivs;
  }

}

data{
  int<lower=1> n_t;
  real<lower=0> ts[n_t];
  int<lower=1> n_states;
  int<lower=0> n_parameters;
  matrix[n_t, n_states] y;
  vector[n_t] uk_production;
  vector[n_t] imports;
  vector[n_t] exports;
  int<lower=0,upper=1> prior_only;
}

transformed data{
  // these are fixed parameters -- parameters that do not need to be estimated
  real b_fix = 130;
  real g_fix = 110 * 0.75;
}

parameters{
  // parameters on the transformed scale, easier to estimate
  real<lower=0> a;
  real<lower=0> e_;
  real<lower=0> f;
  real<lower=0,upper=1> k;
  real<lower=0> w;
  real<lower=0> h;
  real<lower=0> m;
  real<lower=0> q;
  real<lower=0> r;
  real<lower=0> s;
  vector<lower=0>[n_states] states_init;
  vector<lower=0>[n_states] sigma;
  real<lower=0> sigma_production;
  real<lower=0> sigma_trade[2];
}

transformed parameters{
  vector[n_states] states[n_t];
  real p[n_parameters];
  vector[n_states] initial_values;

  // model parameters -- need to re-scale from the prior values to aid the HMC sampler
  p[1] = a / 10.0;
  p[2] = b_fix;
  p[3] = e_ / 10.0;
  p[4] = f;
  p[5] = g_fix;
  p[6] = k;
  p[7] = h * 1e8;
  p[8] = w;
  p[9] = m;
  p[10] = q * 1000;
  p[11] = r;
  p[12] = s;

  initial_values[1] = states_init[1] * 1e5;
  initial_values[2] = states_init[2] * 1e8;
  initial_values[3] = states_init[3] * 1e8;
  initial_values[4] = states_init[4] * 1000;

  // intergate
  states = ode_rk45(cidp_ode, initial_values, 0, ts, p);
}

model{

  // priors, all on a reasonable scale
  a ~ normal(0, 1);
  e_ ~ normal(0, 1);
  f ~ normal(0, 1);
  k ~ beta(2, 2);
  h ~ normal(0, 1);
  w ~ normal(0, 1);
  m ~ normal(0, 1);
  q ~ normal(0, 1);
  r ~ normal(0, 1);
  s ~ normal(0, 1);

  states_init[1] ~ lognormal( log(y[1,1]/1e5), 0.5);
  states_init[2] ~ lognormal( log(y[1,2]/1e8), 0.5);
  states_init[3] ~ lognormal( log(y[1,3]/1e8), 0.5);
  states_init[4] ~ lognormal( log(y[1,4]/1e3), 0.5);

  sigma ~ lognormal(-1, 1);
  sigma_production ~ lognormal(-1, 1);
  sigma_trade ~ lognormal(-1, 1);

  // likelihood
  if(!prior_only){
    for(n in 1:n_states){
      for(i in 1:n_t){
        if(y[i,n] != -100)
          y[i,n] ~ lognormal(log(states[i,n]), sigma[n]);
      }//i
    }//n

    uk_production ~ lognormal( log( to_vector(states[,1]) * p[4] * p[5] ), sigma_production);
    imports ~ lognormal( log( p[6] * p[7]), sigma_trade[1]);
    exports ~ lognormal( log( p[6] * to_vector(states[,1]) * p[4] * p[5] ), sigma_trade[2]);
  }//Lk
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
  // surplus ratio
  real surplus_ratio = alpha*(gamma + 2*omega) / (2*gamma*(1 + beta));
}
