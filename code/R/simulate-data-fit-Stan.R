# Simulate data from the model and fit in Stan

source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/ode-model.R")
source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/linear-stability-analysis-functions.R")

# set the parameters
a <- 1/52
b <- 80
e <- 1/(2.5*52)
f <- 1/52
g <- 110*24*0.75
k <- 0.9
h <- 30e6
w <- 0.3
m <- 1/10
q <- 150
r <- 1/20
s <- 1
C_init <- 100e3
I_init <- 30e6
D_init <- 30e6
P_init <- 140

states <- c(C = C_init, I = I_init, D = D_init, P = P_init)

dt <- 0.01
n_t <- 60
t <- seq(0, n_t, dt)

run_cfs <- as.data.frame(
  deSolve::ode(
    func = m_cfs, 
    y = states, 
    p = c(),
    times = t,
    method = "rk4"
  )
)

par(mfrow=c(2,2))
plot(run_cfs$time, run_cfs$C, type="l", bty="l")
plot(run_cfs$I, type="l", bty="l")
plot(run_cfs$D, type = "l", bty="l")
plot(run_cfs$P, type = "l", bty="l")
par(mfrow=c(1,1))

# select time points (i.e. data is not instantaneous)
real_t <- seq(1, max(t), 1)
d <- run_cfs[ run_cfs$time %in% real_t, ]

# add some noise to the observations
#set.seed(2020)
d$C_noise <- rlnorm( nrow(d), log(d$C), 0.1) #d$C + rnorm(nrow(d), 0, 1000)
d$I_noise <- rlnorm( nrow(d), log(d$I), 0.1) #d$I + rnorm(nrow(d), 0, 1e5)
d$D_noise <- rlnorm( nrow(d), log(d$D), 0.1) #d$D + rnorm(nrow(d), 0, 1e5)
d$P_noise <- rlnorm( nrow(d), log(d$P), 0.1) #d$P + rnorm(nrow(d), 0, 5)
d$domestic_production_noise <- rlnorm( nrow(d), log(d$C * f * g), 0.1)
d$consumption_noise <- rlnorm( nrow(d), log(d$I/(d$D*s + d$I)*d$D), 0.1)

#-------------------------------------------------------------------------------
# fit the model in Stan

p_list <- c(
  a = a,
  b = b,
  e_ = e,
  f = f,
  g = g,
  k = k,
  h = h,
  w = w,
  m = m,
  q = q,
  r = r,
  s = s
)

#initial_values <- c(C = 100e3, I = 30e6, D = 30e6, P = 140)

n_t <- nrow(d)
dt <- 0.5;
n_sims <- n_t * dt^-1 + 1;
sims_seq <- seq(1,n_sims,1/dt)[-1]
n_states <- 4
n_parameters <- 12

# compile the model
cipd_model <- rstan::stan_model("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/cidp-model.stan")

# fit the model
fit_model <- rstan::sampling( 
  object = cipd_model,
  data = list(
    n_t = n_t, 
    ts = 1:n_t,
    n_sims = n_sims, 
    dt = dt,
    sims_seq = sims_seq,
    n_states = n_states, 
    n_parameters = n_parameters,
    y = d[,grep("noise",colnames(d))][,1:n_states],
    domestic_production = d$domestic_production_noise,
    consumption = d$consumption_noise,
    #initial_values = initial_values,
    prior_only=0
    #p = p_list
    ), 
  init = 0, 
  cores = 4, chains = 4, warmup = 200, iter = 400
  #,control = list(adapt_delta=0.9, max_treedepth=15)
)

print(fit_model, probs = c(0.025, 0.975), digits_summary = 2, pars="p")
print(fit_model, probs = c(0.025, 0.975), digits_summary = 2, 
      pars=c("a","k","w","m","q","r","s")
      )
pairs(fit_model, pars=c("a","k","w","m","q","r","s","lp__"))

draws <- as.data.frame(fit_model)

n_draws <- nrow(draws)

# parameter estimates
par(mfrow=c(3, 4))
for(i in 1:length(p_list)){
  plot(density(draws[,grep("\\bp\\b", colnames(draws))][,i]), main=paste0("p", i), xlab="", ylab="", lwd=3)
  abline(v=p_list[i], lty=2, lwd=2, col="blue")
}
par(mfrow=c(1,1))

# posterior predictions
n_samples <- 100
sample_draws <- sample(1:n_draws, n_samples, replace=F)

states_init <- draws[ , grep("initial", colnames(draws))]
states_draws <- draws[ , grep("states", colnames(draws))][,-c(1:4)]
par(mfrow=c(2,2))

for(n in 1:n_states){
  int_states <- seq(1,n_sims,1/dt)
  #the_states <- states_draws[, ( ((n-1)*(n_sims)+1):(n_sims*n))][,int_states]
  the_states <- cbind(states_init[,n], states_draws[ , ( ((n-1)*(n_t)+1):(n_t*n))])
  
  inital <- run_cfs[1,1+n]
  
  plot(1:(n_t+1), c(inital, d[,5+n]), type="n" )
  points(1:(n_t+1), c(inital, d[,5+n]), col=scales::alpha("slateblue",0.6), pch=16)
  points(1, inital, col="slateblue", pch=16)
  for(i in 1:n_samples){
    lines(1:(n_t+1), the_states[sample_draws[i], ], 
          col = scales::alpha("black", 0.2), lwd=0.5)
  }
}
par(mfrow=c(1,1))

#--------------------------------------------------------------------------------------------------
# estimate with missing data

d_missing <- d

# add missing data to the capital variable
p_missing_capital <- 0.0
sample_missing_capital <- sample(2:nrow(d_missing), round(p_missing_capital*nrow(d_missing)), replace=F)
d_missing$C_noise[sample_missing_capital] <- -100
n_missing_C <- length(sample_missing_capital)
missing_idx_C <- sort(sample_missing_capital)

# add missing data to the inventory variable
p_missing_inventory <- 0
sample_missing_inventory <- sample(2:nrow(d_missing), round(p_missing_inventory*nrow(d_missing)), replace=F)
d_missing$I_noise[sample_missing_inventory] <- -100
n_missing_I <- length(sample_missing_inventory)
missing_idx_I <- sort(sample_missing_inventory)

# add missing data to the demand variable -- likely not to have any data on demand
p_missing_demand <- 0.9
sample_missing_demand <- sample(2:nrow(d_missing), round(p_missing_demand*nrow(d_missing)), replace=F)
d_missing$D_noise[sample_missing_demand] <- -100
n_missing_D <- length(sample_missing_demand)
missing_idx_D <- sort(sample_missing_demand)

# add missing data to the price variable
p_missing_price <- 0
sample_missing_price <- sample(2:nrow(d_missing), round(p_missing_price*nrow(d_missing)), replace=F)
d_missing$P_noise[sample_missing_price] <- -100
n_missing_P <- length(sample_missing_price)
missing_idx_P <- sort(sample_missing_price)

n_t <- nrow(d)
dt <- 0.25;
n_sims <- n_t * dt^-1 + 1;
sims_seq <- seq(1,n_sims,1/dt)[-1]
n_states <- 4
n_parameters <- 12

# compile the model
cipd_model_missing <- rstan::stan_model("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/cidp-model-missing-imputation.stan")

# fit the model
fit_model_missing <- rstan::sampling( 
  object = cipd_model_missing,
  data = list(
    n_t = n_t, 
    ts = 1:n_t,
    n_sims = n_sims, 
    dt = dt,
    sims_seq = sims_seq,
    n_states = n_states, 
    n_parameters = n_parameters,
    y = d_missing[,grep("noise",colnames(d))][,1:n_states],
    domestic_production = d_missing$domestic_production_noise,
    consumption = d_missing$consumption_noise,
    b_fix = b, 
    e_fix = e, 
    f_fix = f, 
    g_fix = g, 
    h_fix = h,
    missing_idx_C = missing_idx_C,
    missing_idx_I = missing_idx_I,
    missing_idx_D = missing_idx_D,
    missing_idx_P = missing_idx_P,
    n_missing_C = n_missing_C,
    n_missing_I = n_missing_I,
    n_missing_D = n_missing_D,
    n_missing_P = n_missing_P,
    #initial_values = initial_values,
    prior_only=0,
    rk4 = 0
    #p = p_list
  ), 
  #init = 0, 
  cores = 4, chains = 4, warmup = 250, iter = 500
  #,control = list(adapt_delta=0.85)#, max_treedepth=15)
)

print(fit_model_missing, probs = c(0.025, 0.975), digits_summary = 2, 
      pars=c("a","k","w","m","q","r","s","critical_ratio")
)
pairs(fit_model_missing, pars=c("a","k","w","m","q","r","s","lp__"))
par(mfrow=c(3,3))
rstan::traceplot(fit_model_missing, pars=c("a","k","w","m","q","r","s"))
par(mfrow=c(1,1))

draws_missing <- as.data.frame(fit_model_missing)

n_draws <- nrow(draws_missing)

# posterior predictions
n_samples <- 100
sample_draws <- sample(1:n_draws, n_samples, replace=F)

states_init <- draws_missing[ , grep("initial", colnames(draws_missing))]
states_draws <- draws_missing[ , grep("states", colnames(draws_missing))][,-c(1:4)]
par(mfrow=c(2,2))

for(n in 1:n_states){
  #int_states <- seq(1,n_sims,1/dt)
  #the_states <- states_draws[, ( ((n-1)*(n_sims)+1):(n_sims*n))][,int_states]
  the_states <- cbind(states_init[,n], states_draws[ , ( ((n-1)*(n_t)+1):(n_t*n))])
  
  initial <- run_cfs[1,1+n]
  raw_data_point_type <- ifelse(d_missing[,5+n] != -100, 16, 1)
  
  plot(1:(n_t+1), apply(the_states, 2, mean), type="n",
       ylim=c( min(apply(the_states, 2, rethinking::HPDI)[1,])*0.8, 
               max(apply(the_states, 2, rethinking::HPDI)[2,])*1.2)
       )
  points(1:(n_t+1), c(initial, d[,5+n]), col=scales::alpha("slateblue",0.6), pch=raw_data_point_type)
  points(1, initial, col="slateblue", pch=16)
  
  for(i in 1:n_samples){
    lines(1:(n_t+1), the_states[sample_draws[i], ], 
          col = scales::alpha("black", 0.2), lwd=0.5)
  }
  lines(1:(n_t+1), c(initial, d[,1+n]), type="l", lwd=3, col=scales::alpha("yellow",0.6))
  if(c(n_missing_C, n_missing_I, n_missing_D, n_missing_P)[n] > 0) {
    missing_idx <- list(missing_idx_C,missing_idx_I,missing_idx_D,missing_idx_P)[[n]]
    missing_imputations <- draws_missing[, grep(paste0(c("C","I","D","P")[n], "_impute"), colnames(draws_missing))]
    segments(x0 = missing_idx, x1 = missing_idx, 
             y0 = apply(missing_imputations, 2, rethinking::HPDI)[1,], 
             y1 = apply(missing_imputations, 2, rethinking::HPDI)[2,],
             lwd=0.5, col=scales::alpha("black",0.7)
    )
    points(missing_idx, apply(missing_imputations, 2, mean), 
           col=scales::alpha("red",1), pch=16)
  }
}
par(mfrow=c(1,1))





