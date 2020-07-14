# Simulate data from the model and fit in Stan

source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/ode-model.R")
source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/linear-stability-analysis-functions.R")

# set the parameters
a <- 1/52
b <- 130
e <- 1/(2.5*52)
f <- 1/52
g <- 110*24*0.75
k <- 0.5
h <- 30e6
w <- 0.3
m <- 1/10
q <- 130
r <- 1/20
s <- 1
C_init <- 100e3
I_init <- 30e6
D_init <- 30e6
P_init <- 140

states <- c(C = C_init, I = I_init, D = D_init, P = P_init)

dt <- 0.1
n_t <- 50
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
set.seed(2020)
d$C_noise <- rlnorm( nrow(d), log(d$C), 0.1) #d$C + rnorm(nrow(d), 0, 1000)
d$I_noise <- rlnorm( nrow(d), log(d$I), 0.1) #d$I + rnorm(nrow(d), 0, 1e5)
d$D_noise <- rlnorm( nrow(d), log(d$D), 0.1) #d$D + rnorm(nrow(d), 0, 1e5)
d$P_noise <- rlnorm( nrow(d), log(d$P), 0.1) #d$P + rnorm(nrow(d), 0, 5)

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
dt <- dt;
n_sims <- n_t * dt^-1 + 1;
sims_seq <- seq(1,n_sims,1/dt)[-1]
n_states <- 4
n_parameters <- 12

fit_data <- rstan::stan(file = "~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/cidp-model.stan", 
                        data = list(
                          n_t = n_t, 
                          n_sims = n_sims, 
                          dt = dt,
                          sims_seq = sims_seq,
                          n_states = n_states, 
                          n_parameters = n_parameters,
                          y = d[,grep("noise",colnames(d))],
                          #initial_values = initial_values,
                          prior_only=0
                          #p = p_list
                        ), 
                        init = 0,
                        #algorithm = "Fixed_param", 
                        cores = 4, chains = 4, warmup = 200, iter = 400
                        #,control = list(adapt_delta=0.9, max_treedepth=15)
                        )

print(fit_data, probs = c(0.025, 0.975), digits_summary = 2, pars="p")
pairs(fit_data, pars=names(p_list))

draws <- as.matrix(fit_data)

n_draws <- nrow(draws)

# parameter estimates
par(mfrow=c(3, 4))
for(i in 1:length(p_list)){
  plot(density(draws[,grep("p", colnames(draws))][,i]), main=paste0("p ", i), xlab="", ylab="", lwd=3)
  abline(v=p_list[i], lty=2, lwd=2, col="blue")
}
par(mfrow=c(1,1))

# posterior predictions
n_samples <- 100
sample_draws <- sample(1:n_draws, n_samples, replace=F)

#states_init <- draws[ , grep("initial", colnames(draws))]
states_draws <- draws[ , grep("states", colnames(draws))][,-c(1:4)]
par(mfrow=c(2,2))

for(n in 1:n_states){
  int_states <- seq(1,n_sims,1/dt)
  the_states <- states_draws[, ( ((n-1)*(n_sims)+1):(n_sims*n))][,int_states]
  
  inital <- run_cfs[1,1+n]
  
  plot(1:(n_t+1), c(inital, d[,5+n]), type="n" )
  points(1:(n_t+1), c(inital, d[,5+n]), col="blue")
  points(1, inital, col="blue", pch=16)
  for(i in 1:n_samples){
    lines(1:(n_t+1), the_states[sample_draws[i], ], 
          col = scales::alpha("black", 0.2), lwd=0.5)
  }
}
par(mfrow=c(1,1))

