#####################################################################
#   Fit the complex food systems model to the
#   UK pork industry data set
#   Copyright Conor Goold (2020)
#####################################################################

source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/utilities.R")

get_needed_packages()

# load the data
d_uk_pork <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/data/uk_pork_industry_annual_data.csv")

#--- arrange the data

# time variable
d_uk_pork$time <- 1:nrow(d_uk_pork)

# make estimated total consumption from ONS average weekly per-capita consumption

# UK population
population_size <- 66.65e6 
# 52 weeks in a year
weeks_year <- 52
# grams to thousand tonnes conversion factor = g * 1/1000 kg * 1/1000 * tonnes
grams_thousand_tonnes_conversion_factor <- 1/1000^3
# add to data
d_uk_pork$total_national_consumption_thousand_tonnes <- (
  d_uk_pork$total_ppc_average_g_ONS * weeks_year * population_size * grams_thousand_tonnes_conversion_factor
)

# add an empty demand column to indicate it's missing
d_uk_pork$total_national_demand_thousand_tonnes <- NA

#-----------------------------------------------------------------------------------------------------------
# FIT THE STAN MODEL

# generate the missing numbers and indexes
missing_idx_C <- which(is.na(d_uk_pork$breeding_herd_thousand_head))
n_missing_C <- length(missing_idx_C)

missing_idx_slaughtered <- which(is.na(d_uk_pork$clean_pigs_slaughtered_thousand_head))
n_missing_slaughtered <- length(missing_idx_slaughtered)

missing_idx_I <- which(is.na(d_uk_pork$total_supply_thousand_tonnes_carcass_weight_kg))
n_missing_I <- length(missing_idx_I)

missing_idx_D <- which(is.na(d_uk_pork$total_national_demand_thousand_tonnes))
n_missing_D <- length(missing_idx_D)

missing_idx_consumption <- which(is.na(d_uk_pork$total_national_consumption_thousand_tonnes))
n_missing_consumption <- length(missing_idx_consumption)

missing_idx_P <- which(is.na(d_uk_pork$pig_price_pence_kg))
n_missing_P <- length(missing_idx_P)

# get the state variables, transform NAs to -100 for Stan
state_variables <- apply(with(d_uk_pork, cbind(breeding_herd_thousand_head, 
                                          total_supply_thousand_tonnes_carcass_weight_kg,
                                          total_national_demand_thousand_tonnes,
                                          pig_price_pence_kg
                                          )), 
                         2,
                         function(x) ifelse(is.na(x), -100, x)
)

par(mfrow=c(2,2))
for(n in 1:4) {
  plot(state_variables[,n], type="b")
  point_colours <- ifelse(state_variables[,n] == -100, "black", "slateblue")
  points(1:nrow(d_uk_pork), state_variables[,n], pch=16, col=point_colours)
}
par(mfrow=c(1,1))

slaughtered <- ifelse(is.na(d_uk_pork$clean_pigs_slaughtered_thousand_head), -100, d_uk_pork$clean_pigs_slaughtered_thousand_head)
consumption <- ifelse(is.na(d_uk_pork$total_national_consumption_thousand_tonnes), -100, d_uk_pork$total_national_consumption_thousand_tonnes)

prior_only = 0

stan_data <- list(
  n_t = nrow(d_uk_pork),
  ts = 1:nrow(d_uk_pork),
  n_states = 4, 
  n_parameters = 12, 
  y = state_variables, 
  slaughtered = slaughtered, 
  consumption = consumption, 
  n_missing_C = n_missing_C,
  n_missing_I = n_missing_I,
  n_missing_D = n_missing_D,
  n_missing_P = n_missing_P,
  n_missing_slaughtered = n_missing_slaughtered, 
  n_missing_consumption = n_missing_consumption,
  missing_idx_C = missing_idx_C,
  missing_idx_I = missing_idx_I,
  missing_idx_D = missing_idx_D,
  missing_idx_P = missing_idx_P,
  missing_idx_slaughtered = missing_idx_slaughtered, 
  missing_idx_consumption = missing_idx_consumption,
  prior_only = prior_only
)

# compile the model
cipd_model <- rstan::stan_model("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/cidp-uk-pork-industry.stan")

inits <- function(){
  list(
    a = runif(1,0.11,1), 
    w = runif(1,0.1,1),
    k = runif(1,0.1,1),
    m = runif(1,0.1,1),
    r = runif(1,0.1,1),
    q = runif(1,0.1,2)
  )
}

# fit the model
fit_uk_pork_model <- rstan::sampling( 
  object = cipd_model,
  data = stan_data, 
  init = inits, 
  cores = 4, chains = 4, warmup = 200, iter = 400,
  control = list(max_treedepth=12, stepsize=0.01)
)

print(fit_uk_pork_model, pars=c("a","k","w","m","r","s","critical_ratio"), 
      digits_summary = 2, probs=c(0.025, 0.975)
      )
par(mfrow=c(4,2))
rstan::traceplot(fit_uk_pork_model, pars=c("a","k","w","m","q","r","s"))
par(mfrow=c(1,1))
pairs(fit_uk_pork_model, pars=c("a","k","w","m","q","r","s","lp__"))

draws <- as.data.frame(fit_uk_pork_model)
n_draws <- nrow(draws)

# posterior predictions
n_samples <- 200
sample_draws <- sample(1:n_draws, n_samples, replace=F)

states_init <- draws[ , grep("initial", colnames(draws))]
states_draws <- draws[ , grep("states", colnames(draws))][,-c(1:4)]
state_variables_plot <- apply(state_variables, 2, function(x) ifelse(x==-100, NA, x))
n_t <- nrow(d_uk_pork)

# plot
par(mfrow=c(2,2))
for(n in 1:n_states){
  the_states <- cbind(states_init[,n], states_draws[ , ( ((n-1)*(n_t)+1):(n_t*n))])
  
  initial <- NA

  plot(c(inital, state_variables_plot[,n]), type="n" )
  points(c(inital, state_variables_plot[,n]), col=scales::alpha("slateblue",0.6), pch=16)
  for(i in 1:n_samples){
    lines(1:(n_t+1), the_states[sample_draws[i], ], 
          col = scales::alpha("black", 0.2), lwd=0.5)
  }
}
par(mfrow=c(1,1))

