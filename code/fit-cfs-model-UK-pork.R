#####################################################################
#   Fit the complex food systems model to the UK pork industry data
#   using  Stan
#   Copyright Conor Goold (2020)
#####################################################################
code_path <- "~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/"
data_path <- "~/Dropbox/Leeds_postdoc/Papers/cfs-model/data/"

source(paste0(code_path, "utilities.R"))
source(paste0(code_path, "ode-model.R"))

get_needed_packages()

# load the data
d_monthly_uk_pork <- read.csv(paste0(data_path, "uk_pork_industry_monthly_data.csv"))

# change tonnes to kg
d_monthly_uk_pork$imports_kg <- d_monthly_uk_pork$imports_tonnes * 1000
d_monthly_uk_pork$exports_kg <- d_monthly_uk_pork$exports_tonnes * 1000

# calculate estimated demand as UK production - exports + imports
d_monthly_uk_pork$total_supply_kg <- with(d_monthly_uk_pork, UK_production_kg - exports_kg + imports_kg)

# matrix of state variables
state_variables <- with(d_monthly_uk_pork, 
                        cbind( breeding_herd_head, total_supply_kg, NA, all_pig_price_p_kg)
)

# change NA to -100 so as to be read by Stan
state_variables <- apply(state_variables, 2, function(x) ifelse( is.na(x), -100, x))

stan_data <- list(
  n_t = nrow(state_variables),
  ts = 1:nrow(state_variables),
  n_states = 4, 
  n_parameters = 12, 
  y = state_variables, 
  uk_production = d_monthly_uk_pork$UK_production_kg, 
  trade = d_monthly_uk_pork$imports_kg - d_monthly_uk_pork$exports_kg,
  imports = d_monthly_uk_pork$imports_kg,
  exports = d_monthly_uk_pork$exports_kg,
  prior_only = 0
)

# compile the model
cipd_model <- rstan::stan_model(paste0(code_path, "Stan-cfs-model-UK-pork.stan"))

# fit the model
fit_uk_pork_model <- rstan::sampling( 
  object = cipd_model,
  data = stan_data, 
  init = 0, 
  cores = 4, chains = 4, warmup = 2500, iter = 5000,
  seed = 2020,
  control = list(max_treedepth=10, adapt_delta=0.99)
)

capture.output(rstan::check_hmc_diagnostics(fit_uk_pork_model), 
               file = paste0(code_path, "Stan-MCMC-diagnostics.csv")
               )

draws <- as.data.frame(fit_uk_pork_model)

write.csv(draws, file = paste0(code_path, "Stan-MCMC-results2.csv"), row.names = FALSE)

capture.output(
  print(fit_uk_pork_model, pars=c("p", "a", "e_", "f", "k", "h", "w", "m","r","q","s","critical_ratio",
                                  "states_init",
                                  "sigma", "sigma_trade", "sigma_production"), 
        digits_summary = 4, probs=c(0.025, 0.975)
  ),
  file = paste0(code_path, "Stan-MCMC-summary2.csv")
)
