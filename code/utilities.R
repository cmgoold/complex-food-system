
get_needed_packages <- function(x){
  needed <- c("rstan", "coda", "scales", "RColorBrewer")
  if(sum(installed.packages() %in% needed) < length(needed)){
    install.packages( needed[which( !(needed %in% installed.packages()))] )
  } 
  else{
    "All packages are already installed"
  }
}

#-----------------------------------------------------------------------------------------------
# alpha
col_alpha <- function(col, alpha) scales::alpha(col, alpha = alpha)

#-----------------------------------------------------------------------------------------------
# hdi function
hdi <- function(x, prob=0.95) coda::HPDinterval(coda::as.mcmc(x), prob=prob)
