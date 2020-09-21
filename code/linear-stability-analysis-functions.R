
alpha_critical <- function(p_list) (1+p_list$beta)*p_list$kappa*p_list$gamma/(p_list$omega + p_list$gamma/2)
kappa_critical <- function(p_list) p_list$alpha*(p_list$omega + p_list$gamma/2)/( (1+p_list$beta)*p_list$gamma)

get_fixed_points <- function( p_list ){
  
  alpha <- p_list$alpha
  beta <- p_list$beta
  gamma <- p_list$gamma
  kappa <- p_list$kappa
  omega <- p_list$omega
  delta <- p_list$delta
  
  if(alpha < alpha_critical(p_list)){
    return(list(
      "capital" = 0, 
      "inventory" = kappa*gamma/(omega + gamma/2),
      "demand" = kappa*gamma/(omega + gamma/2),
      "price" = (omega + gamma/2)/(kappa*gamma)
    ))
  }
  if(alpha > alpha_critical(p_list)){
    return(list(
      "capital" = (2*gamma*kappa*(-1-beta) + alpha*(gamma+2*omega))/(2*(1+beta)*(delta*(1-kappa))), 
      "inventory" = alpha/(1+beta),
      "demand" = alpha/(1+beta),
      "price" = (1+beta)/(alpha)
    ))
  }
}

get_jacobian <- function(C, I, D, P){
  matrix(c(
    # row 1
    alpha * P - 1 - beta, 0, 0, alpha*C, 
    # row 2
    C - kappa*delta, 
    -omega - (gamma*D^2)/(D+I)^2,  
    - (gamma*I^2)/(D+I)^2, 
    0,
    # row 3
    0, 0, -mu, -mu*P^2, 
    # row 4
    0, 0, rho*P*I^(-1), -rho*P*I^(-2)
  ), 
  ncol = 4, nrow=4, byrow = T
  )
}

get_eignvalues <- function(x) eigen(x)

# print(lambdas <- c(
#   alpha * (omega + gamma/2)/(kappa*gamma) - 1 - beta, 
#   - omega - gamma/4,
#   1/2*( -1*(rho * ((omega+gamma/2)/(kappa*gamma))^3 + mu) + 
#           sqrt( 
#             ((rho * ((omega+gamma/2)/(kappa*gamma))^3 + mu)^2 - 
#                4*mu*rho*((omega+gamma/2)/(kappa*gamma))^3*(1+((omega+gamma/2)/(kappa*gamma)))) 
#           ) 
#   ),
#   1/2*( -1*(rho * ((omega+gamma/2)/(kappa*gamma))^3 + mu) - 
#           sqrt( 
#             ((rho * ((omega+gamma/2)/(kappa*gamma))^3 + mu)^2 - 
#                4*mu*rho*((omega+gamma/2)/(kappa*gamma))^3*(1+((omega+gamma/2)/(kappa*gamma)))) 
#           ) 
#   )
# ))


# 
# 
# alpha_c <- function(beta, gamma, omega, kappa){
#   return(
#     (1 + beta)*kappa*gamma/(omega+gamma/2)
#   )
# }
# 
# kappa_c <- function(alpha, omega, gamma, beta){
#   alpha*(omega+gamma/2)/((1+beta)*gamma)
# }
# 
# factor_change_grid <- seq(0.1, 5, 0.1)
# res <- data.frame(
#   beta = rep(NA, length(factor_change_grid)),
#   gamma = rep(NA, length(factor_change_grid)),
#   omega = rep(NA, length(factor_change_grid)),
#   kappa = rep(NA, length(factor_change_grid))
# )
# 
# initial_param_vec <- c(beta=beta, gamma=gamma, omega=omega, kappa=0.5)
# 
# for(i in seq_along(factor_change_grid)){
#   for(j in 1:ncol(res)){
#     if(j == which(names(initial_param_vec) == "kappa") && 
#        inital_param_vec["kappa"][[1]]*factor_change_grid[i] > 1.0) break() 
#     param_vec <- initial_param_vec
#     param_vec[j] <- param_vec[j]*factor_change_grid[i]
#     res[i, j] <- alpha_c(beta = param_vec["beta"][[1]], 
#                          gamma = param_vec["gamma"][[1]],
#                          kappa = param_vec["kappa"][[1]],
#                          omega = param_vec["omega"][[1]]
#                          )
#   }
# }
# 
# plot(factor_change_grid, rnorm(length(factor_change_grid)), 
#      type="n", ylim=c(min(unlist(res), na.rm = T), max(unlist(res), na.rm=T)+mean(unlist(res), na.rm=T)), 
#      ylab = "alpha_c")
# for(i in 1:length(initial_param_vec)){
#   lines(factor_change_grid, res[,i], type="l", lty=i, lwd=2)
# }
# legend(
#   "topright",
#   bty="n",
#   legend = names(initial_param_vec),
#   lty = 1:4, 
#   lwd=2
# )
# 
# ## (alpha, kappa) space 
# n_grid <- 100
# kappa_grid <- seq(0.001, 0.999, length.out = n_grid)
# alpha_grid <- seq(0.001, 5, length.out = n_grid)
# kappa_alpha_res <- data.frame(
#   kappa = kappa_grid, 
#   alpha = alpha_grid, 
#   kappa_c = rep(NA, n_grid),
#   alpha_c = rep(NA, n_grid)
# )
# for(i in 1:n_grid){
#   kappa_alpha_res$alpha_c[i] <- alpha_c(
#     beta = 0.2, 
#     gamma = 50, 
#     omega = 1.2, 
#     kappa = kappa_grid[i]
#   )
#   kappa_alpha_res$kappa_c[i] <- kappa_c(
#     alpha = alpha_grid[i], 
#     omega = 1.2, 
#     gamma = 50, 
#     beta = 0.2
#   )
# }
# 
# plot(kappa_alpha_res$kappa, kappa_alpha_res$alpha, type="n", xlim=c(0,1), ylim=c(0,3))
# #lines(kappa_alpha_res$kappa_c, kappa_alpha_res$alpha_c, type="l")
# polygon(x = c(min(kappa_alpha_res$kappa_c), 
#               min(kappa_alpha_res$kappa_c), 
#               max(kappa_alpha_res$kappa_c), 
#               max(kappa_alpha_res$kappa_c)), 
#         y = c(
#           0, kappa_alpha_res[kappa_alpha_res$kappa_c %in% min(kappa_alpha_res$kappa_c), "alpha_c"],
#           kappa_alpha_res[kappa_alpha_res$kappa_c %in% max(kappa_alpha_res$kappa_c), "alpha_c"],
#           0
#         ), 
#         col = scales::alpha("green", 0.5), 
#         border = scales::alpha("green", 0.5)
# )
# 
