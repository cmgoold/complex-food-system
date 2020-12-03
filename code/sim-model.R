code_path <- "~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/"

source(paste0(code_path, "linear-stability-analysis-functions.R"))
source(paste0(code_path, "ode-model.R"))

# set the parameters
a <- 1/52
b <- 120
e <- 1/(2.5*52)
f <- 1/52
g <- 110*24*0.75
k <- 0.6
h <- 30e6
w <- 0.3
m <- 1/2
q <- 150
r <- 1/20
s <- 1
C_init <- 100e3
I_init <- 30e6
D_init <- 30e6
P_init <- 140

states <- c(C = C_init, I = I_init, D = D_init, P = P_init)

dt <- 0.1
t <- seq(0, 52*50, dt)

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

run_cfs[ nrow(run_cfs), ]

# define dimensionless variables
alpha <- q/b
beta <- e/a
delta <- f*g*C_init/(a*h*s)
omega <- w/a
gamma <- 1/(a*s)
kappa <- k
mu <- m/a
rho <- r/a
Chat <- C_init/C_init
Ihat <- I_init/(h*s)
Dhat <- D_init/h
Phat <- P_init/q
that <- t/(1/a) 

states_nd <- c(Chat=Chat, Ihat=Ihat, Dhat=Dhat, Phat=Phat)

#cr(alpha = alpha, kappa = kappa, gamma = gamma, beta = beta, omega = omega)

critical_ratio(p_list = list(alpha = alpha, beta=beta, kappa=kappa, gamma=gamma, omega=omega))
alpha_critical(p_list = list(beta=beta, kappa=kappa, gamma=gamma, omega=omega))
kappa_critical(p_list = list(alpha=alpha, gamma=gamma, beta=beta, omega=omega))
get_fixed_points(p_list = list(alpha=alpha, beta=beta, kappa=kappa, gamma=gamma, omega=omega, delta=delta))
alpha*(gamma + 2*omega)/(2 * gamma * (1 + beta)) # if alpha > alpha_critical, and if alpha > this value, domestic production exceeds reference demand (exports)

run_cfs_nd <- as.data.frame(
  deSolve::ode(
    func = m_cfs_nd, 
    y = states_nd, 
    p = c(),
    times = that,
    method = "rk4"
  )
)

run_cfs_nd[nrow(run_cfs_nd), ]

plot(gamma - delta*run_cfs_nd$Chat, type="l")

par(mfrow=c(2,2))
plot(run_cfs_nd$time, run_cfs_nd$Chat, type="l", bty="l")
plot(run_cfs_nd$Ihat, type="l", bty="l")
plot(run_cfs_nd$Dhat, type = "l", bty="l")
plot(run_cfs_nd$Phat, type = "l", bty="l")
par(mfrow=c(1,1))

