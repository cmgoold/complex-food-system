source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/linear-stability-analysis-functions.R")
source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/utilities.R")
setwd("~/Dropbox/Leeds_postdoc/Papers/cfs-model/paper/figures")

# bifurcation diagram for capital
png("figure_2a.png", width=1000, height=1000, res=200)
alpha_grid <- seq(0.1, 2, length.out = 100)
#cr <- seq(0, 2, length.out = 100)
p_list <- list(beta = 0.2, gamma = 26, omega = 10, kappa = 0.5, delta = 5)
fp <- numeric(length(alpha_grid))
for(i in 1:length(fp)) {
  fp[i] <- get_fixed_points(p_list = list(
      alpha = alpha_grid[i], 
      beta = p_list$beta, 
      gamma = p_list$gamma,
      omega = p_list$omega, 
      kappa = p_list$kappa,
      delta = p_list$delta
    )
  )$capital
}
plot(alpha_grid, fp, type="n", axes=F, 
     ylab = "capital (dimensionless)", 
     xlab = expression(alpha), cex.lab=1.5
)
axis(side=1, lwd=2, at = c(min(alpha_grid), mean(alpha_grid), max(alpha_grid)))
axis(side=2, lwd=2)
lines(alpha_grid[1:(min(which(fp>0))-1)], fp[1:(min(which(fp>0))-1)], lwd=3, 
      col=scales::alpha("slateblue", 0.8))
lines(alpha_grid[(min(which(fp>0))):(length(fp))], rep(0, length(alpha_grid[(min(which(fp>0))):(length(fp))])), 
      lwd=3, lty=2
)
lines(alpha_grid[(min(which(fp>0))):(length(fp))], fp[(min(which(fp>0))):(length(fp))], 
      col=scales::alpha("slateblue", 0.8), 
      lwd=3
      )
#abline(v=alpha_critical(p_list), lty=3, col=scales::alpha("black",0.7))
legend("topleft", 
       bty = "n",
       legend = c("stable", "unstable"), lwd = 3, 
       col = c(scales::alpha("slateblue", 0.8), "black"), 
       lty=c(1,2), 
       cex = 1.3)
mtext("a", line = 0, side = 3, at = 0, font=2, cex=1.25)
dev.off()

# (alpha, kappa) parameter space
png("figure_2b.png", width=1000, height=1000, res=200)
alpha_grid <- seq(0.01, 2, length.out = 100)
kappa_grid <- seq(0.01, 0.99, length.out = 100)
plot(kappa_grid, alpha_grid, type="n", 
     xlab = expression(kappa), ylab=expression(alpha), 
     cex.lab = 1.5, 
     axes=F
     )
ac <- numeric(length(alpha_grid)) 
the_beta <- 0.1
the_gamma <- 26
the_omega <- 10
for(i in 1:length(ac)){
  ac[i] <- alpha_critical(p_list = list(kappa=kappa_grid[i], beta=the_beta, gamma=the_gamma, omega=the_omega))
}
polygon(
  x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
  y = c(min(kappa_grid), ac[1], ac[length(ac)], min(kappa_grid)), 
  col = scales::alpha("red", 0.6), 
  border = scales::alpha("red", 0.2)
)
polygon(x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
        y = c(
          ac[1], 2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega),
          2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega),
          ac[length(ac)]
        ),
        col = scales::alpha("steelblue", 0.8),
        border = scales::alpha("steelblue", 0.1)
)
polygon(x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
        y = c(
          2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega), 2,
          2, 2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega)
        ),
        col = scales::alpha("green", 0.6),
        border = scales::alpha("green", 0.6)
)
axis(side = 1, at = c(0,1), lwd=2, cex.lab=1.2)
axis(side = 2, at = c(0,1,2), lwd=2, cex.lab=1.2)
text(x = 0.75, y = 0.5, labels = "unsustainable")
text(x = 0.3, y = 0.8, labels = "sustainable (imports)")
text(x = 0.3, y = 1.7, labels = "sustainable (exports)")
mtext(text = expression(paste(beta, " = ") ~ 0.1), side = 3, line = 0, at = 0.9, cex=1.25)
mtext("b", line = 0, side = 3, at = 0, font=2, cex=1.25)
dev.off()


# (alpha, kappa) parameter space
png("figure_2c.png", width=1000, height=1000, res=200)
alpha_grid <- seq(0.01, 2, length.out = 100)
kappa_grid <- seq(0.01, 0.99, length.out = 100)
plot(kappa_grid, alpha_grid, type="n", 
     xlab = expression(kappa), ylab=expression(alpha), 
     cex.lab = 1.5, axes=FALSE
)
axis(side = 1, at = c(0,1), lwd=2, cex.lab=1.2)
axis(side = 2, at = c(0,1,2), lwd=2, cex.lab=1.2)
ac <- numeric(length(alpha_grid)) 
the_beta <- 0.6
the_gamma <- 26
the_omega <- 10
for(i in 1:length(ac)){
  ac[i] <- alpha_critical(p_list = list(kappa=kappa_grid[i], beta=the_beta, gamma=the_gamma, omega=the_omega))
}
polygon(
  x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
  y = c(min(kappa_grid), ac[1], ac[length(ac)], min(kappa_grid)), 
  col = scales::alpha("red", 0.6), 
  border = scales::alpha("red", 0.2)
)
polygon(x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
        y = c(
          ac[1], 2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega),
          2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega),
          ac[length(ac)]
        ),
        col = scales::alpha("steelblue", 0.8),
        border = scales::alpha("steelblue", 0.1)
)
polygon(x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
        y = c(
          2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega), 3,
          3, 2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega)
        ),
        col = scales::alpha("green", 0.6),
        border = scales::alpha("green", 0.6)
)
#text(x = 0.75, y = 0.5, labels = "unsustainable")
#text(x = 0.3, y = 1.3, labels = "sustainable (imports)")
#text(x = 0.3, y = 1.9, labels = "sustainable (exports)")
mtext(text = expression(paste(beta, " = ") ~ 0.6), side = 3, line = 0, at = 0.9, cex=1.25)
mtext("c", line = 0, side = 3, at = 0, font=2, cex=1.25)
dev.off()


# (alpha, kappa) parameter space
png("figure_2d.png", width=1000, height=1000, res=200)
alpha_grid <- seq(0.01, 2, length.out = 100)
kappa_grid <- seq(0.01, 0.99, length.out = 100)
plot(kappa_grid, alpha_grid, type="n", 
     xlab = expression(kappa), ylab=expression(alpha), 
     cex.lab = 1.5, axes=FALSE
)
axis(side = 1, at = c(0,1), lwd=2, cex.lab=1.2)
axis(side = 2, at = c(0,1,2), lwd=2, cex.lab=1.2)
ac <- numeric(length(alpha_grid)) 
the_beta <- 1
the_gamma <- 26
the_omega <- 10
for(i in 1:length(ac)){
  ac[i] <- alpha_critical(p_list = list(kappa=kappa_grid[i], beta=the_beta, gamma=the_gamma, omega=the_omega))
}
polygon(
  x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
  y = c(min(kappa_grid), ac[1], ac[length(ac)], min(kappa_grid)), 
  col = scales::alpha("red", 0.6), 
  border = scales::alpha("red", 0.2)
)
polygon(x = c(min(kappa_grid), min(kappa_grid), max(kappa_grid), max(kappa_grid)),
        y = c(
          ac[1], 2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega),
          2*the_gamma*(1 + the_beta)/(the_gamma + 2*the_omega),
          ac[length(ac)]
        ),
        col = scales::alpha("steelblue", 0.8),
        border = scales::alpha("steelblue", 0.1)
)
#text(x = 0.75, y = 0.5, labels = "unsustainable")
#text(x = 0.3, y = 1.3, labels = "sustainable (imports)")
#text(x = 0.3, y = 1.9, labels = "sustainable (exports)")
mtext(text = expression(paste(beta, " = ") ~ 1), side = 3, line = 0, at = 0.9, cex=1.25)
mtext("d", line = 0, side = 3, at = 0, font=2, cex=1.25)
dev.off()

#######################################################################
# model results
# load the data
d_monthly_uk_pork <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/data/uk_pork_industry_monthly_data.csv")

d_monthly_uk_pork$imports_kg <- d_monthly_uk_pork$imports_tonnes * 1000
d_monthly_uk_pork$exports_kg <- d_monthly_uk_pork$exports_tonnes * 1000
d_monthly_uk_pork$total_supply_kg <- with(d_monthly_uk_pork, UK_production_kg - exports_kg + imports_kg)

state_variables <- with(d_monthly_uk_pork, 
                        cbind( breeding_herd_head, NA, total_supply_kg, all_pig_price_p_kg)
)

state_variables <- apply(state_variables, 2, function(x) ifelse( is.na(x), -100, x))

draws <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/Stan-MCMC-results.csv")
n_draws <- nrow(draws)

# posterior predictions
n_samples <- 200
sample_draws <- sample(1:n_draws, n_samples, replace=F)

states_init <- draws[ , grep("initial", colnames(draws))]
states_draws <- draws[ , grep("states", colnames(draws))][,-c(1:4)]
state_variables_plot <- apply(state_variables, 2, function(x) ifelse(x==-100, NA, x))
n_t <- nrow(d_monthly_uk_pork)

# Figure 3a
png("figure_3a.png", width=1300, height=1000, res=200)
the_state <- cbind(states_init[,1], states_draws[ , ( ((1-1)*(n_t)+1):(n_t*1))])
initial <- NA
plot(c(initial, state_variables_plot[,1]), type="n", 
     ylim=c(min(na.omit(state_variables_plot[,1])) * 0.87 , max(na.omit(state_variables_plot[,1])) * 1.13 ),
     axes=F, xlab = "months/years", ylab = "breeding herd (head)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
# axis(side=1, at = seq(2, n_t+1, 6), 
#      labels = rep(c("Jan", "Jul"), times=5), cex.axis=0.8)
# abline(v=seq(2,n_t+1,12), lty=2, col=scales::alpha("black",0.5))
axis(side=2, at = seq(300e3, 550e3,50e3))
points(c(initial, state_variables_plot[,1]), col=scales::alpha("slateblue",0.6), pch=16)
for(i in 1:n_samples){
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.2), lwd=0.5)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=3)
mtext("a", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3b
png("figure_3b.png", width=1300, height=1000, res=200)
the_state <- cbind(states_init[,2], states_draws[ , ( ((2-1)*(n_t)+1):(n_t*2))])
initial <- NA
plot(rnorm(n_t+1), type="n", ylim=c(5e7, 45e7),
     axes=F, xlab = "months/years", ylab = "inventory (kg)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(5e7, 45e7, 10e7))
for(i in 1:n_samples){
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.2), lwd=0.5)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=3)
mtext("b", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3c
png("figure_3c.png", width=1300, height=1000, res=200)
the_state <- cbind(states_init[,3], states_draws[ , ( ((3-1)*(n_t)+1):(n_t*3))])
initial <- NA
plot(c(initial, state_variables_plot[,3]), type="n", 
     ylim=c(min(na.omit(state_variables_plot[,3])) * 0.9 , max(na.omit(state_variables_plot[,3])) * 1.1 ),
     axes=F, xlab = "months/years", ylab = "demand (kg)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(10e7, 20e7, 3e7))
points(c(initial, state_variables_plot[,3]), col=scales::alpha("slateblue",0.6), pch=16, type="b")
for(i in 1:n_samples){
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.2), lwd=0.5)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=3)
mtext("c", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3d
png("figure_3d.png", width=1300, height=1000, res=200)
the_state <- cbind(states_init[,4], states_draws[ , ( ((4-1)*(n_t)+1):(n_t*4))])
initial <- NA
plot(c(initial, state_variables_plot[,4]), type="n", 
     ylim=c(min(na.omit(state_variables_plot[,4])) * 0.9 , max(na.omit(state_variables_plot[,4])) * 1.1 ),
     axes=F, xlab = "months/years", ylab = "price (p/kg)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(110, 180, 10))
points(c(initial, state_variables_plot[,4]), col=scales::alpha("slateblue",0.6), pch=16, type="b")
for(i in 1:n_samples){
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.2), lwd=0.5)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=3)
mtext("d", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3e
png("figure_3e.png", width=1300, height=1000, res=200)
states_bh <- cbind(states_init[,1], states_draws[ , ( ((1-1)*(n_t)+1):(n_t))])
plot(c(NA, d_monthly_uk_pork$UK_production_kg), type="b",
     ylim=c(min(na.omit(d_monthly_uk_pork$UK_production_kg)) * 0.95 , max(na.omit(d_monthly_uk_pork$UK_production_kg)) * 1.05 ) 
     , col = scales::alpha("slateblue",0.6), pch=16,
     axes=F, xlab = "months/years", ylab = "UK production (kg)", cex.lab=1.5)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(65e6, 95e6, 5e6))
for(i in 1:n_samples){
  lines(1:(n_t+1), states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."], 
        col = scales::alpha("black", 0.2), lwd=0.5)
}
lines(1:(n_t+1), apply( states_bh * draws$p.4. * draws$p.5., 2, mean), 
      col = scales::alpha("black", 1), lwd=3)
mtext("e", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3f
png("figure_3f.png", width=1300, height=1000, res=200)
plot(c(NA, d_monthly_uk_pork$exports_kg), type="b", 
     col = scales::alpha("slateblue",0.6), pch=16,
     ylim=c(min(na.omit(d_monthly_uk_pork$exports_kg)) * 0.9 , max(na.omit(d_monthly_uk_pork$exports_kg)) * 1.1 ),
     axes=F, xlab = "months/years", ylab = "Exports (kg)", cex.lab=1.5)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(15e6, 45e6, 5e6))
for(i in 1:n_samples){
  lines(1:(n_t+1), draws[sample_draws[i], "p.6."] * states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."], 
        col = scales::alpha("black", 0.2), lwd=0.5)
}
lines(1:(n_t+1), apply( draws$p.6. * states_bh * draws$p.4. * draws$p.5., 2, mean), 
      col = scales::alpha("black", 1), lwd=3)
mtext("f", line = 1, at = 1, font=2, cex=1.5)
dev.off()
