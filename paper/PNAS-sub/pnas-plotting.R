code_path <- "~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/"
figure_path <- "~/Dropbox/Leeds_postdoc/Papers/cfs-model/paper/PNAS-sub/figures/"

source(paste0(code_path, "linear-stability-analysis-functions.R"))
source(paste0(code_path, "utilities.R"))
setwd(figure_path)

#figure 1a
n_seq <- 20
n_pars <- 7
par_seq <- c("q", "b", "e", "a", "w", "s", "k")
q_ref <- 160
b_ref <- 130
e_ref <- 1/(12*2.5)
a_ref <- 1/5
w_ref <- 1/3
s_ref <- 1
k_ref <- 0.6

p_list_ref <- list(
  q = q_ref, 
  b = b_ref, 
  e = e_ref, 
  a = a_ref, 
  w = w_ref,  
  s = s_ref, 
  k = k_ref
)

change_fac <- 0.9
p_list_seq <- lapply(p_list_ref, function(x) seq(x*(1-change_fac), x*(1+change_fac), length.out = n_seq) )

# container
cr_res <- rep(list(list()), n_pars)

# fill the container with critical ratios, varying the focal parameter each time
for(i in 1:n_pars){
  p_list_ <- p_list_ref
  out <- numeric(n_seq)
  focal_par_ <- par_seq[i]
  for(n in 1:n_seq){
    p_list_[focal_par_][[1]] <- p_list_seq[focal_par_][[1]][[n]]
    out[n] <- critical_ratio_raw(p_list_)
  }
  cr_res[[i]] <- out
}

cols_ <- RColorBrewer::brewer.pal(n_pars, "Paired")

#png("figure_2a.png", width=1000, height=1000, res=200)
png(paste0(figure_path, "figure_2.png"), width=2000, height=2000, res=200)
par(mfrow=c(2,2), mar=c(5.1, 4.7, 2, 1))
plot(
  0, 0, type="n", bty="n",
  xlab = "parameter/reference", ylab = "critical ratio",
  ylim = c(0,3), xlim = c(0,3),
  cex.lab = 1.5
)
abline(h=1, lty=2)
for(i in 1:n_pars){
  lines(p_list_seq[[i]]/jitter(p_list_ref[[i]], factor = runif(1,1.5,3)), cr_res[[i]], col=cols_[i], lwd=2)
}
legend(
  "topright", 
  legend = par_seq,
  col = cols_,
  lwd = 2,
  bty = "n"
)
mtext("a", side = 3, at = 0, font=2, cex=1.25)
#dev.off()

# bifurcation diagram for capital
# alpha_grid <- seq(0.1, 2, length.out = 100)
# #cr <- seq(0, 2, length.out = 100)
# p_list <- list(beta = 0.2, gamma = 26, omega = 10, kappa = 0.5, delta = 5)
# fp <- numeric(length(alpha_grid))
# for(i in 1:length(fp)) {
#   fp[i] <- get_fixed_points(p_list = list(
#       alpha = alpha_grid[i], 
#       beta = p_list$beta, 
#       gamma = p_list$gamma,
#       omega = p_list$omega, 
#       kappa = p_list$kappa,
#       delta = p_list$delta
#     )
#   )$capital
# }
# plot(alpha_grid, fp, type="n", axes=F, 
#      ylab = "capital (dimensionless)", 
#      xlab = expression(alpha), cex.lab=1.5
# )
# axis(side=1, lwd=2, at = c(min(alpha_grid), mean(alpha_grid), max(alpha_grid)))
# axis(side=2, lwd=2)
# lines(alpha_grid[1:(min(which(fp>0))-1)], fp[1:(min(which(fp>0))-1)], lwd=3, 
#       col=scales::alpha("slateblue", 0.8))
# lines(alpha_grid[(min(which(fp>0))):(length(fp))], rep(0, length(alpha_grid[(min(which(fp>0))):(length(fp))])), 
#       lwd=3, lty=2
# )
# lines(alpha_grid[(min(which(fp>0))):(length(fp))], fp[(min(which(fp>0))):(length(fp))], 
#       col=scales::alpha("slateblue", 0.8), 
#       lwd=3
#       )
# #abline(v=alpha_critical(p_list), lty=3, col=scales::alpha("black",0.7))
# legend("topleft", 
#        bty = "n",
#        legend = c("stable", "unstable"), lwd = 3, 
#        col = c(scales::alpha("slateblue", 0.8), "black"), 
#        lty=c(1,2), 
#        cex = 1.3)
# mtext("a", line = 0, side = 3, at = 0, font=2, cex=1.25)

# (alpha, kappa) parameter space
#png("figure_2b.png", width=1000, height=1000, res=200)
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
#dev.off()


# (alpha, kappa) parameter space
#png("figure_2c.png", width=1000, height=1000, res=200)
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
#dev.off()


# (alpha, kappa) parameter space
#png("figure_2d.png", width=1000, height=1000, res=200)
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
#dev.off()

dev.off()

#######################################################################
# model results
# load the data
d_monthly_uk_pork <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/data/uk_pork_industry_monthly_data.csv")

d_monthly_uk_pork$imports_kg <- d_monthly_uk_pork$imports_tonnes * 1000
d_monthly_uk_pork$exports_kg <- d_monthly_uk_pork$exports_tonnes * 1000
d_monthly_uk_pork$total_supply_kg <- with(d_monthly_uk_pork, UK_production_kg - exports_kg + imports_kg)

state_variables <- with(d_monthly_uk_pork, 
                        cbind( breeding_herd_head, total_supply_kg, NA, all_pig_price_p_kg)
)

state_variables <- apply(state_variables, 2, function(x) ifelse( is.na(x), -100, x))

draws <- read.csv(paste0(code_path, "Stan-MCMC-results.csv"))
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
     ylim=c(min(na.omit(state_variables_plot[,1])) * 0.9 , max(na.omit(state_variables_plot[,1])) * 1.1 ),
     axes=F, xlab = "months/years", ylab = "breeding herd (head)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
# axis(side=1, at = seq(2, n_t+1, 6), 
#      labels = rep(c("Jan", "Jul"), times=5), cex.axis=0.8)
# abline(v=seq(2,n_t+1,12), lty=2, col=scales::alpha("black",0.5))
axis(side=2, at = seq(300e3, 550e3,50e3))
for(i in 1:n_samples){
  points(1:(n_t+1), 
         rlnorm(n = n_t+1, 
                meanlog = log(unlist(the_state[sample_draws[i],])), 
                sdlog = draws[sample_draws[i], "sigma.1."]),
         col = scales::alpha("slateblue",0.1)
  )
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.5), lwd=0.3)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=4)
lines(c(initial, state_variables_plot[,1]))
points(c(initial, state_variables_plot[,1]), col=scales::alpha("orange",1), pch=16)
mtext("a", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3b
png("figure_3b.png", width=1300, height=1000, res=200)
the_state <- cbind(states_init[,3], states_draws[ , ( ((3-1)*(n_t)+1):(n_t*3))])
initial <- NA
plot(rnorm(n_t+1), type="n", ylim=c(5e7, 35e7),
     axes=F, xlab = "months/years", ylab = "demand (kg)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(5e7, 50e7, 5e7))
for(i in 1:n_samples){
  points(1:(n_t+1), 
         rlnorm(n = n_t+1, 
                meanlog = log(unlist(the_state[sample_draws[i],])), 
                sdlog = draws[sample_draws[i], "sigma.3."]),
         col = scales::alpha("slateblue",0.1)
  )
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.5), lwd=0.3)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=4)
mtext("b", line = 1, at = 1, font=2, cex=1.5)
dev.off()


# Figure 3c
png("figure_3c.png", width=1300, height=1000, res=200)
the_state <- cbind(states_init[,2], states_draws[ , ( ((2-1)*(n_t)+1):(n_t*2))])
initial <- NA
plot(c(initial, state_variables_plot[,2]), type="n", 
     ylim=c(min(na.omit(state_variables_plot[,2])) * 0.9 , max(na.omit(state_variables_plot[,2])) * 1.1 ),
     axes=F, xlab = "months/years", ylab = "inventory (kg)", cex.lab=1.5
)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(10e7, 20e7, 3e7))
for(i in 1:n_samples){
  points(1:(n_t+1), 
         rlnorm(n = n_t+1, 
                meanlog = log(unlist(the_state[sample_draws[i],])), 
                sdlog = draws[sample_draws[i], "sigma.2."]),
         col = scales::alpha("slateblue",0.1)
  )
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.5), lwd=0.3)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=4)
lines(c(initial, state_variables_plot[,2]))
points(c(initial, state_variables_plot[,2]), col="orange", pch=16)
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
for(i in 1:n_samples){
  points(1:(n_t+1), 
         rlnorm(n = n_t+1, 
                meanlog = log(unlist(the_state[sample_draws[i],])), 
                sdlog = draws[sample_draws[i], "sigma.4."]),
         col = scales::alpha("slateblue",0.1)
  )
  lines(1:(n_t+1), the_state[sample_draws[i], ], 
        col = scales::alpha("black", 0.5), lwd=0.3)
}
lines(1:(n_t+1), apply(the_state, 2, mean), col = scales::alpha("black", 1), lwd=4)
lines(c(initial, state_variables_plot[,4]))
points(c(initial, state_variables_plot[,4]), col="orange", pch=16)
mtext("d", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3e
png("figure_3e.png", width=1300, height=1000, res=200)
states_bh <- cbind(states_init[,1], states_draws[ , ( ((1-1)*(n_t)+1):(n_t))])
initial <- NA
plot(c(initial, d_monthly_uk_pork$UK_production_kg), type="n",
     ylim=c(min(na.omit(d_monthly_uk_pork$UK_production_kg)) * 0.95 , max(na.omit(d_monthly_uk_pork$UK_production_kg)) * 1.05 ) 
     , col = "orange", pch=16,
     axes=F, xlab = "months/years", ylab = "UK production (kg)", cex.lab=1.5)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(65e6, 95e6, 5e6))
for(i in 1:n_samples){
  points(1:(n_t+1), 
         rlnorm(n = n_t+1, 
                meanlog = log(unlist(states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."])), 
                sdlog = draws[sample_draws[i], "sigma_production"]),
         col = scales::alpha("slateblue",0.1)
  )
  lines(1:(n_t+1), states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."], 
        col = scales::alpha("black", 0.5), lwd=0.3)
}
lines(1:(n_t+1), apply( states_bh * draws$p.4. * draws$p.5., 2, mean), 
      col = scales::alpha("black", 1), lwd=4)
lines(c(initial, d_monthly_uk_pork$UK_production_kg))
points(c(initial, d_monthly_uk_pork$UK_production_kg), col="orange", pch=16)
mtext("e", line = 1, at = 1, font=2, cex=1.5)
dev.off()

# Figure 3f
# png("figure_3f.png", width=1300, height=1000, res=200)
# initial <- NA
# plot(c(initial, d_monthly_uk_pork$exports_kg), type="n", 
#      col = scales::alpha("slateblue",0.6), pch=16,
#      ylim=c(min(na.omit(d_monthly_uk_pork$exports_kg)) * 0.9 , max(na.omit(d_monthly_uk_pork$exports_kg)) * 1.1 ),
#      axes=F, xlab = "months/years", ylab = "Exports (kg)", cex.lab=1.5)
# axis(side=1, at = seq(2, n_t+1, 6),
#      labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
#      cex.axis=0.7)
# axis(side=2, at = seq(15e6, 45e6, 5e6))
# for(i in 1:n_samples){
#   points(1:(n_t+1), 
#          rlnorm(n = n_t+1, 
#                 meanlog = log(unlist( draws[sample_draws[i], "p.6."] * states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."])), 
#                 sdlog = draws[sample_draws[i], "sigma_production"]),
#          col = scales::alpha("slateblue",0.1)
#   )
#   lines(1:(n_t+1), draws[sample_draws[i], "p.6."] * states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."], 
#         col = scales::alpha("black", 0.5), lwd=0.3)
# }
# lines(1:(n_t+1), apply( draws$p.6. * states_bh * draws$p.4. * draws$p.5., 2, mean), 
#       col = scales::alpha("black", 1), lwd=4)
# lines(c(initial, d_monthly_uk_pork$exports_kg))
# points(c(initial, d_monthly_uk_pork$exports_kg), col="orange", pch=16)
# mtext("f", line = 1, at = 1, font=2, cex=1.5)
# dev.off()

png("figure_3f.png", width=1300, height=1000, res=200)
initial <- NA
imports_exports <- d_monthly_uk_pork$imports_kg - d_monthly_uk_pork$exports_kg
states_bh <- cbind(states_init[,1], states_draws[ , ( ((1-1)*(n_t)+1):(n_t))])
plot(c(initial, imports_exports), type="n", 
     col = scales::alpha("slateblue",0.6), pch=16,
     ylim=c(min(na.omit(imports_exports)) * 0.9 , max(na.omit(imports_exports)) * 1.1 ),
     axes=F, xlab = "months/years", ylab = "imports - exports (kg)", cex.lab=1.5)
axis(side=1, at = seq(2, n_t+1, 6),
     labels = paste( rep(c("Jan", "Jul"), times=5), rep( c("'15", "'16", "'17", "'18", "'19"), each=2)),
     cex.axis=0.7)
axis(side=2, at = seq(40e6, 70e6, 5e6))
for(i in 1:n_samples){
  points(1:(n_t+1), 
         rlnorm(n = n_t+1, 
                meanlog = log(unlist( draws[sample_draws[i], "p.6."] * (draws[sample_draws[i], "p.7."] - states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."]))), 
                sdlog = draws[sample_draws[i], "sigma_production"]),
         col = scales::alpha("slateblue",0.1)
  )
  lines(1:(n_t+1), draws[sample_draws[i], "p.6."] * (draws[sample_draws[i], "p.7."] - states_bh[sample_draws[i], ] * draws[sample_draws[i], "p.4."] * draws[sample_draws[i], "p.5."]), 
        col = scales::alpha("black", 0.5), lwd=0.3)
}
lines(1:(n_t+1), apply( draws$p.6. * (draws$p.7. - states_bh * draws$p.4. * draws$p.5.), 2, mean), 
      col = scales::alpha("black", 1), lwd=4)
lines(c(initial, imports_exports))
points(c(initial, imports_exports), col="orange", pch=16)
mtext("f", line = 1, at = 1, font=2, cex=1.5)
dev.off()