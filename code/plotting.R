source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/linear-stability-analysis-functions.R")
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
chain_1 <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/R/chain_1_final.csv")
chain_2 <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/R/chain_2_final.csv")
chain_3 <- read.csv("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/R/chain_3_final.csv")

chains <- rbind(chain_1, chain_2, chain_3)

# Figure 3a
plot(d_uk_pork_fit$breeding_herd_head, type="b", ylim=c(250e3, 550e3), 
     axes=F, ylab="breeding sows (head)", xlab = "years", cex.lab=1.5)
axis(side=1, at=c(1:5), labels = c(2015, 2016, 2017, 2018, 2019))
axis(side=2, at=seq(100e3, 600e3, 100e3))
mtext("a", font=2, line=0, at = 1, cex=1.5)
weeks_year <- seq(1, N*1/dt, 1/dt)
for(i in sample(1:nrow(chains), 100, replace=F)){
  res <- cidp_model(as.matrix(chains[i,]))[weeks_year,1]
  lines( aggregate(res ~ rep(1:n_years, each=52), FUN = mean)[,2], col=scales::alpha("blue",0.8), lwd=0.2 )
}
lines( aggregate(cidp_model( as.matrix( apply(chains, 2, mean) ) )[weeks_year,1] ~ 
                   rep(1:n_years, each=52), FUN = mean)[,2], lwd=3)

# Figure 3b
plot(d_uk_pork_fit$pig_price_pence_kg, type="b", ylim=c(100,200), 
     axes=F, xlab = "years", ylab="price (p/kg)", cex.lab=1.5)
axis(side=1, at=c(1:5), labels = c(2015, 2016, 2017, 2018, 2019))
axis(side=2, at=seq(100, 300, 50))
mtext("b", font=2, line=0, at = 1, cex=1.5)
for(i in sample(1:nrow(chains), 100, replace=T)){
  res <- cidp_model(as.matrix(chains[i,]))[weeks_year,4]
  lines( aggregate(res ~ rep(1:n_years, each=52), FUN = mean)[,2], col=scales::alpha("blue",0.8), lwd=0.2 )
}
lines( aggregate(cidp_model( as.matrix( apply(chains, 2, mean) ) )[weeks_year,4] ~ 
                   rep(1:n_years, each=52), FUN = mean)[,2], lwd=3)

# Figure 3c
plot(d_uk_weekly_price[-1], type="l", ylim=c(100,200),
     axes=F, xlab = "weeks", ylab="price (p/kg)", cex.lab=1.5)
for(i in sample(1:nrow(out), 100, replace=T)){
  res <- cidp_model(as.matrix(chains[i,]))[weeks_year,4]
  lines( res, col=scales::alpha("blue",0.8), lwd=0.2 )
}
lines(cidp_model(apply(out, 3, mean))[,4], lwd=3)