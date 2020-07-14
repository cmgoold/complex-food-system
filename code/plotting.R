source("~/Dropbox/Leeds_postdoc/Papers/cfs-model/code/linear-stability-analysis-functions.R")
setwd("~/Dropbox/Leeds_postdoc/Papers/cfs-model/paper/figures")

# bifurcation diagram for capital
png("~/Documents/bifur-fig.png", width=1000, height=1000, res=200)
alpha_grid <- seq(0.1, 2, length.out = 100)
p_list <- list(beta = 0.2, gamma = 26, omega = 3, kappa = 0.5, delta = 5)
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
     xlab = expression(alpha)
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
abline(v=alpha_critical(p_list), lty=3, col=scales::alpha("black",0.7))
legend("topleft", 
       bty = "n",
       legend = c("stable", "unstable"), lwd = 3, 
       col = c(scales::alpha("slateblue", 0.8), "black"), 
       lty=c(1,2))
dev.off()

# (alpha, kappa) parameter space
png("~/Documents/kappa_alpha_beta_1.png", width=1000, height=1000, res=200)
alpha_grid <- seq(0.01, 2, length.out = 100)
kappa_grid <- seq(0.01, 0.99, length.out = 100)
plot(kappa_grid, alpha_grid, type="n", 
     xlab = expression(kappa), ylab=expression(alpha), 
     cex.lab = 1.5
     )
ac <- numeric(length(alpha_grid)) 
the_beta <- 1
for(i in 1:length(ac)){
  ac[i] <- alpha_critical(p_list = list(kappa=kappa_grid[i], beta=the_beta, gamma=26, omega=3, delta = 5))
}
polygon(
  x = c(-1, -1, 3, 3),
  y = c(-1, ac[1], ac[length(ac)], -1), 
  col = scales::alpha("steelblue", 0.8), 
  border = scales::alpha("steelblue", 0.8)
)
abline(h=1, lty=2, col=scales::alpha("black", 0.5))
text(x = 0.5, y = 0.5, labels = "unsustainable")
text(x = 0.5, y = 1.75, labels = "sustainable")
legend("topleft", bty="n", 
       legend = expression(paste(beta, " = ") ~ 1)
       )
dev.off()
