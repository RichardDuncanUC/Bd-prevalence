
# function to plot thermal curve
thermal <- function(t.min, t.max, a, b) {
  xx <- seq(t.min, t.max, 0.1)
  yy <- (a*(xx - t.min)*(1 - exp(b*(xx - t.max))))^2
  yy <- yy / max(yy)
  return(list(y = yy, x = xx))
}

# chytrid
a <- thermal(4, 27, 1, 0.3)
pathogen.opt <- a$x[which.max(a$y)]
pathogen.opt
  
  
# function to draw a plot
plot.fig <- function(tmin, tmax, b, cl, ytitle1 = "", ytitle2 = "", xtitle = "", m1 = "", m2 = "") {
  para <- thermal(4, 27, 1, 0.3)
  host <- thermal(tmin, tmax, b, 0.3)
  host$y <- host$y*1
  host.opt <- host$x[which.max(host$y)]
  hp <- which(host$y > 0.3)
  amin <- min(host$x[hp])
  amax <- max(host$x[hp])

  plot(para$y ~ para$x, type = "l", bty = "l", lwd = 2, yaxt = "n", 
       ylab = "", xlab = "", ylim = c(0, 1), xlim = c(4, 30))
  lines(host$y ~ host$x, col = cl, lwd = 2)
  mtext(m1, line = 0, adj = 0)
  
  
  title(ylab = ytitle1, line = 0.5, cex.lab = 1.2)
  title(xlab = xtitle, line = 2.2, cex.lab = 1.2)
  arrows(amin, 1.15, amax, 1.15, col = "black", xpd = NA, length = 0.05,
         angle = 90, code = 3)
  points(host.opt, max(host$y), pch = 16, col = cl, cex = 1.5, xpd = NA)
  points(pathogen.opt, max(para$y), pch = 16, cex = 1.5, xpd = NA)
  
  host.overlap <- which(host$x >= 4 & host$x <= 27)
  para.overlap <- which(para$x >= tmin & para$x <= tmax)
  dif <- rep(0, length(host$x))
  dif[host.overlap] <- para$y[para.overlap] / host$y[host.overlap]
  
  plot(dif[hp][dif[hp] > 0] ~ host$x[hp][dif[hp] > 0], type = "l", lwd = 2, bty = "l", 
       ylab = "", xlab = "", xlim = c(8, 30), ylim = c(0, 2),
       yaxt = "n", col = cl)
  lines(dif[hp][dif[hp] == 0] ~ host$x[hp][dif[hp] == 0], col = "grey80", lwd = 2)
  title(ylab = ytitle2, line = 0.5, cex.lab = 1.2)
  title(xlab = xtitle, line = 2.2, cex.lab = 1.2)
  mtext(m2, line = 0, adj = 0)
  
  # host optimum
  points(host.opt, dif[hp][which.min(abs(host.opt - host$x[hp]))], col = cl, cex = 1.5, pch = 19)
  
  # mean prevalence
  mp <- mean(dif[hp]) / 2
  text(28, 1.8, paste("Mean =", round(mp, 2)), cex = 1.2)
}

cl <- c("red3", "orange2", "green4", "blue2", "mediumorchid4")

pdf("./figures/Conceptual figure 1.pdf", height = 12, width = 8)

par(mar = c(4, 4, 2, 1))
layout(matrix(1:12, ncol=2, byrow=TRUE), 
  widths=c(1,1), 
  heights=c(1, 1, 1, 1, 1, 1.5))

plot.fig(tmin = 25, tmax = 30, b = 1, cl = cl[1], ytitle1 = "Performance", ytitle2 = "Prevalence", m1 = "A", m2 = "B")
#text(-15, 2.3, "Within host species", xpd = NA, cex = 1.3)
plot.fig(tmin = 20, tmax = 30, b = 1, cl = cl[2], ytitle1 = "Performance", ytitle2 = "Prevalence", m1 = "C", m2 = "D")
plot.fig(tmin = 10, tmax = 28, b = 1, cl = cl[3], ytitle1 = "Performance", ytitle2 = "Prevalence", m1 = "E", m2 = "F")
  text(-24, 1, expression(bold("Within host")), xpd = NA, cex = 1.6, srt = 90)
  
plot.fig(tmin = 5, tmax = 27, b = 1, cl = cl[4], ytitle1 = "Performance", ytitle2 = "Prevalence", m1 = "G", m2 = "H")
plot.fig(tmin = 5, tmax = 18, b = 0.01, cl = cl[5], ytitle1 = "Performance", ytitle2 = "Prevalence",
         xtitle = "Temperature", m1 = "I", m2 = "J")

# overall outcomes
# slope versus temperature
xx <- seq(8, 30, 0.01)
yy <- 1 - xx

plot(yy ~ xx, type = "l", lwd = 2, bty = "l", yaxt = "n", xlim = c(8, 30), xlab = "", ylab = "")
  abline(h = mean(yy))
  title(ylab = "Prevalence-Temp slope", line = 0.5, cex.lab = 1.2)
  title(xlab = "Thermal optimum", line = 2.2, cex.lab = 1.2)
  text(30, mean(yy) + 5, "+ve", xpd = NA)
  text(30, mean(yy) - 5, "-ve", xpd = NA)
  mtext("K", line = 0, adj = 0)

yy <- -(xx-17)^2

plot(yy ~ xx, type = "l", lwd = 2, bty = "l", yaxt = "n", xlim = c(8, 30), xlab = "", ylab = "",
     cex.lab = 1.2)
  title(ylab = "Prevalence at thermal optimum", line = 0.5, cex.lab = 1.2)
  title(xlab = "Thermal optimum", line = 2.2, cex.lab = 1.2)
  mtext("L", line = 0, adj = 0)
  
text(-24, -65, expression(bold("Among hosts")), xpd = NA, cex = 1.6, srt = 90)
#lines(c(-24, 31), c(35, 35), xpd = NA, lwd = 3, lty = 3)

dev.off()

