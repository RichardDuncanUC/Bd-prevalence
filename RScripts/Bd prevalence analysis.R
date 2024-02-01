
library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(jagsUI)
library(ggrepel)

# clear workspace
rm(list = ls())

# data files required to run the analysis:
# Tidy Bd prevalence data.csv
# Australia shape vector.gpkg (to draw the map of Australia)

# logit, unlogit and scale functions
logit <- function(x) log(x / (1-x))
unlogit <- function(x) exp(x) / (1 + exp(x))
sc <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)

# read in shapefile of Australian mainland
mainvec <- vect("./data/Australia shape vector.gpkg")

#--------------------------------------------------------------------------------------------------------------
# read in tidy chytrid data
all.dat <- read.csv("./data/Tidy Bd prevalence data.csv", strip.white = T)

# check for duplicates
table(duplicated(cbind(all.dat$species, all.dat$loc, all.dat$year)))

# total number of individuals per species and number of locations
all.dat <- all.dat |>
  group_by(species) |>
  mutate(total.ind = sum(n.ind),
         total.loc = n()) |>
  glimpse()
  
################################################################################
# first year infection recorded
first.year <- min(all.dat$year[all.dat$n.pos > 0])
first.year

# year range
range(all.dat$year)

# infection prevalence over time
y <- all.dat |>
  group_by(year) |>
  summarise(nind = sum(n.ind),
            npos = sum(n.pos)) |>
  mutate(prop = npos/nind) |>
  glimpse()

# prevalence per species
ps <- all.dat |>
  mutate(prev = n.pos / n.ind) |>
  group_by(species) |>
  summarise(mean.prev = mean(prev),
            n = n()) |>
  glimpse() 

#-------------------------------------------------------------------------------
# map of infections, temporal pattern and prevalence distribution
pdf("./Figures/Figure 2.pdf")

par(mfrow = c(2, 2))
plot(mainvec, ylab = "Latitude", xlab = "Longitude")
points(lat ~ long, data = filter(all.dat, n.pos == 0), pch = 16, col = rgb(0, 0, 1, 0.3), cex = 1)
points(lat ~ long, data = filter(all.dat, n.pos > 0), pch = 16, col = rgb(1, 0, 0, 0.3), cex = 1)
mtext("A", adj = 0)

# mean prevalence per species
par(mar = c(4, 4, 1, 1))
hist(ps$mean.prev, breaks = seq(0, 1, 0.1), xlab = "Mean prevalence", main = "",
     cex.lab = 1.2)
text(0.8, 12, "Species, n = 42", cex = 1.2, xpd = NA)
mtext("C", adj = 0, xpd = NA)

# infection prevalence over time
# scale size proportional to log number
ss <- log(y$nind) / 4
plot(prop ~ year, pch = 19, data = y, cex = ss, bty = "l", 
     ylab = "Mean prevalence", xlab = "Year", cex.lab = 1.2)

p <- c(10, 100, 1000)
lp <- log(p) / 4
legend(1980, 1, legend = c("10", "100", "1000"), pch = 19, pt.cex = lp, cex = 0.8,
       title = "Number of individuals")
mtext("B", adj = 0, xpd = NA)

# prevalence per location
hist(all.dat$n.pos / all.dat$n.ind, breaks = seq(0, 1, 0.1), xlab = "Prevalence", main = "",
     cex.lab = 1.2)
text(0.8, 800, "Samples, n = 1455", cex = 1.2, xpd = NA)
mtext("D", adj = 0, xpd = NA)

dev.off()

################################################################################
# use the model that allows for annual temperature fluctuations
all.dat$temp <- all.dat$temp.year

#-------------------------------------------------------------------------------
# summaries of data

# number of individuals sampled
sum(all.dat$n.ind)

# number of positives
sum(all.dat$n.pos)
sum(all.dat$n.pos) / sum(all.dat$n.ind)

# number of species
length(table(all.dat$species))

# number of sampling events (species x location x year)
nrow(all.dat)

# number of locations
length(table(all.dat$loc))

sum.dat <- all.dat |>
  group_by(species) |>
  summarise(n.loc = n(),
            n.ind = sum(n.ind)) |>
  arrange(-n.ind)

sum.dat

#-------------------------------------------------------------------------------
# extract data for JAGS model

y <- all.dat$n.pos                           # number of Bd positive individuals in a species-by-location-by-year sample
n.ind <- all.dat$n.ind                       # total number of individuals tested in a sample
spp <- as.numeric(factor(all.dat$species))   # indicator for species
loc <- as.numeric(factor(all.dat$loc))       # indicator for sample location

temp <- all.dat$temp                         # mean annual temperature per location-year
temp[is.na(temp)] <- mean(temp, na.rm = T)   # set missing values to the mean

# centre and scale temperature
sc.temp <- sc(temp)

N <- length(y)
N.spp <- max(spp)
N.loc <- max(loc)

# continuous value for year scaled
year <- sc(all.dat$year)
table(year)

# indicator for year
fac.year <- as.numeric(factor(year))
N.year <- max(fac.year)

# observation level random effect to account for overdispersion
od <- 1:length(y)

#-------------------------------------------------------------------------------
# get thermal tolerance data
mt <- read.csv("./data/Thermal optima.csv") |>
  glimpse()

# add decline status
dec <- all.dat |>
  dplyr::select(species, decline) |>
  unique()

mt <- left_join(mt, dec)

# species-level median temp
mtemp <- mt$med.temp

# put on same scale as temp
sc.mtemp <- (mtemp - mean(temp)) / sd(temp)

################################################################################
# model with b.temp slope as a function of mean temp

mod <- "model {
  for(i in 1:N) {
    # likelihood for prevalence
    y[i] ~ dbin(p[i], n.ind[i])                                                   # eq 1
    logit(p[i]) <- int + b.temp[spp[i]]*temp[i] + b.year.trend*year[i] +
                   b.year[fac.year[i]] + b.spp[spp[i]] + b.loc[loc[i]] + b.od[i]  # eq 2
    
    # overdispersion
    b.od[i] ~ dnorm(0, tau.od)                                                    # eq 8
  }

  for(i in 1:N.spp) {
    mtemp[i] ~ dnorm(0, 0.01)
    b.spp[i] ~ dnorm(0, tau.spp)                                                  # eq 3
    b.temp[i] ~ dnorm(mu.temp[i], tau.temp)                                       # eq 4
    mu.temp[i] <- int.temp + slope.temp * mtemp[i]                                # eq 5
  }

  for(i in 1:N.loc) {
    b.loc[i] ~ dnorm(0, tau.loc)                                                  # eq 6
  }
  
  for(i in 1:N.year) {
    b.year[i] ~ dnorm(0, tau.year)                                                # eq 7
  }

  # priors
  int ~ dnorm(0, 0.01)
  int.temp ~ dnorm(0, 0.01)
  slope.temp ~ dnorm(0, 0.01)
  b.year.trend ~ dnorm(0, 0.01)

  tau.spp <- 1/(sigma.spp * sigma.spp)
  sigma.spp ~ dunif(0, 10)

  tau.loc <- 1/(sigma.loc * sigma.loc)
  sigma.loc ~ dunif(0, 10)
  
  tau.year <- 1 / (sigma.year * sigma.year)
  sigma.year ~ dunif(0, 10)

  tau.od <- 1 / (sigma.od * sigma.od)
  sigma.od ~ dunif(0, 10)

  tau.temp <- 1 / (sigma.temp * sigma.temp)
  sigma.temp ~ dunif(0, 10)

}"  

write(mod, "model.txt")

# Sample model in JAGS with three chains
set.seed(234)
mod.jags <- jags(model = "model.txt",
                  data = list(y = y, n.ind = n.ind, N = N, spp = spp, loc = loc, N.spp = N.spp, 
                              N.loc = N.loc, year = year, fac.year = fac.year, N.year = N.year,
                              temp = sc.temp, mtemp = sc.mtemp),
                  #              inits = function() list(), 
                  param = c("int", "int.temp", "slope.temp", "b.year.trend", "b.year", 
                            "sigma.od", "sigma.year", "sigma.spp", "sigma.loc", "sigma.temp", 
                            "b.temp", "b.spp","b.loc"),
                  n.chains = 3,
                  n.iter = 12000,
                  n.burnin = 2000,
                  parallel = T)

jags.sum <- mod.jags$summary
jags.sum[1:40, c(1, 2, 3, 7, 8, 9)]

tail(jags.sum)[, 1:9]

#-------------------------------------------------------------------------------
# standard deviation parameters
sigma.spp <- mod.jags$sims.list$sigma.spp
sigma.loc <- mod.jags$sims.list$sigma.loc
sigma.year <- mod.jags$sims.list$sigma.year
sigma.od <- mod.jags$sims.list$sigma.od
sigma.temp <- mod.jags$sims.list$sigma.temp

median(sigma.spp)
median(sigma.loc)
median(sigma.year)
median(sigma.od)
median(sigma.temp)

#-------------------------------------------------------------------------------
# estmate mean prevalence for each species as a function of temperature
# extract overall intercept, species random effect and species slopes
overall.int <- mod.jags$sims.list$int
spp.int <- mod.jags$sims.list$b.spp
slope <- mod.jags$sims.list$b.temp
b.year <- mod.jags$sims.list$b.year.trend

# predicted prevalence for each species
# do this for the year 2006, which is the median year across all individuals tested
median(rep(all.dat$year, all.dat$n.ind))

# convert to standardised year value
yr <- unique(year[all.dat$year == 2006])
# range of standardised temperature values to model
tt <- seq(-3, 2.5, 0.01)

# p is prevalence on the logit scale
# prev is prevalence on the probability scale
# m.prev is mean prevalence across the species temperature range
# opt.prev is prevalence at each species thermal optimum
p <- prev <- array(dim = c(ncol(spp.int), length(tt), nrow(spp.int)))
m.prev <- opt.prev <- matrix(nrow = nrow(spp.int), ncol = ncol(spp.int))

# convert tt to actual temperatures
s <- sd(temp)
m <- mean(temp)
tmp <- tt*s + m

# calculation for each draw from the posterior (rows)
dim(spp.int)
for(i in 1:nrow(spp.int)) {
  for(j in 1:ncol(spp.int)) {
    p[j, ,i] <- overall.int[i] + spp.int[i, j] + slope[i, j] * tt + b.year[i] * yr 
    prev[j, ,i] <- exp(p[j, ,i]) / (1 + exp(p[j, ,i]))
    m.prev[i, j] <- mean(prev[j, ,i][tmp > mt$min.temp[j] & tmp < mt$max.temp[j]])
    opt.temp <- which.min(abs(tmp - mt$med.temp[j]))
    opt.prev[i, j] <- prev[j, ,i][opt.temp]
  }
}

#-------------------------------------------------------------------------------
# draw figure using the thermal range for each species

pdf("./Figures/Figure 3.pdf")

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

cl <- ifelse(mt$decline == "yes", rgb(1, 0, 0, 0.7), rgb(0.1, 0.1, 0.1, 0.7))
cx <- sqrt(ps$n/20)

# sens is prevalence at mean temp
sens <- numeric()

# mean prevalence across the thermal range
mean.prev <- numeric()

# for first species
stmp <- tmp[tmp >= mt$min.temp[1] & tmp <= mt$max.temp[1]]
sprev <- apply(prev[1, , ], 1, median)[tmp > mt$min.temp[1] & tmp < mt$max.temp[1]]
mean.prev[1] <- mean(sprev)
plot(sprev ~ stmp, type = "l", ylim = c(0, 1), col = cl[i], xlim = c(5, 30),
     bty = "l", xlab = "", ylab = "")
sens[1] <- sprev[which.min(abs(mt$med.temp[1] - stmp))]
points(mt$med.temp[1], sens[1], pch = 16, col = cl[1])
title(ylab = "Prevalence", xlab = "Mean annual temperature (\u00B0C)", line = 2.5)
mtext("A", adj = 0)

# subsequent species
for(i in 2:nrow(mt)) {
  stmp <- tmp[tmp >= mt$min.temp[i] & tmp <= mt$max.temp[i]]
  sprev <- apply(prev[i, , ], 1, median)[tmp > mt$min.temp[i] & tmp < mt$max.temp[i]]
  mean.prev[i] <- mean(sprev)
  lines(sprev ~ stmp, col = cl[i])
  sens[i] <- sprev[which.min(abs(mt$med.temp[i] - stmp))]
  points(mt$med.temp[i], sens[i], pch = 16, col = cl[i])
}

# slope as a function of mean temperature
int.temp <- median(mod.jags$sims.list$int.temp)
slope.temp <- median(mod.jags$sims.list$slope.temp)
plot(apply(mod.jags$sims.list$b.temp, 2, median) ~ sc.mtemp, pch = 16, col = cl, bty = "l",
     ylab = "", xlab = "", xaxt = "n", cex = cx)
abline(int.temp, slope.temp)
abline(h = 0, lty = 3)
# add x axis
locx <- seq(10, 25, 5)
loc.sc <- (locx - m) / s
axis(1, at = loc.sc, labels = locx)
title(ylab = "Prevalence-Temperature slope", xlab = "Thermal optimum(\u00B0C)", line = 2.5)
mtext("B", adj = 0)

dev.off()

#-------------------------------------------------------------------------------
# fit models for relationship between prevalence and thermal optimum and propagate uncertainty
nsim <- 100
cf <- cf.opt <- array(dim = c(nrow(m.prev), 5, nsim))
med.temp <- mt$med.temp
for(i in 1:nrow(m.prev)) {
  # mean prevalence
  m0 <- lm(logit(m.prev[i, ]) ~ med.temp)
  cf[i, 1, ] <- rnorm(nsim, mean = coef(m0)[1], sd = sqrt(diag(vcov(m0)))[1])
  cf[i, 2, ] <- rnorm(nsim, mean = coef(m0)[2], sd = sqrt(diag(vcov(m0)))[2])
  m1 <- lm(logit(m.prev[i, ]) ~ med.temp + I(med.temp^2))
  cf[i, 3, ] <- rnorm(nsim, mean = coef(m1)[1], sd = sqrt(diag(vcov(m1)))[1])
  cf[i, 4, ] <- rnorm(nsim, mean = coef(m1)[2], sd = sqrt(diag(vcov(m1)))[2])
  cf[i, 5, ] <- rnorm(nsim, mean = coef(m1)[3], sd = sqrt(diag(vcov(m1)))[3])
  # prevalence at optimum
  m0 <- lm(logit(opt.prev[i, ]) ~ med.temp)
  cf.opt[i, 1, ] <- rnorm(nsim, mean = coef(m0)[1], sd = sqrt(diag(vcov(m0)))[1])
  cf.opt[i, 2, ] <- rnorm(nsim, mean = coef(m0)[2], sd = sqrt(diag(vcov(m0)))[2])
  m1 <- lm(logit(opt.prev[i, ]) ~ med.temp + I(med.temp^2))
  cf.opt[i, 3, ] <- rnorm(nsim, mean = coef(m1)[1], sd = sqrt(diag(vcov(m1)))[1])
  cf.opt[i, 4, ] <- rnorm(nsim, mean = coef(m1)[2], sd = sqrt(diag(vcov(m1)))[2])
  cf.opt[i, 5, ] <- rnorm(nsim, mean = coef(m1)[3], sd = sqrt(diag(vcov(m1)))[3])
}

#-------------------------------------------------------------------------------
# plot Bd prevalence versus thermal optima

pdf("./Figures/Figure 5.pdf")

par(mfrow = c(2, 2))
cl <- ifelse(mt$decline == "yes", rgb(1, 0, 0, 0.7), rgb(0.1, 0.1, 0.1, 0.7))
cx <- sqrt(ps$n/20)

# prevalence at thermal optimum versus thermal optimum
# fitted relationship
xx <- seq(7.3, 26, 0.1)
obs <- logit(colMeans(opt.prev))
yy <- median(cf.opt[, 3, ]) + median(cf.opt[, 4, ])*xx + median(cf.opt[, 5, ])*xx^2
pred <- median(cf.opt[, 3, ]) + median(cf.opt[, 4, ])*mt$med.temp + median(cf.opt[, 5, ])*mt$med.temp^2
rss <- sum((obs - pred)^2)
ss <- sum((obs - mean(obs))^2)
r2 <- round(1 - rss/ss, 2)
eq <- substitute(italic(r)^2~"="~r2, 
         list(r2 = r2))
plot(logit(colMeans(opt.prev)) ~ mt$med.temp, pch = 19, col = cl, cex = cx, bty = "l", yaxt = "n",
     ylim = c(-6.91, 1.39), ylab = "Predicted prevalence at thermal optimum", xlab = "Thermal optimum(\u00B0C)",
     main = eq)
axis(2, at = logit(c(0.001, 0.01, 0.1, 0.4, 0.8)), labels = c("0.001", "0.01", "0.1", "0.4", "0.8"))
lines(yy ~ xx, lwd = 2, col = "blue")
mtext("A", adj = 0)

# mean prevalence versus thermal optimum
# fitted relationship
obs <- logit(colMeans(m.prev))
yy <- median(cf[, 3, ]) + median(cf[, 4, ])*xx + median(cf[, 5, ])*xx^2
pred <- median(cf[, 3, ]) + median(cf[, 4, ])*mt$med.temp + median(cf[, 5, ])*mt$med.temp^2
rss <- sum((obs - pred)^2)
ss <- sum((obs - mean(obs))^2)
r2 <- round(1 - rss/ss, 2)
eq <- substitute(italic(r)^2~"="~r2, 
         list(r2 = r2))
plot(logit(colMeans(m.prev)) ~ mt$med.temp, pch = 19, col = cl, cex = cx, bty = "l", yaxt = "n",
     ylim = c(-6.91, 1.39), ylab = "Mean predicted prevalence", xlab = "Thermal optimum(\u00B0C)",
     main = eq)
axis(2, at = logit(c(0.001, 0.01, 0.1, 0.4, 0.8)), labels = c("0.001", "0.01", "0.1", "0.4", "0.8"))
lines(yy ~ xx, lwd = 2, col = "blue")
mtext("B", adj = 0)

# mean observed prevalence versus thermal optimum
summary(m1 <- lm(logit(ps$mean.prev + 0.001) ~ mt$med.temp + I(mt$med.temp^2)))
eq <- substitute(italic(r)^2~"="~r2, 
         list(r2 = format(summary(m1)$r.squared, digits = 2)))
plot(logit(ps$mean.prev + 0.001) ~ mt$med.temp, pch = 19, col = cl, cex = cx, bty = "l", yaxt = "n",
     ylim = c(-6.91, 1.39), ylab = "Mean observed prevalence", xlab = "Thermal optimum(\u00B0C)",
     main = eq)
axis(2, at = logit(c(0.001, 0.01, 0.1, 0.4, 0.8)), labels = c("0.001", "0.01", "0.1", "0.4", "0.8"))
yy <- coef(m1)[1] + coef(m1)[2]*xx + coef(m1)[3]*xx^2
lines(yy ~ xx, lwd = 2, col = "blue")
mtext("C", adj = 0)

# mean observed prevalence versus mean predicted prevalence
summary(m1 <- lm(ps$mean.prev ~ colMeans(m.prev)))
eq <- substitute(italic(r)^2~"="~r2, 
         list(r2 = format(summary(m1)$r.squared, digits = 2)))
plot(ps$mean.prev ~ colMeans(m.prev), pch = 19, col = cl, cex = cx, bty = "l",
     ylab = "Mean observed prevalence", xlab = "Mean predicted prevalence", main = eq, yaxs = "i", xaxs = "i",
     xpd = NA, ylim = c(0, 0.8), xlim = c(0, 0.8))
abline(0, 1, lwd = 2, col = "blue")
mtext("D", adj = 0)

dev.off()

#-------------------------------------------------------------------------------
# plot parameter estimates
# function to scale axes for plot
pl <- function(x) {
  sc <<- max(abs(quantile(x, c(0.01, 0.99))))
  mid <- median(x)
  ucl <- quantile(x, 0.975)
  lcl <- quantile(x, 0.025)
  ucl50 <- quantile(x, 0.75)
  lcl50 <- quantile(x, 0.25)
  plot(mid, ylim = c(-sc, sc), yaxs = "i", pch = 19, cex = 2, xaxt = "n", ylab = "", xlab = "", bty = "l",
       cex.lab = 1.2)
  arrows(1, ucl, 1, lcl, length = 0)
  arrows(1, ucl50, 1, lcl50, length = 0, lwd = 3)
  abline(h = 0)
}

pdf("./Figures/Figure 4.pdf")

# Year trend
par(mfrow = c(1, 4), mar = c(6, 4, 8, 1))
pl(mod.jags$sims.list$b.year.trend)
  title(ylab = "Parameter estimate", line = 2.5, cex.lab = 1.2)
  text(1, sc + sc/3, "Year trend", xpd = NA, cex = 1.2)
  mtext("A", adj = -0.2)

# Prevalence-temperature slope
pl(mod.jags$sims.list$slope.temp)
  text(1, sc + sc/3, "Prevalence-temp slope vs", xpd = NA, cex = 1.2)
  text(1, sc + sc/4, "thermal optimum", xpd = NA, cex = 1.2)
  mtext("B", adj = -0.2)

# linear and quadratic terms for predicted prevalence at thermal optimum versus thermal optimum
pl(cf.opt[, 4, ])
  text(1, sc + sc/8, "Linear", xpd = NA, cex = 1.2)
  lines(c(0.6, 2.7), c(sc + sc/5, sc + sc/5), xpd = NA, lty = 2)
  mtext("C", adj = -0.2)

pl(cf.opt[, 5, ])
  text(0.3, sc + sc/3, "Prevalence vs", xpd = NA, cex = 1.2)
  text(0.3, sc + sc/4, "thermal optimum", xpd = NA, cex = 1.2)
  text(1, sc + sc/8, "Quadratic", xpd = NA, cex = 1.2)
  mtext("D", adj = -0.2)

dev.off()

#-------------------------------------------------------------------------------
# results for mean and observed prevalence
pdf("./Figures/Suppl Figure prevalence parameters.pdf")

par(mfrow = c(1, 4), mar = c(6, 4, 8, 1))

# Mean prevalence
pl(cf[, 4, ])
title(ylab = "Parameter estimate", line = 2.5, cex.lab = 1.2)
text(1, sc + sc/8, "Linear", xpd = NA, cex = 1.2)
lines(c(0.6, 2.7), c(sc + sc/5, sc + sc/5), xpd = NA, lty = 2)
mtext("A", adj = -0.2)

pl(cf[, 5, ])
text(0.3, sc + sc/3, "Mean predicted prevalence vs", xpd = NA, cex = 1.2)
text(0.3, sc + sc/4, "thermal optimum", xpd = NA, cex = 1.2)
text(1, sc + sc/8, "Quadratic", xpd = NA, cex = 1.2)
mtext("B", adj = -0.2)

# Observed prevalence
m1 <- lm(logit(ps$mean.prev + 0.001) ~ mt$med.temp + I(mt$med.temp^2))
ci <- confint(m1)
sc <- max(abs(ci[2, ]))
plot(coef(m1)[2], ylim = c(-sc, sc), yaxs = "i", pch = 19, cex = 2, xaxt = "n", ylab = "", xlab = "", bty = "l",
     cex.lab = 1.2)
arrows(1, ci[2, 2], 1, ci[2, 1], length = 0)
abline(h = 0)
text(1, sc + sc/8, "Linear", xpd = NA, cex = 1.2)
lines(c(0.6, 2.7), c(sc + sc/5, sc + sc/5), xpd = NA, lty = 2)
mtext("C", adj = -0.2)

sc <- max(abs(ci[3, ]))
plot(coef(m1)[3], ylim = c(-sc, sc), yaxs = "i", pch = 19, cex = 2, xaxt = "n", ylab = "", xlab = "", bty = "l",
     cex.lab = 1.2)
arrows(1, ci[3, 2], 1, ci[3, 1], length = 0)
abline(h = 0)
text(0.3, sc + sc/3, "Mean observed prevalence vs", xpd = NA, cex = 1.2)
text(0.3, sc + sc/4, "thermal optimum", xpd = NA, cex = 1.2)
text(1, sc + sc/8, "Quadratic", xpd = NA, cex = 1.2)
mtext("D", adj = -0.2)

dev.off()

#-------------------------------------------------------------------------------
# labelled observed prevalence per species with predicted mean prevalence

ps$pred.prev <- colMeans(m.prev)

ggplot(ps, aes(y = mean.prev, x = pred.prev)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = species), max.overlaps = 50) +
  geom_abline() +
  theme_classic()

#-------------------------------------------------------------------------------
# linking prevalence to decline

dec <- full_join(mt, ps) |>
  glimpse()

summary(m1 <- lm(logit(pred.prev) ~ decline, data = dec))
summary(m2 <- lm(logit(mean.prev + 0.001) ~ decline, data = dec))

pdf("./Figures/Figure 6.pdf")

par(mfrow = c(1, 2), mar = c(6, 4, 8, 1))
mtt <- paste0(round(m1$coef[2], 3), "  (", round(confint(m1)[2, 1], 3), " - ",
              round(confint(m1)[2, 2], 3), ")")
plot(logit(pred.prev) ~ factor(decline), data = dec, col = "white",
     xlab = "Decline", ylab = "Mean predicted prevalence", ylim = c(-6.91, 1.39),
     yaxt = "n")
text(1.5, 3, "Difference in mean", xpd = NA)
text(1.5, 2.5, "on logit scale", xpd = NA)
text(1.5, 2, mtt, xpd = NA)
axis(2, at = logit(c(0.001, 0.01, 0.1, 0.4, 0.8)), labels = c("0.001", "0.01", "0.1", "0.4", "0.8"))
points(logit(pred.prev) ~ ifelse(decline == "no", 1, 2), data = dec, 
       pch = 19, col = cl, cex = cx)
mtext("A", adj = -0.2)

mtt <- paste0(round(m2$coef[2], 3), "  (", round(confint(m2)[2, 1], 3), " - ",
              round(confint(m2)[2, 2], 3), ")")
plot(logit(mean.prev + 0.001) ~ factor(decline), data = dec, col = "white",
     xlab = "Decline", ylab = "Mean observed prevalence", ylim = c(-6.91, 1.39),
     yaxt = "n")
text(1.5, 3, "Difference in mean", xpd = NA)
text(1.5, 2.5, "on logit scale", xpd = NA)
text(1.5, 2, mtt, xpd = NA)
axis(2, at = logit(c(0.001, 0.01, 0.1, 0.4, 0.8)), labels = c("0.001", "0.01", "0.1", "0.4", "0.8"))
points(logit(mean.prev + 0.001) ~ ifelse(decline == "no", 1, 2), data = dec, 
       pch = 19, col = cl, cex = cx)
mtext("B", adj = -0.2)

dev.off()
