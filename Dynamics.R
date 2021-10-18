#####################################################################################################
##### CODE TO GENERATE FIGURES
#####################################################################################################

### Dependencies

source("Functions.R")

##############################################################################
#### Dynamics figures
##############################################################################

# Setting parameters

pars <- BuildPars(S = 15, Ctype = "Random", Cd = 20, C0 = 2, bC = 1, Ptype = "Random", eps_val = 0.05, n = 0, r = 0)

eta <- rep(1, times = pars$S)

rho_mult <- 0.1
rho_const <- rep(1, times = pars$S)
rho_rand <- runif(pars$S)
rho_one <- c(1, rep(0, times = pars$S - 1))
rho <- rho_rand * rho_mult

pars$R <- GetR(eta, pars)
pars$N <- GetN(rho, pars)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))

pars$eta <- eta
pars$rho <- rho
pars$dyn <- "Minimum"

# Integrating dynamics for simplification plot

ini_mult <- 0.5
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 0.5), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult
inistate[inistate < 0] <- 0.1

endtime <- 150
timestep <- 0.5
nonconst_out <- integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep)
pars_const <- pars
pars_const$dyn <- "Constant"
const_out <- integrate_dyn(inistate, pars_const, endtime = endtime, fn = dynamics, timestep = timestep)


dyn_data <- rbind(cbind(const_out, data.frame(Dynamics = "Constant")),
                  cbind(nonconst_out, data.frame(Dynamics = "Liebig")))

rdyn <- dyn_data[,1:(pars$S+1)]
rdyn$Dynamics <- dyn_data$Dynamics
rdyn$Type <- "Resources"
cdyn <- dyn_data[,-(2:(pars$S+1))] %>% mutate(Type = "Consumers")
names(cdyn) <- names(rdyn)
dyn_data <- rbind(rdyn, cdyn)

# Simplification plot

melt_dyn <- melt(dyn_data, id.vars = c("time", "Type", "Dynamics"))
plDyn <- ggplot() +
  geom_line(data = melt_dyn, aes(x = time, y = value, color = variable),
            size = 1.5, alpha = 0.75) +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + scale_y_log10() +
  facet_grid(Type~Dynamics, scales = "free_y") +
  theme(legend.position = "none") +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        legend.text=element_text(size=15))
plDyn

tiff("../figs/SIFigSimple.tiff", width = 1000, height = 800, res = 300)
plDyn
dev.off()

# Setting parameters again

pars <- BuildPars(S = 15, Ctype = "Sparse", Cd = 5, C0 = 1, bC = 0.5, Ptype = "Constant", eps_val = 0.05, n = 0, r = 0)

eta <- rep(1, times = pars$S)

rho_mult <- 1
rho_const <- rep(1, times = pars$S)
rho_rand <- runif(pars$S)
rho_one <- c(1, rep(0, times = pars$S - 1))
rho <- rho_const * rho_mult

pars$R <- GetR(eta, pars)
pars$N <- GetN(rho, pars)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))
print(pars$N)

pars$eta <- eta
pars$rho <- rho
pars$dyn <- "SubMin"

# Integrating stable dynamics

ini_mult <- 1
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult
inistate[inistate < 0] <- 0.1

endtime <- 50
timestep <- 0.1
stable_out <- cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                    data.frame(Rho = rho_mult))

# Integrating unstable dynamics

rho_mult <- 0.0001
pars$rho <- rho_mult * pars$rho
pars$N <- GetN(pars$rho, pars)
print(pars$N)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))

ini_mult <- 0.01
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult

endtime <- endtime / rho_mult / 100
timestep <- timestep / rho_mult / 100
unstable_out <- cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                      data.frame(Rho = rho_mult))

dyn_data <- rbind(stable_out, unstable_out)

cdyn_data <- cbind(dyn_data[,-(2:(pars$S+1))], data.frame(Abd = "Consumers"))
rdyn_data <- cbind(dyn_data[,-((pars$S+2):(ncol(dyn_data)-1))], data.frame(Abd = "Resources"))

names(cdyn_data) <- names(rdyn_data)
dyn_data <- rbind(cdyn_data, rdyn_data)

# Plotting the bifurcation dynamics

melt_dyn <- melt(dyn_data, id.vars = c("time", "Rho", "Abd")) %>%
  mutate(Inflow = ifelse(Rho == min(Rho), "Unstable", "Stable"))

plBifDyn <- ggplot() +
  geom_line(data = melt_dyn, aes(x = time, y = value, color = variable),
            size = 3, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + scale_y_log10() +
  facet_wrap(Abd~Inflow, scales = "free") +
  theme(legend.position = "none") +
  theme(text = element_text(size=30),
        axis.text=element_text(size=15),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size=25))
plBifDyn

tiff("../figs/FigBifDyn.tiff", width = 1250, height = 1000, res = 125)
plBifDyn
dev.off()

# Initial conditions figure

ini_mult <- 0.1
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult

endtime <- 1e3
timestep <- 1
IC_out <- cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                      data.frame(Rho = rho_mult))

unstable_out$IC <- "(A) Close to Equilibrium"
IC_out$IC <- "(B) Far from Equilibrium"

dyn_data <- rbind(unstable_out, IC_out)
dyn_data <- dyn_data[,-(2:(pars$S+1))]

# Plotting the bifurcation dynamics

melt_dyn <- melt(dyn_data, id.vars = c("time", "Rho", "IC"))

plIC <- ggplot() +
  geom_line(data = melt_dyn, aes(x = time, y = value, color = variable),
            size = 3, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + scale_y_log10() +
  facet_wrap(~IC, scales = "free") +
  theme(legend.position = "none") +
  theme(text = element_text(size=30),
        axis.text=element_text(size=15),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size=25))
plIC

tiff("../figs/SIFigIC.tiff", width = 1800, height = 800, res = 125)
plIC
dev.off()

# Setting parameters again

pars <- BuildPars(S = 15, Ctype = "Sparse", Cd = 5, C0 = 1, bC = 0.5, Ptype = "Constant", eps_val = 0.05, n = 0, r = 0)
Cunstable <- pars$C
Cstable   <- BuildC(S = 15, MatType = "Banded", Cd = 5, C0 = 1, bC = 0.5)
eta <- rep(1, times = pars$S)

rho_mult <- 1
rho_const <- rep(1, times = pars$S)
rho_rand <- runif(pars$S)
rho_one <- c(1, rep(0, times = pars$S - 1))
rho <- rho_one * rho_mult

pars$R <- GetR(eta, pars)
pars$N <- GetN(rho, pars)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))
print(pars$N)

pars$eta <- eta
pars$rho <- rho
pars$dyn <- "SubMin"

# Integrating stable dynamics

ini_mult <- 1
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult
inistate[inistate < 0] <- 0.1

endtime <- 1000
timestep <- 1
rho1_out <- cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                    data.frame(Rho = rho_mult, Network = "Sparse"))
pars$C <- Cstable
pars$R <- GetR(eta, pars)
pars$N <- GetN(rho, pars)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))
print(pars$N)

endtime <- 1000
timestep <- 1

rho1_out <- rbind(rho1_out, cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                  data.frame(Rho = rho_mult, Network = "Banded")))

# Integrating unstable dynamics
pars$C <- Cunstable
rho_mult <- 0.1
pars$rho <- rho_mult * pars$rho
pars$R <- GetR(eta, pars)
pars$N <- GetN(pars$rho, pars)
print(pars$N)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))

ini_mult <- 0.1
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult

endtime <- 1e4
timestep <- 10
rho2_out <- cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                  data.frame(Rho = rho_mult, Network = "Sparse"))
pars$C <- Cstable
pars$R <- GetR(eta, pars)
pars$N <- GetN(rho, pars)
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))
print(pars$N)

endtime <- 1e4
timestep <- 10
rho2_out <- rbind(rho2_out, cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                  data.frame(Rho = rho_mult, Network = "Banded")))

# Integrating unstable dynamics
pars$C <- Cunstable
rho_mult <- 0.01
pars$rho <- rho_mult * pars$rho
pars$R <- GetR(eta, pars)
pars$N <- GetN(pars$rho, pars)
print(pars$N)

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))

ini_mult <- 0.01
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult

endtime <- 5e3
timestep <- 10
rho3_out <- cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                  data.frame(Rho = rho_mult, Network = "Sparse"))
pars$C <- Cstable
pars$R <- GetR(eta, pars)
pars$N <- GetN(rho, pars)
inistate <- c(pars$R, pars$N)
inistate <- inistate + c(rnorm(pars$S, sd = 1), rnorm(pars$S, sd = 1)) * c(pars$R, pars$N) * ini_mult

J <- BuildJ(pars)
plotEigs(J, Title = paste("Spectrum: Rho Multiplied by", rho_mult))
print(pars$N)

endtime <- 5e5
timestep <- 5e2
rho3_out <- rbind(rho3_out, cbind(integrate_dyn(inistate, pars, endtime = endtime, fn = dynamics, timestep = timestep),
                  data.frame(Rho = rho_mult, Network = "Banded")))

# merging different resource inflows

dyn_data <- rbind(rho1_out, rho2_out, rho3_out)
dyn_data <- dyn_data[,-(2:(pars$S+1))]

# Plotting the bifurcation dynamics

melt_dyn <- melt(dyn_data, id.vars = c("time", "Rho", "Network")) %>%
  mutate(Inflow = ifelse(Rho == min(Rho), "Low Inflow",
                         ifelse(Rho == max(Rho),
                                "High Inflow", "Medium Inflow"))) %>%
  mutate(Inflow = factor(Inflow, levels = c("Low Inflow", "Medium Inflow", "High Inflow")))

plInflows <- ggplot() +
  geom_line(data = melt_dyn, aes(x = time, y = value, color = variable),
            size = 3, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  theme_bw() + scale_y_log10() +
  facet_wrap(Network~Inflow, scales = "free") +
  theme(legend.position = "none") +
  theme(text = element_text(size=30),
        axis.text=element_text(size=15),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size=25))
plInflows

tiff("../figs/FigInflows.tiff", width = 1750, height = 1000, res = 125)
plInflows
dev.off()



