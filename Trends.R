#####################################################################################################
##### CODE TO GENERATE FIGURES
#####################################################################################################

### Dependencies

source("Functions.R")

##############################################################################
#### Trend figures
##############################################################################

maxCd_mult <- 500
num_repl <- 100

S <- 15
Ctype <- "Circulant"
C0s <- 1
bCs <- 1
Ptype <- "Circulant"
eps <- 0.9
n <- 10^-seq(2.5, 0.5, length.out = 10)
r <- 1
xname <- "N"
xvar <- n

in_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                      Ptype = Ptype, eps_val = eps, n = n, r = r)
in_params$xname <- xname
in_params$xvar <- xvar

Ctype <- "Tradeoff"
Ptype <- "Constant"

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Random"
Ptype <- "Random"

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Circulant"
Ptype <- "Circulant"
r <- n
n <- 1
xname <- "R"
xvar <- r

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Tradeoff"
Ptype <- "Constant"

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Random"
Ptype <- "Random"

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Circulant"
Ptype <- "Circulant"
eps <- 10^-seq(1, 0, length.out = 10)
n <- 1
r <- 1
xname <- "Epsilon"
xvar <- eps

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Tradeoff"
Ptype <- "Constant"

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

Ctype <- "Random"
Ptype <- "Random"

new_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs,
                       Ptype = Ptype, eps_val = eps, n = n, r = r)
new_params$xname <- xname
new_params$xvar <- xvar
in_params <- rbind(in_params, new_params)

trend_raw_data <- RunSimulation(params = in_params, maxCd_mult = maxCd_mult, num_repl = num_repl)

trend_data <- trend_raw_data %>%
  mutate(C0 = as.factor(C0)) %>%
  group_by(n, r, bC, eps_val, Ctype, C0) %>%
  summarise(AvgStab = mean(EigCd), AvgCd = mean(Cd), Bound = mean(Bound), Pred = mean(Bpred),
            StDevStab = sd(EigCd), StDevCd = sd(Cd), xname = unique(xname), xvar = unique(xvar))

empty_str <- "                           "

plStabTrends <- ggplot(trend_data, aes(x = xvar, y = AvgStab, shape = Ctype, color = Ctype)) +
  geom_errorbar(aes(ymin = AvgStab - StDevStab, ymax = AvgStab + StDevStab), width = 0, size = 1) +
  geom_point(size = 3) + theme_bw() + geom_line(aes(x = xvar, y = Bound, color = Ctype),
                                                alpha = 0.75, linetype = "dashed", size = 1) +
  geom_line(aes(x = xvar, y = Pred, color = Ctype), alpha = 0.75, size = 1) +
  xlab(paste("Epsilon", empty_str, "      n   ", empty_str, "    r      ")) +
  ylab(expression(C[d])) + labs(color = "", shape = "") + scale_x_log10() +
  facet_wrap(Ctype ~ xname, scales = "free") +
  theme(text = element_text(size=20),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.position = "top",
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin=unit(c(0,0.5,0.1,0.1),"cm"))
plStabTrends

tiff("../figs/SIFigAllStabTrends.tiff", width = 1200, height = 800, res = 125)
plStabTrends
dev.off()

plTrends <- ggplot(trend_data, aes(x = xvar, y = AvgCd, shape = Ctype, color = Ctype)) +
  geom_errorbar(aes(ymin = AvgCd - StDevCd, ymax = AvgCd + StDevCd), width = 0, size = 1) +
  geom_point(size = 3) + theme_bw() + geom_line(aes(x = xvar, y = Bound, color = Ctype),
                                                alpha = 0.75, linetype = "dashed", size = 1) +
  geom_line(aes(x = xvar, y = Pred, color = Ctype), alpha = 0.75, size = 1) +
  xlab(paste("Epsilon", empty_str, "      n   ", empty_str, "    r      ")) +
  ylab(expression(C[d])) + labs(color = "", shape = "") + scale_x_log10() +
  facet_wrap(Ctype ~ xname, scales = "free") +
  theme(text = element_text(size=20),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.position = "top",
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin=unit(c(0,0.5,0.1,0.1),"cm"))
plTrends

tiff("../figs/SIFigAllTrends.tiff", width = 1200, height = 800, res = 125)
plTrends
dev.off()

n_trend_data <- trend_data %>%
  filter(xname == "N") %>%
  mutate(Cname = ifelse(Ctype == "Tradeoff", "(A) Tradeoff",
                        ifelse(Ctype == "Circulant", "(B) Circulant", "(C) Random")))

plNTrends <- ggplot(n_trend_data, aes(x = n, y = AvgStab, shape = Ctype, color = Ctype)) +
  geom_errorbar(aes(ymin = AvgStab - StDevStab, ymax = AvgStab + StDevStab), width = 0, size = 1) +
  geom_point(size = 3) + theme_bw() + geom_line(aes(x = xvar, y = Pred, color = Ctype), size = 1, alpha = 0.75) +
  ylab(expression(C[d])) + labs(color = "", shape = "") + scale_x_log10() +
  facet_wrap(~ Cname, scales = "free") +
  theme(text = element_text(size=16),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.position = "top",
        panel.spacing = unit(2, "lines"),
        plot.margin=unit(c(0,0.5,0.1,0.1),"cm"))
plNTrends

tiff("../figs/FigNTrends.tiff", width = 1200, height = 400, res = 125)
plNTrends
dev.off()

### Correlation SI figure

S <- 15
Ctype <- "Correlated"
C0 <- 1
bC <- 0.5
ps <- seq(0, 1, length.out = 15)
Ptype <- "Constant"
eps_val <- 0.9
ns <- 10^-(2:6)
r <- 1

in_params <- crossing(S = S, Ctype = Ctype, C0 = C0s, bC = bCs, p = ps,
                      Ptype = Ptype, eps_val = eps, n = ns, r = r)

out_corr_data <- RunCorrSim(in_params, maxCd_mult = 1e4, num_repl = 5)

corr_data <- out_corr_data %>% group_by(p, n) %>%
  summarise(MeanCd = mean(Cd), ImEig = mean(ImEig)) %>%
  melt(id.vars = c("p", "n")) %>%
  mutate(n = as.factor(n)) %>%
  mutate(variable = ifelse(variable == "MeanCd", "(A) Diagonal Consumption", "(B) Imaginary Eigenvalue"))

plCorr <- ggplot(corr_data, aes(x = p, y = value, color = n)) +
  geom_point(size = 4, alpha = 0.75) +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        panel.spacing = unit(2, "lines"))
plCorr

tiff("../figs/SIFigCorrTrends.tiff", width = 1200, height = 500, res = 125)
plCorr
dev.off()


### Spectrum SI Figures
pars <- BuildPars(100, "Circulant", 10, 1, 1, 0, "Circulant", 0.05, 0.001, 1)

J <- BuildJ(pars)
plotJ <- plotEigs(J)
plotJ

C <- pars$C
Ctilde <- C
diag(Ctilde) <- (1 - diag(pars$eps)) * diag(Ctilde)
P <- pars$P

A <- J[1:pars$S, 1:pars$S]
B <- J[1:pars$S, (pars$S + 1):(2*pars$S)]
Bprime <- unique(diag(J[(pars$S + 1):(2*pars$S), 1:pars$S])) * B

omega_A <- eigen(A)$values[2]
omega_B <- eigen(Bprime)$values[4]

A_vecs <- eigen(A)
omega_As <- A_vecs$values
A_vecs <- A_vecs$vectors

B_vecs <- eigen(Bprime)
B_eigs <- B_vecs$values
B_vecs <- B_vecs$vectors
omega_Bs <- c()

for(i in 1:pars$S) {
  cur_A_vec <- A_vecs[,i]
  cur_A_vec <- abs(cur_A_vec)
  cur_diffs <- abs(colSums(cur_A_vec - abs(B_vecs)))
  cur_ind <- which(min(cur_diffs) == cur_diffs)
  omega_Bs <- c(omega_Bs, B_eigs[cur_ind])
}

lambda <- 0.5 * (omega_As + sqrt(as.complex(omega_As^2 + 4 * omega_Bs)))
lambda <- c(lambda, 0.5 * (omega_As - sqrt(as.complex(omega_As^2 + 4 * omega_Bs))))

omega_Bs <- B_eigs

lambdas <- c()
for(i in 1:pars$S) {
  for(j in 1:pars$S) {
    
    lambdas <- c(lambdas, 0.5 * (omega_As[i] + sqrt(as.complex(omega_As[i]^2 + 4 * omega_Bs[j]))))
    lambdas <- c(lambdas, 0.5 * (omega_As[i] - sqrt(as.complex(omega_As[i]^2 + 4 * omega_Bs[j]))))
    
  }
}

J_eigs <- eigen(J)$values
errors <- c()

for(lambda in lambdas) {
  errors <- c(errors, min(abs(lambda - J_eigs)))
}

sort_lambdas <- data.frame(Lambda = lambdas, Error = errors) %>%
  arrange(Error) %>%
  mutate(Lambda = round(Lambda, 10)) %>%
  distinct(Lambda)
sort_eigs <- data.frame(Lambda = J_eigs) %>%
  mutate(Lambda = round(Lambda, 10)) %>%
  distinct(Lambda)
sort_eigs <- data.frame(Real = Re(sort_eigs$Lambda), Imaginary = Im(sort_eigs$Lambda))


lambdas <- sort_lambdas$Lambda[1:nrow(sort_eigs)]

plSpectrum <- ggplot(sort_eigs, aes(x = Real, y = Imaginary)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_bw() +
  geom_point(aes(x = Re(lambdas), y = Im(lambdas)),
             color = "blue", shape = 4, size = 4) +
  ggtitle("Spectrum of the Jacobian with Predictions")


tiff("../figs/SIFigSpectrum.tiff", width = 600, height = 500, res = 125)
plSpectrum
dev.off()

# small n spectrum intuition

pars <- BuildPars(100, "Tradeoff", 15, 1, 1, 0, "Constant", 0.05, 0.01, 1)

J <- BuildJ(pars)
plJ11 <- plotEigs(J, Title = "Tradeoff Interactions:\nLarge Abundance")
plJ11

pars <- BuildPars(100, "Tradeoff", 15, 1, 1, 0, "Constant", 0.05, 0.0001, 1)
J <- BuildJ(pars)

plJ12 <- plotEigs(J, Title = "Tradeoff Interactions:\nSmall Abundance")
plJ12

pars <- BuildPars(100, "Random", 15, 1, 1, 0, "Random", 0.05, 0.01, 1)
J <- BuildJ(pars)

plJ21 <- plotEigs(J, Title = "Random Interactions:\nLarge Abundance")
plJ21

pars <- BuildPars(100, "Random", 15, 1, 1, 0, "Random", 0.05, 0.0001, 1)
J <- BuildJ(pars)

plJ22 <- plotEigs(J, Title = "Random Interactions:\nSmall Abundance")
plJ22

tiff("../figs/SIFigAbdSpectra.tiff", width = 600, height = 500, res = 125)
grid.arrange(plJ11, plJ12, plJ21, plJ22, nrow = 2)
dev.off()

# Spectrum of J as a function of the spectrum of B

n <- 1
Ss <- c(3, 5, 6, 10)
num_trials <- 1000
eig_data <- data.frame(S = c(), EigB = c(), EigJ = c())

for(S in Ss) {
  for(i in 1:num_trials) {
    
    A <- matrix(runif(S^2), S, S)
    new_diag <- eigen(0.5 * (A + t(A)), only.values = TRUE)$values[1]
    diag(A) <- diag(A) - 1.001 * new_diag
    
    B <- matrix(runif(S^2, -1, 1), S, S)
    B <- 0.5 * (B + t(B))
    new_diag <- eigen(B, only.values = TRUE)$values[1]
    new_eig <- runif(1, 0.25, 2)
    diag(B) <- diag(B) - new_eig * new_diag
    
    eig_B <- max(eigen(B, only.values = TRUE)$values)
    
    J <- matrix(0, 2*S, 2*S)
    
    J[1:S, 1:S] <- A
    J[1:S, (S+1):(2*S)] <- B
    J[(S+1):(2*S), 1:S] <- diag(n, S, S)
    
    eig_J <- max(Re(eigen(J, only.values = TRUE)$values))
    
    eig_data <- rbind(eig_data, data.frame(S = S, EigB = eig_B, EigJ = eig_J))
  }
}

eig_data$S <- as.factor(eig_data$S)
plIff <- ggplot(eig_data, aes(x = EigB, y = EigJ, color = S, shape = S)) +
  geom_point(alpha = 0.5, size = 1.5) + theme_bw() +
  geom_hline(yintercept = 0, alpha = 0.5) + geom_vline(xintercept = 0, alpha = 0.5) +
  xlab("Leading Eigenvalue of B") +
  ylab("Leading Eigenvalue of J")
plIff


tiff("../figs/SIFigIff.tiff", width = 600, height = 500, res = 125)
plIff
dev.off()



