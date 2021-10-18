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
n <- 10^-seq(3, 0.5, length.out = 15)
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
eps <- 10^-seq(1, 0, length.out = 15)
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
  summarise(Average = mean(Cd), Bound = mean(Bound), StDev = sd(Cd), xname = unique(xname), xvar = unique(xvar))

empty_str <- "                           "

plTrends <- ggplot(trend_data, aes(x = xvar, y = Average, shape = Ctype, color = Ctype)) +
  geom_errorbar(aes(ymin = Average - StDev, ymax = Average + StDev), width = 0, size = 1) +
  geom_point(size = 3) + theme_bw() + geom_line(aes(x = xvar, y = Bound, color = Ctype)) +
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

plNTrends <- ggplot(n_trend_data, aes(x = n, y = Average, shape = Ctype, color = Ctype)) +
  geom_errorbar(aes(ymin = Average - StDev, ymax = Average + StDev), width = 0, size = 1) +
  geom_point(size = 3) + theme_bw() + geom_line(aes(x = xvar, y = Bound, color = Ctype)) +
  ylab(expression(C[d])) + labs(color = "", shape = "") + scale_x_log10() +
  facet_wrap(~ Cname, scales = "free") +
  theme(text = element_text(size=20),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.position = "top",
        panel.spacing = unit(2, "lines"),
        plot.margin=unit(c(0,0.5,0.1,0.1),"cm"))
plNTrends

tiff("../figs/FigNTrends.tiff", width = 2500, height = 1000, res = 300)
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






