#####################################################################################################
##### CODE TO GENERATE FIGURES
#####################################################################################################

### Dependencies

source("Functions.R")

##############################################################################
#### Bound figures
##############################################################################

maxCd_mult <- 500
num_repl <- 100

S <- 15
Ctype <- "Tradeoff"
C0 <- 1
in_bCs <- seq(0.001, C0, length.out = 10)
Ptype <- "Constant"
eps <- 0.5
n <- 1
r <- 1

in_params <- crossing(S = S, Ctype = Ctype, C0 = C0, bC = in_bCs, Ptype = Ptype, eps_val = eps, n = n, r = r)

C0 <- 3
bCs <- C0 * in_bCs

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

C0 <- 5
bCs <- C0 * in_bCs

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

Ctype <- "Circulant"
Ptype <- "Circulant"

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

C0 <- 3
bCs <- C0 * in_bCs

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

C0 <- 1
bCs <- C0 * in_bCs

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

Ctype <- "Random"
Ptype <- "Random"

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

C0 <- 3
bCs <- C0 * in_bCs

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

C0 <- 5
bCs <- C0 * in_bCs

in_params <- rbind(in_params, crossing(S = S, Ctype = Ctype, C0 = C0, bC = bCs, Ptype = Ptype, eps_val = eps, n = n, r = r))

out_data <- RunSimulation(params = in_params, maxCd_mult = maxCd_mult, num_repl = num_repl)

proc_data <- out_data %>%
  filter(!is.na(Feasibility)) %>%
  mutate(CV = bC / C0 / sqrt(3)) %>%
  mutate(C0 = as.factor(C0))

mean_data <- proc_data %>%
  group_by(S, Ctype, C0, bC, Ptype, eps_val, n, r) %>%
  summarise(MeanCd = mean(EigCd), MeanPred = mean(Bpred), MeanBound = mean(Bound), CV = unique(CV)) %>%
  mutate(Cname = ifelse(Ctype == "Tradeoff", "(A) Tradeoff",
                        ifelse(Ctype == "Circulant", "(B) Circulant", "(C) Random")))

plSuffBounds <- ggplot(mean_data, aes(x = CV, y = MeanCd, color = C0, shape = C0)) +
  geom_point(size = 4) + theme_bw() + facet_grid(. ~ Cname) +
  geom_line(aes(x = CV, y = MeanPred, color = C0), size = 1.5) +
  xlab("Coefficient of Variation") + ylab(expression(C[d])) +
  labs(color = expression(C[0]), shape = expression(C[0])) +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        panel.spacing = unit(2, "lines"))
plSuffBounds

tiff("../figs/FigBounds.tiff", width = 1200, height = 400, res = 125)
plSuffBounds
dev.off()

diff_data <- proc_data %>%
  mutate(Diff = Bpred - EigCd) %>%
  mutate(CV = as.factor(round(CV, 2))) %>%
  mutate(Cname = ifelse(Ctype == "Tradeoff", "(A) Tradeoff",
                        ifelse(Ctype == "Circulant", "(B) Circulant", "(C) Random")))

plDiff <- ggplot(diff_data, aes(x = CV, y = Diff, color = C0)) +
  geom_boxplot(size = 1.25) + theme_bw() + facet_grid(.~ Cname) + geom_hline(yintercept = 0) +
  xlab("Coefficient of Variation") + ylab("Difference") + labs(color = expression(C[0])) +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(2, "lines"))
plDiff

tiff("../figs/SIFigDiffs.tiff", width = 1200, height = 400, res = 125)
plDiff
dev.off()