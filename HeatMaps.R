#####################################################################################################
##### CODE TO GENERATE FIGURES
#####################################################################################################

### Dependencies

source("Functions.R")

##############################################################################
#### Heat map figures
##############################################################################

S <- 15

Ctype <- "Random"
Ptype <- "Random"

C0s <- 1
bCs <- 0.5
eps <- 0.05

Cds <- seq(1, 20, length.out = 25)
rho_mults <- 10^-seq(-1, 6, length.out = 25)

in_params <- crossing(S = S, Ctype = Ctype, Cd = Cds, C0 = C0s, bC = bCs,
                      Ptype = Ptype, eps_val = eps, rho_mult = rho_mults)
new_params <- in_params
new_params$Ctype <- "Non-Sym Tradeoff"
new_params$Ptype <- "Constant"

in_params <- rbind(in_params, new_params)

new_params$Ctype <- "Banded"
new_params$Ptype <- "Constant"

in_params <- rbind(in_params, new_params)

new_params$Ctype <- "Lower Triangular"
new_params$Ptype <- "Upper Triangular"

in_params <- rbind(in_params, new_params)

new_params$Ctype <- "Sparse"
new_params$Ptype <- "Constant"

in_params <- rbind(in_params, new_params)

new_params$Ctype <- "Correlated"
new_params$Ptype <- "Random"

in_params <- rbind(in_params, new_params)

new_params$Ctype <- "Tradeoff"
new_params$Ptype <- "Constant"

in_params <- rbind(in_params, new_params)

out_data <- RunHeatMap(in_params = in_params, num_repl = 5)

hm_data <- out_data %>%
  group_by(Ctype, Cd, Ptype, rho_mult, Rho) %>%
  summarise(ProbFeasible = mean(Feas),
            ProbStable = mean(Stab),
            ProbLLFeas = mean(LLFeas),
            ProbLLSuccess = mean(LLSuccess),
            ProbSubSuccess = mean(SubSuccess))

plLLHeatMap <- ggplot(hm_data, aes(x = rho_mult, y = Cd, fill = ProbLLSuccess)) +
  geom_tile() + scale_fill_viridis(end = 0.9, na.value="darkred") +
  theme(strip.background = element_rect(fill="lightgray", color = "black"),
        aspect.ratio = 1,
        axis.text = element_text( size = 10 ),
        panel.background = element_rect(fill = NA),
        text = element_text(size=20),
        legend.key.size = unit(3, "cm"),
        legend.key.width = unit(1,"cm")) +
  facet_grid(Rho~Ctype) +
  scale_x_log10(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression("Resource Inflow:"~rho),
       y = expression("Diagonal Consumption Coefficient:"~C[d]),
       fill = expression("Probability")) +
  ggtitle("Feasible and stable equilibria under Liebig's law")
show(plLLHeatMap)

tiff("../figs/SIFigLLHeatMap.tiff", width = 2000, height = 1000, res = 125)
plLLHeatMap
dev.off()

plSubHeatMap <- ggplot(hm_data, aes(x = rho_mult, y = Cd, fill = ProbSubSuccess)) +
  geom_tile() + scale_fill_viridis(end = 0.9, na.value="darkred") +
  theme(strip.background = element_rect(fill="lightgray", color = "black"),
        aspect.ratio = 1,
        axis.text = element_text( size = 10 ),
        panel.background = element_rect(fill = NA),
        text = element_text(size=20),
        legend.key.size = unit(3, "cm"),
        legend.key.width = unit(1,"cm")) +
  facet_grid(Rho~Ctype) +
  scale_x_log10(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression("Resource Inflow:"~rho),
       y = expression("Diagonal Consumption Coefficient:"~C[d]),
       fill = expression("Probability")) +
  ggtitle("Feasible and stable equilibria under structured Liebig's law")
show(plSubHeatMap)

tiff("../figs/SIFigSubHeatMap.tiff", width = 2000, height = 1000, res = 125)
plSubHeatMap
dev.off()

new_params <- in_params %>%
  filter(Ctype %in% c("Banded", "Random", "Tradeoff"))

out_data <- RunHeatMap(in_params = new_params, num_repl = 50)

filt_data <- out_data %>%
  filter(Rho == "One Resource") %>%
  group_by(Ctype, Cd, Ptype, rho_mult) %>%
  summarise(ProbFeasible = mean(Feas),
            ProbStable = mean(Stab),
            ProbLLFeas = mean(LLFeas),
            ProbLLSuccess = mean(LLSuccess),
            ProbSubSuccess = mean(SubSuccess)) %>%
  mutate(Label = ifelse(Ctype == "Tradeoff", "(A) Tradeoff",
                       ifelse(Ctype == "Random", "(B) Random", "(C) Banded")))

plHeatMap <- ggplot(filt_data, aes(x = rho_mult, y = Cd, fill = ProbSubSuccess)) +
  geom_tile() + scale_fill_viridis(end = 0.9, na.value="darkred") +
  theme(strip.background = element_rect(fill="lightgray", color = "black"),
        aspect.ratio = 1,
        axis.text = element_text( size = 12.5 ),
        panel.background = element_rect(fill = NA),
        text = element_text(size=30),
        legend.key.size = unit(2, "cm"),
        legend.key.width = unit(1,"cm")) +
  facet_grid(~Label) +
  scale_x_log10(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression("Resource Inflow:"~rho),
       y = expression("Diagonal Consumption:"~C[d]),
       fill = expression("Probability")) +
  ggtitle("Feasible and stable equilibria under structured Liebig's law")
show(plHeatMap)

tiff("../figs/FigHeatMap.tiff", width = 2000, height = 800, res = 125)
plHeatMap
dev.off()
