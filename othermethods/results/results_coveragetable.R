
rm(list=ls())

library(tidyverse)
library(cowplot)

EpiFilterSmoothed = read.csv("othermethods/EpiFilterSmoothing/results.csv") %>% filter(posterior=="Smoothing") %>% dplyr::select(-posterior)
EpiNow2 = read.csv("othermethods/EpiNow2/results.csv")
EpiNow2InvGamma = read.csv("othermethods/EpiNow2/results_invgamma.csv") %>% dplyr::select(-date)
EpiLPS = read.csv("othermethods/EpiLPSMAP/results.csv")
EpiLPSMALA = read.csv("othermethods/EpiLPSMALA/results.csv")
rtestim = read.csv("othermethods/rtestim/results.csv")

all_methods = rbind(
  EpiFilterSmoothed %>% mutate(model = "EpiFilter (smoothing)"),
  EpiNow2 %>% mutate(model = "EpiNow2 (default)"),
  EpiNow2InvGamma %>% mutate(model = "EpiNow2 (inv-gamma)"),
  EpiLPS %>% mutate(model = "EpiLPS (MAP)"),
  EpiLPSMALA %>% mutate(model = "EpiLPS (MALA)"),
  rtestim %>% mutate(model = "rtestim")
)

all_methods_coverage = all_methods %>%
  filter(t>=10) %>%
  group_by(simulation, variable, model) %>%
  summarise(value = mean(TrueValue >= Lower & TrueValue <= Upper)) %>%
  ungroup()

all_methods_crps = all_methods %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, model) %>%
  summarise(value = mean(CRPS)) %>%
  ungroup()

all_methods_stats = rbind(
  all_methods_coverage %>% mutate(metric = case_when(variable=="I" ~ "Ct coverage", variable=="R" ~ "Rt coverage")) %>% dplyr::select(simulation, model, metric, value),
  all_methods_crps %>% mutate(metric = "CRPS")
) %>% mutate(model = factor(model, levels=c("EpiFilter (smoothing)", "EpiNow2 (default)", "EpiNow2 (inv-gamma)", "EpiLPS (MAP)", "EpiLPS (MALA)", "rtestim")))

# Change line 2 to select the section of table you wish to extract
test = all_methods_stats %>%
  filter(simulation=="Step change") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(model) %>%
  mutate(
    `Rt coverage` = sprintf("%.1f\\%%", `Rt coverage` * 100),
    `Ct coverage` = sprintf("%.1f\\%%", `Ct coverage` * 100),
    `CRPS` = sprintf("%.2f", `CRPS`)
  ) %>%
  mutate(row = sprintf("%s & %s & %s & %s \\\\",
                       model,
                       `Rt coverage`,
                       `Ct coverage`,
                       `CRPS`)) %>%
  pull(row) %>%
  cat(sep = "\n")
