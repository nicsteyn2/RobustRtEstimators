
rm(list=ls())

library(tidyverse)
library(cowplot)

# Load the data
df = read.csv("paper/outputs/s7_swll_data.csv")
df_ee = read.csv("paper/outputs/s7_swll_epiestim.csv") %>% mutate(Model="EpiEstim")
df_ef = read.csv("paper/outputs/s7_swll_epifilter.csv") %>% mutate(Model="EpiFilter")

df_ee$swll[df_ee$t <= 10] = NaN
df_ef$swll[df_ef$t <= 10] = NaN

custom_theme <- theme_bw() +
  theme(
    axis.title = element_text(size=8),
    axis.text = element_text(size=7.5),
    strip.background = element_blank(),
    strip.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=8),
    legend.title = element_text(size=8)
  )



# Plot random walk simulation results
plt_Rt_rw = ggplot(df %>% filter(simulation=="Random walk"), aes(x=t, y=Rt)) +
  geom_line() +
  custom_theme +
  theme_bw() + 
  ylab("Reproduction number") + xlab("") +
  ggtitle("Random walk Rt simulation")

plt_cases_rw = ggplot(df %>% filter(simulation=="Random walk"), aes(x=t, y=Ct)) +
  geom_col() +
  custom_theme +
  theme_bw() +
  ylab("Reported cases") + xlab("")

plt_ee_rw = ggplot(df_ee %>% filter(simulation=="Random walk"), aes(x=t, y=swll, color=k, group=k)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("EpiEstim log-likelihood") + xlab("") +
  theme(legend.position="none")

plt_ef_rw = ggplot(df_ef %>% filter(simulation=="Random walk"), aes(x=t, y=swll, color=eta, group=eta)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("EpiFilter log-likelihood") + xlab("Time (days)") +
  theme(legend.position="none")

plt_rw = plot_grid(plt_Rt_rw, plt_cases_rw, plt_ee_rw, plt_ef_rw, ncol=1, align="v")



# Plot sinusoidal simulation results
plt_Rt_sin = ggplot(df %>% filter(simulation=="Sinusoidal"), aes(x=t, y=Rt)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("Reproduction number") + xlab("") +
  ggtitle("Sinusoidal Rt simulation")

plt_cases_sin = ggplot(df %>% filter(simulation=="Sinusoidal"), aes(x=t, y=Ct)) +
  geom_col() +
  custom_theme +
  theme_bw() +
  ylab("Reported cases") + xlab("")

plt_ee_sin = ggplot(df_ee %>% filter(simulation=="Sinusoidal"), aes(x=t, y=swll, color=k, group=k)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("EpiEstim log-likelihood") + xlab("") +
  theme(legend.position="none")

plt_ef_sin = ggplot(df_ef %>% filter(simulation=="Sinusoidal"), aes(x=t, y=swll, color=eta, group=eta)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("EpiFilter log-likelihood") + xlab("Time (days)") +
  theme(legend.position="none")

plt_sin = plot_grid(plt_Rt_sin, plt_cases_sin, plt_ee_sin, plt_ef_sin, ncol=1, align="v")



# Plot step-change simulation results
plt_Rt_sc = ggplot(df %>% filter(simulation=="Step-change"), aes(x=t, y=Rt)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("Reproduction number") + xlab("")+
  ggtitle("Step change Rt simulation")

plt_cases_sc = ggplot(df %>% filter(simulation=="Step-change"), aes(x=t, y=Ct)) +
  geom_col() +
  custom_theme +
  theme_bw() +
  ylab("Reported cases") + xlab("")

plt_ee_sc = ggplot(df_ee %>% filter(simulation=="Step-change"), aes(x=t, y=swll, color=k, group=k)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("EpiEstim log-likelihood") + xlab("")

plt_ef_sc = ggplot(df_ef %>% filter(simulation=="Step-change"), aes(x=t, y=swll, color=eta, group=eta)) +
  geom_line() +
  custom_theme +
  theme_bw() +
  ylab("EpiFilter log-likelihood") + xlab("Time (days)")

plt_sc = plot_grid(plt_Rt_sc, plt_cases_sc, plt_ee_sc, plt_ef_sc, ncol=1, align="v")


plt = plot_grid(plt_rw, plt_sin, plt_sc, ncol=3, rel_widths=c(1,1,1.1))
plt
ggsave("paper/figures/s7_stepwiselikelihoods.png", plt, width=25, height=20, dpi=300, units="cm")
ggsave("paper/figures/s7_stepwiselikelihoods.pdf", plt, width=25, height=20, dpi=300, units="cm")





