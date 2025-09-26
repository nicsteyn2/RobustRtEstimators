
# Produce plots of the example results in individual model folders.
# This script assumes there exists an appropriately formatted results.csv file within each subdirectory.

rm(list=ls())

library(tidyverse)
library(cowplot)

# General options
inplot_textsize = 3.5
ct_pointsize = 0.5

# EpiLPS (MAP)

EpiLPS = read.csv("othermethods/EpiLPSMAP/results.csv")

EpiLPS_coverage = EpiLPS %>%
  filter(t>=10) %>%
  group_by(simulation, variable) %>%
  summarise(Coverage = mean(Lower<=TrueValue & TrueValue<=Upper))%>%
  mutate(coverage_text = paste0("Coverage = ", round(Coverage*100, 1), "%"))

EpiLPS_crps = EpiLPS %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, variable) %>%
  summarise(CRPS = mean(CRPS))%>%
  mutate(crps_text = paste0("CRPS = ", round(CRPS, 2)))

EpiLPS_text = left_join(EpiLPS_coverage, EpiLPS_crps, by=c("simulation", "variable")) %>%
  mutate(text = case_when(
    variable=="I" ~ paste0(coverage_text, "\n", crps_text),
    variable=="R" ~ paste0(coverage_text)
  ))

EpiLPS_pltR = ggplot(EpiLPS %>% filter(variable=="R"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_line(aes(y=TrueValue), color="red") +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Reproduction number") + ggtitle("EpiLPS (MAP)") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiLPS_text %>% filter(variable=="R"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.5, size = inplot_textsize, color = "black")

EpiLPS_pltI = ggplot(EpiLPS %>% filter(variable=="I"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_point(aes(y=TrueValue), color="red", size=ct_pointsize) +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Infection incidence") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiLPS_text %>% filter(variable=="I"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.2, size = inplot_textsize, color = "black") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))

EpiLPS_plt = plot_grid(EpiLPS_pltR, EpiLPS_pltI, ncol=2, align="h")
EpiLPS_plt



# EpiLPS (MALA)

EpiLPSMALA = read.csv("othermethods/EpiLPSMALA/results.csv")

EpiLPSMALA_coverage = EpiLPSMALA %>%
  filter(t>=10) %>%
  group_by(simulation, variable) %>%
  summarise(Coverage = mean(Lower<=TrueValue & TrueValue<=Upper))%>%
  mutate(coverage_text = paste0("Coverage = ", round(Coverage*100, 1), "%"))

EpiLPSMALA_crps = EpiLPSMALA %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, variable) %>%
  summarise(CRPS = mean(CRPS))%>%
  mutate(crps_text = paste0("CRPS = ", round(CRPS, 2)))

EpiLPSMALA_text = left_join(EpiLPSMALA_coverage, EpiLPSMALA_crps, by=c("simulation", "variable")) %>%
  mutate(text = case_when(
    variable=="I" ~ paste0(coverage_text, "\n", crps_text),
    variable=="R" ~ paste0(coverage_text)
  ))

EpiLPSMALA_pltR = ggplot(EpiLPSMALA %>% filter(variable=="R"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_line(aes(y=TrueValue), color="red") +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Reproduction number") + ggtitle("EpiLPS (MALA)") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiLPSMALA_text %>% filter(variable=="R"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.5, size = inplot_textsize, color = "black")

EpiLPSMALA_pltI = ggplot(EpiLPSMALA %>% filter(variable=="I"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_point(aes(y=TrueValue), color="red", size=ct_pointsize) +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Infection incidence") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiLPSMALA_text %>% filter(variable=="I"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.2, size = inplot_textsize, color = "black") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))

EpiLPSMALA_plt = plot_grid(EpiLPSMALA_pltR, EpiLPSMALA_pltI, ncol=2, align="h")
EpiLPSMALA_plt



# EpiFilter (smoothed version)

EpiFilterSmoothed = read.csv("othermethods/EpiFilterSmoothing/results.csv") %>% filter(posterior=="Smoothing")

EpiFilterSmoothed_coverage = EpiFilterSmoothed %>%
  filter(t>=10) %>%
  group_by(simulation, variable) %>%
  summarise(Coverage = mean(Lower<=TrueValue & TrueValue<=Upper))%>%
  mutate(coverage_text = paste0("Coverage = ", round(Coverage*100, 1), "%"))

EpiFilterSmoothed_crps = EpiFilterSmoothed %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, variable) %>%
  summarise(CRPS = mean(CRPS))%>%
  mutate(crps_text = paste0("CRPS = ", round(CRPS, 2))) 

EpiFilterSmoothed_text = left_join(EpiFilterSmoothed_coverage, EpiFilterSmoothed_crps, by=c("simulation", "variable")) %>%
  mutate(text = case_when(
    variable=="I" ~ paste0(coverage_text, "\n", crps_text),
    variable=="R" ~ paste0(coverage_text)
  ))

EpiFilterSmoothed_pltR = ggplot(EpiFilterSmoothed %>% filter(variable=="R"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_line(aes(y=TrueValue), color="red") +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Reproduction number") + ggtitle("EpiFilter (smoothed)") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiFilterSmoothed_text %>% filter(variable=="R"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.5, size = inplot_textsize, color = "black")

EpiFilterSmoothed_pltI = ggplot(EpiFilterSmoothed %>% filter(variable=="I"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_point(aes(y=TrueValue), color="red", size=ct_pointsize) +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Infection incidence") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiFilterSmoothed_text %>% filter(variable=="I"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.2, size = inplot_textsize, color = "black") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))

EpiFilterSmoothed_plt = plot_grid(EpiFilterSmoothed_pltR, EpiFilterSmoothed_pltI, ncol=2, align="h")
EpiFilterSmoothed_plt



# EpiNow2 (default)

EpiNow2 = read.csv("othermethods/EpiNow2/results.csv")

EpiNow2_coverage = EpiNow2 %>%
  filter(t>=10) %>%
  group_by(simulation, variable) %>%
  summarise(Coverage = mean(Lower<=TrueValue & TrueValue<=Upper))%>%
  mutate(coverage_text = paste0("Coverage = ", round(Coverage*100, 1), "%"))

EpiNow2_crps = EpiNow2 %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, variable) %>%
  summarise(CRPS = mean(CRPS))%>%
  mutate(crps_text = paste0("CRPS = ", round(CRPS, 2))) 

EpiNow2_text = left_join(EpiNow2_coverage, EpiNow2_crps, by=c("simulation", "variable")) %>%
  mutate(text = case_when(
    variable=="I" ~ paste0(coverage_text, "\n", crps_text),
    variable=="R" ~ paste0(coverage_text)
  ))

EpiNow2_pltR = ggplot(EpiNow2 %>% filter(variable=="R"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_line(aes(y=TrueValue), color="red") +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Reproduction number") + ggtitle("EpiNow2 (default)") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiNow2_text %>% filter(variable=="R"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.5, size = inplot_textsize, color = "black")

EpiNow2_pltI = ggplot(EpiNow2 %>% filter(variable=="I"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_point(aes(y=TrueValue), color="red", size=ct_pointsize) +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Infection incidence") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiNow2_text %>% filter(variable=="I"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.2, size = inplot_textsize, color = "black") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))

EpiNow2_plt = plot_grid(EpiNow2_pltR, EpiNow2_pltI, ncol=2, align="h")
EpiNow2_plt


# EpiNow2 (inverse gamma prior distribution for GP lengthcsale)

EpiNow2 = read.csv("othermethods/EpiNow2/results_invgamma.csv")

EpiNow2_coverage = EpiNow2 %>%
  filter(t>=10) %>%
  group_by(simulation, variable) %>%
  summarise(Coverage = mean(Lower<=TrueValue & TrueValue<=Upper))%>%
  mutate(coverage_text = paste0("Coverage = ", round(Coverage*100, 1), "%"))

EpiNow2_crps = EpiNow2 %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, variable) %>%
  summarise(CRPS = mean(CRPS))%>%
  mutate(crps_text = paste0("CRPS = ", round(CRPS, 2))) 

EpiNow2_text = left_join(EpiNow2_coverage, EpiNow2_crps, by=c("simulation", "variable")) %>%
  mutate(text = case_when(
    variable=="I" ~ paste0(coverage_text, "\n", crps_text),
    variable=="R" ~ paste0(coverage_text)
  ))

EpiNow2_pltR = ggplot(EpiNow2 %>% filter(variable=="R"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_line(aes(y=TrueValue), color="red") +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Reproduction number") + ggtitle("EpiNow2 (inverse gamma)") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiNow2_text %>% filter(variable=="R"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.5, size = inplot_textsize, color = "black")

EpiNow2_pltI = ggplot(EpiNow2 %>% filter(variable=="I"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_point(aes(y=TrueValue), color="red", size=ct_pointsize) +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Infection incidence") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = EpiNow2_text %>% filter(variable=="I"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.2, size = inplot_textsize, color = "black") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))

EpiNow2_plt_invgam = plot_grid(EpiNow2_pltR, EpiNow2_pltI, ncol=2, align="h")
EpiNow2_plt_invgam


# EpiNow2 (inverse gamma prior distribution for GP lengthcsale)

rtestim = read.csv("othermethods/rtestim/results.csv")

rtestim_coverage = rtestim %>%
  filter(t>=10) %>%
  group_by(simulation, variable) %>%
  summarise(Coverage = mean(Lower<=TrueValue & TrueValue<=Upper))%>%
  mutate(coverage_text = paste0("Coverage = ", round(Coverage*100, 1), "%"))

rtestim_crps = rtestim %>%
  filter(t>=10, variable=="I") %>%
  group_by(simulation, variable) %>%
  summarise(CRPS = mean(CRPS))%>%
  mutate(crps_text = paste0("CRPS = ", round(CRPS, 2))) 

rtestim_text = left_join(rtestim_coverage, rtestim_crps, by=c("simulation", "variable")) %>%
  mutate(text = case_when(
    variable=="I" ~ paste0(coverage_text, "\n", crps_text),
    variable=="R" ~ paste0(coverage_text)
  ))

rtestim_pltR = ggplot(rtestim %>% filter(variable=="R"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_line(aes(y=TrueValue), color="red") +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Reproduction number") + ggtitle("rtestim") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = rtestim_text %>% filter(variable=="R"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.5, size = inplot_textsize, color = "black")

rtestim_pltI = ggplot(rtestim %>% filter(variable=="I"), aes(x=t)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), data=.%>%filter(t>=10), alpha=0.3) +
  geom_line(aes(y=Mean), data=.%>%filter(t>=10)) +
  geom_point(aes(y=TrueValue), color="red", size=ct_pointsize) +
  facet_wrap(~simulation, ncol=1, scales="free_y") +
  xlab("Time (days)") + ylab("Infection incidence") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  geom_text(data = rtestim_text %>% filter(variable=="I"),
            aes(x = 100, y = Inf, label = text),
            hjust = 1, vjust = 1.2, size = inplot_textsize, color = "black") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))

rtestim_plt = plot_grid(rtestim_pltR, rtestim_pltI, ncol=2, align="h")
rtestim_plt



plt = plot_grid(EpiFilterSmoothed_plt, rtestim_plt, EpiLPS_plt, EpiLPSMALA_plt, EpiNow2_plt, EpiNow2_plt_invgam, ncol=2, align="v")
plt
# ggsave("othermethods/results/othermethods.png", plt, dpi=300, units="cm", width=30, height=30)
ggsave("othermethods/results/othermethods.pdf", plt, dpi=300, units="cm", width=30, height=30)
