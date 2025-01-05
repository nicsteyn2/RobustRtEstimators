
rm(list=ls())

library(tidyverse)
library(cowplot)

custom_theme <- theme_bw() +
  theme(
    axis.title = element_text(size=8),
    axis.text = element_text(size=8),
    strip.background = element_blank(),
    strip.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5)
  )

# Load the data
df = rbind(
  read.csv("paper/outputs/s5b_varying_rw.csv"),
  read.csv("paper/outputs/s5b_varying_sinusoidal.csv"),
  read.csv("paper/outputs/s5b_varying_stepchanges.csv")
)  %>% mutate(TrueParam = case_when(TrueParam==51 ~ 1,
                                    TrueParam==34 ~ 2,
                                    TrueParam==26 ~ 3,
                                    TrueParam==21 ~ 4,
                                    TrueParam==18 ~ 5,
                                    simulation=="Step change" & TrueParam==15 ~ 6,
                                    TrueParam==14 ~ 7,
                                    TrueParam==12 ~ 8,
                                    TrueParam==11 ~ 9,
                                    TrueParam==10 ~ 10,
                                    TRUE ~ TrueParam))

# Tidy up
df$fit = ifelse(grepl("Conditional", df$fit), "Default", df$fit)
df = df %>% mutate(NameForColors=paste0(model, fit),
                   fit = factor(fit, levels=c("Marginalised", "Default")))

# Append colors
colors <- c(EpiEstim="#006666", EpiFilter="#DE7139")
pointalpha = 0.25

# Set other plotting params
ylims_correct = c(0, 100)
ylims_incorrect = c(0, 45)

# Calculate practical statistics
df_sign = df %>%
  mutate(is_confident = (lower > 1 | upper < 1),
         correctly_confident = is_confident & ((TrueRt < 1 & upper < 1) | (TrueRt > 1 & upper > 1)),
         incorrectly_confident = is_confident & ((TrueRt <= 1 & lower > 1) | (TrueRt >= 1 & upper < 1))) %>%
  filter(t > 10) %>%
  group_by(iter, simulation, model, fit, TrueParam) %>%
  summarise(propTimeConfident = 100*mean(is_confident),
            propTimeCorrectlyConfident = 100*mean(correctly_confident),
            propTimeIncorrectlyConfident = 100*mean(incorrectly_confident))


# Plot random walk results
plt_rw_cc = ggplot(df_sign %>% filter(simulation=="Random walk")) +
  geom_point(aes(x=TrueParam, y=propTimeCorrectlyConfident, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(x=TrueParam, y=y, color=model, linetype=fit), data=.%>%group_by(TrueParam, model, fit) %>% summarise(y=mean(propTimeCorrectlyConfident))) +
  ylab("Time correctly confident (%)") + xlab("") + ggtitle("Random walk Rt simulation") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ylim(ylims_correct)

plt_rw_ic = ggplot(df_sign %>% filter(simulation=="Random walk")) +
  geom_point(aes(x=TrueParam, y=propTimeIncorrectlyConfident, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(x=TrueParam, y=y, color=model, linetype=fit), data=.%>%group_by(TrueParam, model, fit) %>% summarise(y=mean(propTimeIncorrectlyConfident))) +
  ylab("Time incorrectly confident (%)") + xlab("Standard deviation of random walk") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ylim(ylims_incorrect)

plt_rw = plot_grid(plt_rw_cc, plt_rw_ic, ncol=1)


# Plot sinusoidal results
plt_sin_cc = ggplot(df_sign %>% filter(simulation=="Sinusoidal")) +
  geom_point(aes(x=TrueParam, y=propTimeCorrectlyConfident, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(x=TrueParam, y=y, color=model, linetype=fit), data=.%>%group_by(TrueParam, model, fit) %>% summarise(y=mean(propTimeCorrectlyConfident))) +
  ylab("Time correctly confident (%)") + xlab("") + ggtitle("Sinusoidal Rt simulation") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ylim(ylims_correct)

plt_sin_ic = ggplot(df_sign %>% filter(simulation=="Sinusoidal")) +
  geom_point(aes(x=TrueParam, y=propTimeIncorrectlyConfident, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(x=TrueParam, y=y, color=model, linetype=fit), data=.%>%group_by(TrueParam, model, fit) %>% summarise(y=mean(propTimeIncorrectlyConfident))) +
  ylab("Time incorrectly confident (%)") + xlab("Period of Rt simulation") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ylim(ylims_incorrect)

plt_sin = plot_grid(plt_sin_cc, plt_sin_ic, ncol=1)


# Plot step changes results
plt_sc_cc = ggplot(df_sign %>% filter(simulation=="Step change")) +
  geom_point(aes(x=TrueParam, y=propTimeCorrectlyConfident, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(x=TrueParam, y=y, color=model, linetype=fit), data=.%>%group_by(TrueParam, model, fit) %>% summarise(y=mean(propTimeCorrectlyConfident))) +
  ylab("Time correctly confident (%)") + xlab("") + ggtitle("Step-change Rt simulation") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ylim(ylims_correct)

plt_sc_ic = ggplot(df_sign %>% filter(simulation=="Step change")) +
  geom_point(aes(x=TrueParam, y=propTimeIncorrectlyConfident, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(x=TrueParam, y=y, color=model, linetype=fit), data=.%>%group_by(TrueParam, model, fit) %>% summarise(y=mean(propTimeIncorrectlyConfident))) +
  ylab("Time incorrectly confident (%)") + xlab("Number of step-changes in Rt") +
  custom_theme +
  theme(legend.direction="horizontal") + labs(color="Model", shape="Fit", linetype="Fit") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ylim(ylims_incorrect)

legend = get_legend(plt_sc_ic)

plt_sc = plot_grid(plt_sc_cc, plt_sc_ic  + theme(legend.position="none"), ncol=1)


plt = plot_grid(plot_grid(plt_rw, plt_sin, plt_sc, ncol=3), legend, ncol=1, rel_heights=c(1,0.2))
plt
ggsave("paper/figures/s5c_confidence.png", plt, width=25, height=12, dpi=300, units="cm")
