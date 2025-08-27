

rm(list=ls())

library(tidyverse)
library(cowplot)

custom_theme <- theme_bw() +
  theme(
    axis.title.y = element_text(size=9),
    axis.title.x = element_text(size=8),
    axis.text = element_text(size=8),
    strip.background = element_blank(),
    strip.text = element_text(size=9),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Load the data
df = rbind(read.csv("paper/outputs/s4a_epiestim_siuncertainty.csv"), read.csv("paper/outputs/s4a_epifilter_siuncertainty.csv"))
df_params = rbind(
  read.csv("paper/outputs/s4a_epifilter_siuncertainty_params.csv") %>% rename(x=eta, y=peta),
  read.csv("paper/outputs/s4a_epiestim_siuncertainty_params.csv") %>% rename(x=k, y=pk)
)

# Calculate coverage statistics
df_coverage = df %>%
  filter(t>10) %>%
  group_by(simulation, model, fit) %>%
  summarise(Rt_coverage = 100*round(mean((TrueRt <= upper) & (TrueRt >= lower)), digits=3),
            Ct_coverage = 100*round(mean((Cases <= upper_cases) & (Cases >= lower_cases)), digits=3),
            score = signif(mean(crps), 3))


# Make individual plots
plt_param_ef = ggplot(df_params %>% filter(model=="EpiFilter")) +
  geom_line(aes(x=x, y=y, color=fit)) +
  xlim(c(0, 0.3)) +
  theme_bw() +
  labs(color="") + theme(legend.direction="horizontal") +
  ggtitle("EpiFilter") +
  xlab("eta") + ylab("")

plt_Rt_ef = ggplot(df %>% filter(model=="EpiFilter")) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=fit), alpha=0.3, data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=mean, color=fit), data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=TrueRt), color="black", lwd=0.3) +
  geom_text(data = df_coverage %>% filter(fit=="No Uncertainty", model=="EpiFilter"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%"), color=fit),
            hjust = 0, vjust = 1.5, size = 3) +
  geom_text(data = df_coverage %>% filter(fit=="SI Uncertainty", model=="EpiFilter"),
            aes(x = 100, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%"), color=fit),
            hjust = 1, vjust = 1.5, size = 3) +
  custom_theme +
  xlab("Time (days)") + ylab("") +
  ylim(0.35, 1.7)

plt_Ct_ef = ggplot(df %>% filter(model=="EpiFilter")) +
  geom_ribbon(aes(x=t, ymin=lower_cases, ymax=upper_cases, fill=fit), alpha=0.3, data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=mean_cases, color=fit), data=.%>%filter(t>10)) +
  geom_point(aes(x=t, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(fit=="No Uncertainty", model=="EpiFilter"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score), color=fit),
            hjust = 0, vjust = 1.2, size = 3) +
  geom_text(data = df_coverage %>% filter(fit=="SI Uncertainty", model=="EpiFilter"),
            aes(x = 100, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score), color=fit),
            hjust = 1, vjust = 1.2, size = 3) +
  custom_theme +
  xlab("Time (days)") + ylab("") +
  ylim(0,400)
    
plt_param_ee = ggplot(df_params %>% filter(model=="EpiEstim")) +
  geom_step(aes(x=x, y=y, color=fit)) +
  custom_theme +
  ggtitle("EpiEstim") +
  xlab("k") + ylab("Param. posterior")

plt_Rt_ee = ggplot(df %>% filter(model=="EpiEstim")) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=fit), alpha=0.3, data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=mean, color=fit), data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=TrueRt), color="black", lwd=0.3) +
  geom_text(data = df_coverage %>% filter(fit=="No Uncertainty", model=="EpiEstim"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%"), color=fit),
            hjust = 0, vjust = 1.5, size = 3) +
  geom_text(data = df_coverage %>% filter(fit=="SI Uncertainty", model=="EpiEstim"),
            aes(x = 100, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%"), color=fit),
            hjust = 1, vjust = 1.5, size = 3) +
  custom_theme +
  xlab("Time (days)") + ylab("Reproduction number") +
  ylim(0.35, 1.7)

plt_Ct_ee = ggplot(df %>% filter(model=="EpiEstim")) +
  geom_ribbon(aes(x=t, ymin=lower_cases, ymax=upper_cases, fill=fit), alpha=0.3, data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=mean_cases, color=fit), data=.%>%filter(t>10)) +
  geom_point(aes(x=t, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(fit=="No Uncertainty", model=="EpiEstim"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score), color=fit),
            hjust = 0, vjust = 1.2, size = 3) +
  geom_text(data = df_coverage %>% filter(fit=="SI Uncertainty", model=="EpiEstim"),
            aes(x = 100, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score), color=fit),
            hjust = 1, vjust = 1.2, size = 3) +
  custom_theme +
  xlab("Time (days)") + ylab("Reported cases") +
  ylim(0,400)


# Put it all together
legend = get_legend(plt_param_ef)
plt = plot_grid(plt_param_ee, plt_param_ef + custom_theme, plt_Rt_ee, plt_Rt_ef, plt_Ct_ee, plt_Ct_ef, ncol=2)
plt_with_legend = plot_grid(plt, legend, ncol=1, rel_heights=c(1, 0.05))
plt_with_legend
ggsave(paste0("paper/figures/s4a_gentimeuncertainty.png"), plt_with_legend, dpi=300, width=20, height=14, units="cm")
ggsave(paste0("paper/figures/s4a_gentimeuncertainty.pdf"), plt_with_legend, dpi=600, width=20, height=14, units="cm")



