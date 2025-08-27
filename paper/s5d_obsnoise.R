
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
    panel.grid = element_blank()
  )

# Load the data
df = read.csv("paper/outputs/s5d_obsnoise.csv")
df_eepost = read.csv("paper/outputs/s5d_obsnoise_epiestimposteriors.csv")
df_efpost = read.csv("paper/outputs/s5d_obsnoise_epifilterposteriors.csv")

# Clean the data
df$fit = ifelse(grepl("Conditional", df$fit), "Conditional", df$fit)
df = df %>% mutate(NameForColors=paste0(model, fit))

# Append colors
colors <- c(EpiEstimMarginalised="#006666", EpiEstimConditional="#20B2AA",
            EpiFilterMarginalised="#DE7139", EpiFilterConditional="#E5785A")

# Calculate coverage statistics
df_coverage = df %>%
  filter(t>10) %>%
  group_by(simulation, model, fit) %>%
  summarise(Rt_coverage = 100*round(mean((TrueRt <= upper) & (TrueRt >= lower)), digits=3),
            Ct_coverage = 100*round(mean((Cases <= upper_cases) & (Cases >= lower_cases)), digits=3),
            score = signif(mean(crps), 3))



# Plot the raw data
plt_data_rw = ggplot(df %>% filter(simulation=="Random walk", model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=t, y=Cases), width = 0.8) +
  scale_x_continuous(limits = c(0, 100)) +
  custom_theme +
  xlab("") + ylab("Reported cases") +
  ggtitle("Random walk Rt simulation")

plt_data_sinusoidal = ggplot(df %>% filter(simulation=="Sinusoidal", model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=t, y=Cases), width = 0.8) +
  custom_theme +
  xlab("Time (days)") + ylab("") +
  ggtitle("Sinusoidal Rt simulation")

plt_data_steps = ggplot(df %>% filter(simulation=="Step change", model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=t, y=Cases), width = 0.8) +
  custom_theme +
  xlab("") + ylab("") +
  ggtitle("Step change Rt simulation")

plt_cases = plot_grid(plt_data_rw, plt_data_sinusoidal, plt_data_steps, nrow=1, align="h")



# Plot posteriors
plt_post_rw = ggplot(df %>% filter(simulation=="Random walk")) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=t, y=mean, color=NameForColors), data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=TrueRt), color="black", lwd=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == "Random walk"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%")),
            hjust = 0, vjust = 1.5, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0.2, 2.0)) +
  ylab("Reproduction number") + xlab("") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_post_sinusoidal = ggplot(df %>% filter(simulation=="Sinusoidal")) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=t, y=mean, color=NameForColors), data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=TrueRt), color="black", lwd=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == "Sinusoidal"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%")),
            hjust = 0, vjust = 1.5, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0.2, 2.0)) +
  ylab("") + xlab("Time (days)") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  custom_theme +
  theme(legend.position="none")

plt_post_steps = ggplot(df %>% filter(simulation=="Step change")) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=t, y=mean, color=NameForColors), data=.%>%filter(t>10)) +
  geom_line(aes(x=t, y=TrueRt), color="black", lwd=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == "Step change"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Rt_coverage, "%")),
            hjust = 0, vjust = 1.5, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0.2, 2.0)) +
  ylab("") + xlab("") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  custom_theme +
  theme(legend.position="none")




# Plot predictive posteriors
plt_predpost_rw = ggplot(df %>% filter(simulation=="Random walk")) +
  geom_ribbon(aes(x=t, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=t, y=mean_cases, color=NameForColors), data=.%>%filter(t>10)) +
  geom_point(aes(x=t, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == "Random walk"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0,470)) +
  ylab("Reported cases") + xlab("") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_predpost_sinusoidal = ggplot(df %>% filter(simulation=="Sinusoidal")) +
  geom_ribbon(aes(x=t, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=t, y=mean_cases, color=NameForColors), data=.%>%filter(t>10)) +
  geom_point(aes(x=t, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == "Sinusoidal"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0,350)) +
  ylab("") + xlab("Time (days)") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  custom_theme +
  theme(legend.position="none")

plt_predpost_steps = ggplot(df %>% filter(simulation=="Step change")) +
  geom_ribbon(aes(x=t, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=t, y=mean_cases, color=NameForColors), data=.%>%filter(t>10)) +
  geom_point(aes(x=t, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == "Step change"),
            aes(x = 0, y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  ylab("") + xlab("") +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  custom_theme +
  theme(legend.position="none")






# Plot parameter estimates
xlims = c(0, 0.3)

plt_param_rw_ef = ggplot(df_efpost %>% filter(simulation=="Random walk")) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_rw_ee = ggplot(df_eepost %>% filter(simulation=="Random walk")) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("Posterior density") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(1,10)) +
  geom_vline(xintercept=7, color="black", linetype="dashed")

plt_param_rw = plot_grid(plt_param_rw_ee, plt_param_rw_ef, ncol=2)


plt_param_sin_ef = ggplot(df_efpost %>% filter(simulation=="Sinusoidal")) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dashed") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_sin_ee = ggplot(df_eepost %>% filter(simulation=="Sinusoidal")) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(1,10)) +
  geom_vline(xintercept=7, color="black", linetype="dashed")

plt_param_sinusoidal = plot_grid(plt_param_sin_ee, plt_param_sin_ef, ncol=2)


plt_param_sin_ef = ggplot(df_efpost %>% filter(simulation=="Step change")) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dashed") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_sin_ee = ggplot(df_eepost %>% filter(simulation=="Step change")) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(0,10)) +
  geom_vline(xintercept=7, color="black", linetype="dashed")

plt_param_steps = plot_grid(plt_param_sin_ee, plt_param_sin_ef, ncol=2)


plt = plot_grid(plt_data_rw, plt_data_sinusoidal, plt_data_steps,
                plt_param_rw, plt_param_sinusoidal, plt_param_steps,
                plt_post_rw, plt_post_sinusoidal, plt_post_steps,
                plt_predpost_rw, plt_predpost_sinusoidal, plt_predpost_steps,
                ncol=3, rel_heights = c(0.35,0.25,0.6,0.6)) 
plt
ggsave(paste0("paper/figures/s5d_obsnoise.png"), plt, dpi=300, width=25, height=21, units="cm")
ggsave(paste0("paper/figures/s5d_obsnoise.pdf"), plt, dpi=600, width=25, height=21, units="cm")
 

