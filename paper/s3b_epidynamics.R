
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
  read.csv("paperv2/outputs/s3b_varying_rw.csv"),
  read.csv("paperv2/outputs/s3b_varying_sinusoidal.csv"),
  read.csv("paperv2/outputs/s3b_varying_stepchanges.csv")
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

df_eepost = rbind(
  read.csv("paperv2/outputs/s3b_varying_rw_eepost.csv"),
  read.csv("paperv2/outputs/s3b_varying_sinusoidal_eepost.csv"),
  read.csv("paperv2/outputs/s3b_varying_stepchanges_eepost.csv")
) %>% mutate(TrueParam = case_when(TrueParam==51 ~ 1,
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

df_efpost = rbind(
  read.csv("paperv2/outputs/s3b_varying_rw_efpost.csv"),
  read.csv("paperv2/outputs/s3b_varying_sinusoidal_efpost.csv"),
  read.csv("paperv2/outputs/s3b_varying_stepchanges_efpost.csv")
) %>% mutate(TrueParam = case_when(TrueParam==51 ~ 1,
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


df %>% group_by(simulation) %>% summarise(maxiter=max(iter))
df %>% filter(iter==10, simulation=="Sinusoidal") %>% summarise(maxparam=max(TrueParam))

# Tidy up
df$fit = ifelse(grepl("Conditional", df$fit), "Default", df$fit)
df = df %>% mutate(NameForColors=paste0(model, fit),
                   fit = factor(fit, levels=c("Marginalised", "Default")))


# Append colors
colors <- c(EpiEstim="#006666", EpiFilter="#DE7139")
pointalpha = 0.25
linewidth = 1


# Calculate statistics
df_coverage = df %>%
  filter(t > 10) %>%
  group_by(iter, simulation, model, fit, TrueParam) %>%
  summarise(Rt_coverage = 100*round(mean((TrueRt <= upper) & (TrueRt >= lower)), digits=3),
            Ct_coverage = 100*round(mean((Cases <= upper_cases) & (Cases >= lower_cases)), digits=3)) %>%
  pivot_longer(cols=c(Rt_coverage, Ct_coverage), names_to="Variable", values_to="Value") %>%
  mutate(Variable=factor(case_when(Variable=="Rt_coverage" ~ "Reproduction number",
                            Variable=="Ct_coverage" ~ "Predictive cases"),
                         levels=c("Reproduction number", "Predictive cases")))

df_coverage_mean = df_coverage %>%
  group_by(simulation, model, fit, Variable, TrueParam) %>%
  summarise(Value = mean(Value))

df_scores = df %>%
  filter(t > 10) %>%
  group_by(iter, simulation, model, fit, TrueParam) %>%
  summarise(Score = mean(crps))

df_scores_tmp = df_scores %>%
  ungroup() %>%
  filter(model=="EpiEstim", fit=="Default") %>%
  select(-model, -fit) %>%
  rename(EpiEstimScore=Score)

df_scores = left_join(df_scores, df_scores_tmp, by=c("iter", "simulation", "TrueParam")) %>%
  mutate(RelativeScore = Score/EpiEstimScore)



# Set plotting options
parampointsize = 0.67
paramlinewidth = 0.25

# PLOT RANDOM WALK RESULTS
df_params_rw = df %>%
  filter(t==100, simulation=="Random walk", fit=="Marginalised") %>%
  mutate(TrueParamJit = TrueParam + (iter-5.5)*0.001)

plt_rw_param_ef = ggplot(df_params_rw %>% filter(model=="EpiFilter")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#DE7139", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#DE7139", linewidth=paramlinewidth) +
  custom_theme +
  xlab("") + ylab("Estimated eta") +
  geom_abline(intercept=0, slope=1) +
  geom_hline(yintercept=0.1, linetype="dotted")

plt_rw_param_ee = ggplot(df_params_rw %>% filter(model=="EpiEstim")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#006666", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#006666", linewidth=paramlinewidth) +
  custom_theme +
  scale_y_continuous(breaks=seq(1, 30, by=2), labels=seq(1,30, by=2), minor_breaks=NULL) +
  xlab("") + ylab("Estimated k") +
  geom_hline(yintercept=7, linetype="dotted")

plt_rw_params = plot_grid(plt_rw_param_ef, plt_rw_param_ee, ncol=2)

plt_rw_cov = ggplot(df_coverage %>% filter(simulation=="Random walk"), aes(x=TrueParam)) +
  geom_point(aes(y=Value, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=Value, color=model, linetype=fit), data=df_coverage_mean %>% filter(simulation=="Random walk"), lwd=linewidth) +
  geom_hline(yintercept=95, linetype="dashed") +
  facet_wrap(~Variable, ncol=2) +
  ylab("Coverage (%)") + xlab("") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) +  scale_shape_manual(values=c(16,4)) +
  ylim(0, 100)  + ggtitle("Random walk Rt simulation")

plt_rw_scores = ggplot(df_scores %>% filter(simulation=="Random walk"), aes(x=TrueParam)) +
  geom_point(aes(y=RelativeScore, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=RelativeScore, color=model, linetype=fit), data= . %>% group_by(TrueParam, model, fit) %>% summarise(RelativeScore=mean(RelativeScore)), lwd=linewidth) +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  custom_theme +
  theme(legend.position="none") +
  ylab("Relative CRPS") + xlab("Standard deviation of random walk")

plt_rw_stats = plot_grid(plt_rw_cov, plt_rw_params, plt_rw_scores, ncol=1, rel_heights = c(1.4,1,1))
plt_rw_stats




# PLOT SINUSOIDAL RESULTS
df_params_sin = df %>%
  filter(t==100, simulation=="Sinusoidal", fit=="Marginalised") %>%
  mutate(TrueParamJit = TrueParam + (iter-5.5)*0.1)

plt_sin_param_ef = ggplot(df_params_sin %>% filter(model=="EpiFilter")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#DE7139", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#DE7139", linewidth=paramlinewidth) +
  custom_theme +
  xlab("") + ylab("Estimated eta") +
  geom_hline(yintercept=0.1, linetype="dotted")

plt_sin_param_ee = ggplot(df_params_sin %>% filter(model=="EpiEstim")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#006666", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#006666", linewidth=paramlinewidth) +
  custom_theme +
  scale_y_continuous(breaks=seq(1, 30), labels=seq(1,30), minor_breaks=NULL, limits = c(1, 10)) +
  xlab("") + ylab("Estimated k") +
  geom_hline(yintercept=7, linetype="dotted")

plt_sin_params = plot_grid(plt_sin_param_ef, plt_sin_param_ee, ncol=2)

plt_sin_cov = ggplot(df_coverage %>% filter(simulation=="Sinusoidal"), aes(x=TrueParam)) +
  geom_point(aes(y=Value, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=Value, color=model, linetype=fit), data=df_coverage_mean %>% filter(simulation=="Sinusoidal"), lwd=linewidth) +
  geom_hline(yintercept=95, linetype="dashed") +
  facet_wrap(~Variable, ncol=2) +
  ylab("Coverage (%)") + xlab("") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  ggtitle("Sinusoidal Rt simulation") +
  ylim(0, 100)

plt_sin_scores = ggplot(df_scores %>% filter(simulation=="Sinusoidal"), aes(x=TrueParam)) +
  geom_point(aes(y=RelativeScore, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=RelativeScore, color=model, linetype=fit), data= . %>% group_by(TrueParam, model, fit) %>% summarise(RelativeScore=mean(RelativeScore)), lwd=linewidth) +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  custom_theme +
  theme(legend.position="none") +
  scale_y_log10() + ylab("Relative CRPS") + xlab("Period of Rt simulation")

plt_sin_stats = plot_grid(plt_sin_cov, plt_sin_params, plt_sin_scores, ncol=1, rel_heights = c(1.4,1,1))
plt_sin_stats



# PLOT STEP-CHANGE RESULTS
df_params_sc = df %>%
  filter(t==100, simulation=="Step change", fit=="Marginalised") %>%
  mutate(TrueParamJit = TrueParam + (iter-5.5)*0.02)

plt_sc_param_ef = ggplot(df_params_sc %>% filter(model=="EpiFilter")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#DE7139", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#DE7139", linewidth=paramlinewidth) +
  custom_theme +
  xlab("") + ylab("Estimated eta") +
  geom_hline(yintercept=0.1, linetype="dotted")

plt_sc_param_ee = ggplot(df_params_sc %>% filter(model=="EpiEstim")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#006666", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#006666", linewidth=paramlinewidth) +
  custom_theme +
  scale_y_continuous(breaks=seq(1, 30), labels=seq(1,30), minor_breaks=NULL, limits = c(1, 10)) +
  xlab("") + ylab("Estimated k") +
  geom_hline(yintercept=7, linetype="dotted")

plt_sc_params = plot_grid(plt_sc_param_ef, plt_sc_param_ee, ncol=2)

plt_sc_cov = ggplot(df_coverage %>% filter(simulation=="Step change"), aes(x=TrueParam)) +
  geom_point(aes(y=Value, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=Value, color=model, linetype=fit), data=df_coverage_mean %>% filter(simulation=="Step change"), lwd=linewidth) +
  geom_hline(yintercept=95, linetype="dashed") +
  facet_wrap(~Variable, ncol=2) +
  ylab("Coverage (%)") + xlab("") +
  custom_theme +
  theme(legend.position="none") +
  scale_color_manual(values=colors) +  scale_shape_manual(values=c(16,4)) +
  ggtitle("Step-change Rt simulation") +
  ylim(0, 100)

plt_sc_scores = ggplot(df_scores %>% filter(simulation=="Step change"), aes(x=TrueParam)) +
  geom_point(aes(y=RelativeScore, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=RelativeScore, color=model, linetype=fit), data= . %>% group_by(TrueParam, model, fit) %>% summarise(RelativeScore=mean(RelativeScore)), lwd=linewidth) +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  custom_theme +
  theme(legend.direction="horizontal") + labs(color="Model", shape="Fit", linetype="Fit") +
  ylab("Relative CRPS") + xlab("Number of step-changes in Rt")

legend = get_legend(plt_sc_scores)


plt_sc_stats = plot_grid(plt_sc_cov, plt_sc_params, plt_sc_scores + theme(legend.position="none"), ncol=1, rel_heights = c(1.4,1,1))
plt_sc_stats


plt = plot_grid(plt_rw_stats, plt_sin_stats, plt_sc_stats, NULL, legend, NULL, rel_heights=c(1,0.2), ncol=3)
plt
ggsave("paperv2/figures/s3b_epidynamics.png", plt, width=25, height=18, dpi=300, units="cm")

