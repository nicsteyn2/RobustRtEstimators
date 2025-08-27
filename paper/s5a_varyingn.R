
rm(list=ls())

library(tidyverse)
library(cowplot)

custom_theme <- theme_bw() +
  theme(
    axis.title = element_text(size=8),
    axis.text = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=8),
    legend.title = element_text(size=8)
  )

# Load the data
df = read.csv("paper/outputs/s5a_varying_C0_sinusoidal.csv")


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
  summarise(Score = mean(crps))#

df_scores_tmp = df_scores %>%
  ungroup() %>%
  filter(model=="EpiEstim", fit=="Default") %>%
  select(-model, -fit) %>%
  rename(EpiEstimScore=Score)

df_scores = left_join(df_scores, df_scores_tmp, by=c("iter", "simulation", "TrueParam")) %>%
  mutate(RelativeScore = Score/EpiEstimScore)

# 
# # Plot estimates (just to check)
# CASES = 3200
# ITER = 1
# df_tmp = df %>% filter(simulation=="Sinusoidal", TrueParam==CASES, iter==ITER)  %>% mutate(NameForColors=paste0(model, fit))
# 
# colors <- c(EpiEstimMarginalised="#006666", EpiEstimDefault="#20B2AA",
#             EpiFilterMarginalised="#DE7139", EpiFilterDefault="#E5785A")
# 
# 
# plt_data_sinusoidal = ggplot(df_tmp %>% filter(model=="EpiEstim", fit=="Marginalised")) +
#   geom_col(aes(x=t, y=Cases), width = 0.8) +
#   custom_theme +
#   xlab("Time (days)") + ylab("") +
#   ggtitle("Sinusoidal Rt simulation")
# 
# plt_post_sinusoidal = ggplot(df_tmp) +
#   geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
#   geom_line(aes(x=t, y=mean, color=NameForColors), data=.%>%filter(t>10)) +
#   geom_line(aes(x=t, y=TrueRt), color="black", lwd=0.3) +
#   facet_grid(fit~model) +
#   ylab("") + xlab("Time (days)") +
#   scale_fill_manual(values=colors) +
#   scale_color_manual(values=colors) +
#   custom_theme +
#   theme(legend.position="none")
# 
# plt_predpost_sinusoidal = ggplot(df_tmp) +
#   geom_ribbon(aes(x=t, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
#   geom_line(aes(x=t, y=mean_cases, color=NameForColors), data=.%>%filter(t>10)) +
#   geom_point(aes(x=t, y=Cases), color="black", size=0.3) +
#   facet_grid(fit~model) +
#   ylab("") + xlab("Time (days)") +
#   scale_fill_manual(values=colors) +
#   scale_color_manual(values=colors) +
#   custom_theme +
#   theme(legend.position="none")
# 
# plt_check = plot_grid(plt_data_sinusoidal, plt_post_sinusoidal, plt_predpost_sinusoidal, nrow=3)
# plt_check


# Set plotting options
parampointsize = 0.67
paramlinewidth = 0.25

# PLOT SINUSOIDAL RESULTS
df_params_sin = df %>%
  filter(t==100, simulation=="Sinusoidal", fit=="Marginalised") %>%
  mutate(TrueParamJit = exp(log(TrueParam) + (iter-5.5)*0.01))

plt_sin_param_ef = ggplot(df_params_sin %>% filter(model=="EpiFilter")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#DE7139", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#DE7139", linewidth=paramlinewidth) +
  custom_theme +
  xlab("") + ylab("Estimated eta") +
  geom_hline(yintercept=0.1, linetype="dotted") +
  scale_x_log10(breaks=c(25, 50, 100, 200, 400, 800, 1600, 3200), 
                labels=c(25, 50, 100, 200, 400, 800, 1600, 3200)) +
  ylim(c(0,0.3))

plt_sin_param_ee = ggplot(df_params_sin %>% filter(model=="EpiEstim")) +
  geom_point(aes(x=TrueParamJit, y=mean_param), color="#006666", size=parampointsize) +
  geom_linerange(aes(x=TrueParamJit, ymin=lower_param, ymax=upper_param), color="#006666", linewidth=paramlinewidth) +
  custom_theme +
  scale_y_continuous(breaks=seq(1, 30), labels=seq(1,30), minor_breaks=NULL, limits = c(1, 10)) +
  xlab("") + ylab("Estimated k") +
  geom_hline(yintercept=7, linetype="dotted") +
  scale_x_log10(breaks=c(25, 50, 100, 200, 400, 800, 1600, 3200), 
                labels=c(25, 50, 100, 200, 400, 800, 1600, 3200))

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
  ylim(0, 100) +
  scale_x_log10(breaks=c(25, 50, 100, 200, 400, 800, 1600, 3200), 
                labels=c(25, 50, 100, 200, 400, 800, 1600, 3200))

plt_sin_scores = ggplot(df_scores %>% filter(simulation=="Sinusoidal"), aes(x=TrueParam)) +
  geom_point(aes(y=RelativeScore, color=model, shape=fit), alpha=pointalpha) +
  geom_line(aes(y=RelativeScore, color=model, linetype=fit), data= . %>% group_by(TrueParam, model, fit) %>% summarise(RelativeScore=mean(RelativeScore)), lwd=linewidth) +
  scale_color_manual(values=colors) + scale_shape_manual(values=c(16,4)) +
  custom_theme +
  theme(legend.position="bottom") +
  labs(color="Model", shape="Fit", linetype="Fit") +
  scale_y_log10() + ylab("Relative CRPS") + xlab("Number of initial cases (log-scale)") +
  scale_x_log10(breaks=c(25, 50, 100, 200, 400, 800, 1600, 3200), 
                labels=c(25, 50, 100, 200, 400, 800, 1600, 3200))
 
plt_sin_stats = plot_grid(plt_sin_cov, plt_sin_params, plt_sin_scores, ncol=1, rel_heights = c(1.4,1,1.2))
plt_sin_stats
ggsave("paper/figures/s5a_varyingn.png", plt_sin_stats, width=12, height=15, dpi=300, units="cm")
ggsave("paper/figures/s5a_varyingn.pdf", plt_sin_stats, width=12, height=15, dpi=600, units="cm")

