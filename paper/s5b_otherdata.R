
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
    plot.title = element_text(hjust = 0.5)#,
    # panel.grid = element_blank()
  )

# Load the data
df = read.csv("paper/outputs/s5otherdatasets.csv")
df_eepost = read.csv("paper/outputs/s5otherdatasets_epiestimposteriors.csv")
df_efpost = read.csv("paper/outputs/s5otherdatasets_epifilterposteriors.csv")

# Set current simulation to plot
CURRENT_SIM = "sars"
labels = c("flufilt"="1918 Influenza (moving average)", "flu"="1918 Influenza", "sars"="2003 SARS", "sarsfilt"="2003 SARS (moving average)")

# Clean the data
df$fit = ifelse(grepl("Conditional", df$fit), "Default", df$fit)
df = df %>% mutate(NameForColors=paste0(model, fit),
                   date = as.Date("2020-02-28") - 1 + t)

# Append colors
colors <- c(EpiEstimMarginalised="#006666", EpiEstimDefault="#20B2AA",
            EpiFilterMarginalised="#DE7139", EpiFilterDefault="#E5785A")

# Calculate coverage statistics
df_coverage = df %>%
  filter(t>10) %>%
  group_by(simulation, model, fit) %>%
  summarise(Rt_coverage = 100*round(mean((TrueRt <= upper) & (TrueRt >= lower)), digits=3),
            Ct_coverage = 100*round(mean((Cases <= upper_cases) & (Cases >= lower_cases)), digits=3),
            score = signif(mean(crps), 3))



# Plot the raw data
plt_data_rw = ggplot(df %>% filter(simulation==CURRENT_SIM, model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=date, y=Cases), width = 0.8) +
  # scale_x_continuous(limits = c(0, 100)) +
  custom_theme +
  xlab("Date") + ylab("Reported cases") +
  ggtitle(labels[CURRENT_SIM])


# Plot posteriors
plt_post_rw = ggplot(df %>% filter(simulation==CURRENT_SIM)) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=date, y=mean, color=NameForColors), data=.%>%filter(t>10)) +
  geom_line(aes(x=date, y=TrueRt), color="black", lwd=0.3) +
  facet_grid(fit~model) +
  # coord_cartesian(ylim=c(0.2, 2.0)) +
  ylab("Reproduction number") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")


# Plot predictive posteriors
plt_predpost_rw = ggplot(df %>% filter(simulation==CURRENT_SIM)) +
  geom_ribbon(aes(x=date, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>10), alpha=0.5) +
  geom_line(aes(x=date, y=mean_cases, color=NameForColors), data=.%>%filter(t>10)) +
  geom_point(aes(x=date, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(simulation == CURRENT_SIM),
            aes(x = min(df$date), y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  # coord_cartesian(ylim=c(0,200)) +scal
  ylab("Reported cases") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")


# Plot parameter estimates
xlims = c(0, 1.0)

plt_param_rw_ef = ggplot(df_efpost %>% filter(simulation==CURRENT_SIM)) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dotted") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_rw_ee = ggplot(df_eepost %>% filter(simulation==CURRENT_SIM)) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("Posterior density") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(0,10), minor_breaks=NULL) +
  geom_vline(xintercept=7, color="black", linetype="dotted")

plt_param_rw = plot_grid(plt_param_rw_ee, plt_param_rw_ef, ncol=2)



plt = plot_grid(plot_grid(plt_data_rw, plt_param_rw, ncol=1, rel_heights=c(1, 0.7)),
                plt_post_rw,
                plt_predpost_rw,
                ncol=3) 
plt
ggsave(paste0("paper/figures/s5", CURRENT_SIM, ".png"), plt, dpi=300, width=25, height=8, units="cm")


