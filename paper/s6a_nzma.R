
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
    plot.title = element_text(hjust = 0.5, size=10)#,
    # panel.grid = element_blank()
  )

windin=3
xlims = c(0, 0.4)

# Load the data
df = read.csv("paper/outputs/s6a_nzmovingaverage.csv")
df_eepost = read.csv("paper/outputs/s6a_nzmovingaverage_epiestimposteriors.csv")
df_efpost = read.csv("paper/outputs/s6a_nzmovingaverage_epifilterposteriors.csv")

# Clean the data
df$fit = ifelse(grepl("Conditional", df$fit), "Default", df$fit)
df = df %>%
  mutate(NameForColors=paste0(model, fit),
         date = as.Date("2021-08-18") - 1 + t) %>%
  filter(date <= as.Date("2021-12-01")) %>%
  mutate(fit2 = factor(fit, levels=c("Marginalised", "Default")),
         fit=factor(fit, levels=c("Default", "Marginalised")))

# Append colors
colors <- c(EpiEstimMarginalised="#006666", EpiEstimDefault="#20B2AA",
            EpiFilterMarginalised="#DE7139", EpiFilterDefault="#E5785A")

# Calculate coverage statistics
df_coverage = df %>%
  filter(t>windin) %>%
  group_by(movingaverage, model, fit) %>%
  summarise(Rt_coverage = 100*round(mean((TrueRt <= upper) & (TrueRt >= lower)), digits=3),
            Ct_coverage = 100*round(mean((Cases <= upper_cases) & (Cases >= lower_cases)), digits=3),
            score = signif(mean(crps), 3))





# Plot results from raw data
plt_data_raw = ggplot(df %>% filter(movingaverage=="Raw", model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=date, y=Cases), width = 0.8) +
  custom_theme +
  xlab("Date") + ylab("Reported cases") +
  ggtitle("Raw data")

plt_post_raw = ggplot(df %>% filter(movingaverage=="Raw")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_line(aes(x=date, y=TrueRt), color="black", lwd=0.3) +
  facet_grid(fit~model) +
  ylab("Reproduction number") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_predpost_raw = ggplot(df %>% filter(movingaverage=="Raw")) +
  geom_ribbon(aes(x=date, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean_cases, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_point(aes(x=date, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(movingaverage == "Raw"),
            aes(x = min(df$date), y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0,300)) +
  ylab("Reported cases") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_param_rw_ef = ggplot(df_efpost %>% filter(movingaverage=="Raw")) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dotted") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_rw_ee = ggplot(df_eepost %>% filter(movingaverage=="Raw")) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("Posterior density") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(1,10), minor_breaks=NULL) +
  geom_vline(xintercept=7, color="black", linetype="dotted")

plt_param_raw = plot_grid(plt_param_rw_ee, plt_param_rw_ef, ncol=2)

plt_raw = plot_grid(plt_data_raw, plt_param_raw, plt_post_raw, plt_predpost_raw, ncol=1 , rel_heights=c(0.6, 0.4, 1, 1))


# Plot results from 5-day ma
plt_data_5 = ggplot(df %>% filter(movingaverage=="5-day", model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=date, y=Cases), width = 0.8) +
  custom_theme +
  xlab("Date") + ylab("Reported cases") +
  ggtitle("5-day MA")

plt_post_5 = ggplot(df %>% filter(movingaverage=="5-day")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_line(aes(x=date, y=TrueRt), color="black", lwd=0.3) +
  facet_grid(fit~model) +
  ylab("Reproduction number") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_predpost_5 = ggplot(df %>% filter(movingaverage=="5-day")) +
  geom_ribbon(aes(x=date, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean_cases, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_point(aes(x=date, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(movingaverage == "5-day"),
            aes(x = min(df$date), y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0,300)) +
  ylab("Reported cases") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_param_5_ef = ggplot(df_efpost %>% filter(movingaverage=="5-day")) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dotted") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_5_ee = ggplot(df_eepost %>% filter(movingaverage=="5-day")) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("Posterior density") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(1,10), minor_breaks=NULL) +
  geom_vline(xintercept=7, color="black", linetype="dotted")

plt_param_5 = plot_grid(plt_param_5_ee, plt_param_5_ef, ncol=2)

plt_5 = plot_grid(plt_data_5, plt_param_5, plt_post_5, plt_predpost_5, ncol=1 , rel_heights=c(0.6, 0.4, 1, 1))



# Plot results from 10-day ma
plt_data_10 = ggplot(df %>% filter(movingaverage=="10-day", model=="EpiEstim", fit=="Marginalised")) +
  geom_col(aes(x=date, y=Cases), width = 0.8) +
  custom_theme +
  xlab("Date") + ylab("Reported cases") +
  ggtitle("10-day MA")

plt_post_10 = ggplot(df %>% filter(movingaverage=="10-day")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_line(aes(x=date, y=TrueRt), color="black", lwd=0.3) +
  facet_grid(fit~model) +
  ylab("Reproduction number") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_predpost_10 = ggplot(df %>% filter(movingaverage=="10-day")) +
  geom_ribbon(aes(x=date, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean_cases, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_point(aes(x=date, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage %>% filter(movingaverage == "10-day"),
            aes(x = min(df$date), y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0,300)) +
  ylab("Reported cases") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

plt_param_10_ef = ggplot(df_efpost %>% filter(movingaverage=="10-day")) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dotted") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_10_ee = ggplot(df_eepost %>% filter(movingaverage=="10-day")) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("Posterior density") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(1,10), minor_breaks=NULL) +
  geom_vline(xintercept=7, color="black", linetype="dotted")

plt_param_10 = plot_grid(plt_param_10_ee, plt_param_10_ef, ncol=2)

plt_10 = plot_grid(plt_data_10, plt_param_10, plt_post_10, plt_predpost_10, ncol=1 , rel_heights=c(0.6, 0.4, 1, 1))



plt = plot_grid(plt_raw, plt_5, plt_10, ncol=3) 
plt
ggsave(paste0("paper/figures/s6a_nzma.png"), plt, dpi=300, width=25, height=21, units="cm")
