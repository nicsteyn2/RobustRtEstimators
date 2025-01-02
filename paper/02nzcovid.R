
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

perioda_st = as.Date("2021-08-29")
perioda_en = as.Date("2021-09-06")
periodb_st = as.Date("2021-09-10")
periodb_en = as.Date("2021-10-06")
periodc_st = as.Date("2021-11-13")
periodc_en = as.Date("2021-12-01")

# Load the data
df = read.csv("paper/outputs/02nzcovid.csv")
df_eepost = read.csv("paper/outputs/02nzcovid_epiestimposteriors.csv")
df_efpost = read.csv("paper/outputs/02nzcovid_epifilterposteriors.csv")

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
  group_by(simulation, model, fit) %>%
  summarise(Rt_coverage = 100*round(mean((TrueRt <= upper) & (TrueRt >= lower)), digits=3),
            Ct_coverage = 100*round(mean((Cases <= upper_cases) & (Cases >= lower_cases)), digits=3),
            score = signif(mean(crps), 3))



# Plot the raw data
plt_data = ggplot(df %>% filter(model=="EpiEstim", fit=="Marginalised")) +
  geom_rect(aes(xmin=perioda_st, xmax=perioda_en, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.4, data=.%>%filter(t==1)) +
  geom_rect(aes(xmin=periodb_st, xmax=periodb_en, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.4, data=.%>%filter(t==1)) +
  geom_rect(aes(xmin=periodc_st, xmax=periodc_en, ymin=-Inf, ymax=Inf), fill="gray", alpha=0.4, data=.%>%filter(t==1)) +
  annotate("text", x=perioda_st+1, y=210, label="A", hjust=0, size=3) +
  annotate("text", x=periodb_st+1, y=210, label="B", hjust=0, size=3) +
  annotate("text", x=periodc_st+1, y=210, label="C", hjust=0, size=3) +
  geom_col(aes(x=date, y=Cases), width = 0.8) +
  # scale_x_continuous(limits = c(0, 100)) +
  custom_theme +
  xlab("Date") + ylab("Reported cases") +
  ggtitle("")

# Plot posteriors
plt_post = ggplot(df) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_line(aes(x=date, y=TrueRt), color="black", lwd=0.3) +
  facet_grid(fit~model) +
  # coord_cartesian(ylim=c(0.2, 2.0)) +
  ylab("Reproduction number") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")

# Plot predictive posteriors
plt_predpost = ggplot(df) +
  geom_ribbon(aes(x=date, ymin=lower_cases, ymax=upper_cases, fill=NameForColors), data=.%>%filter(t>windin), alpha=0.5) +
  geom_line(aes(x=date, y=mean_cases, color=NameForColors), data=.%>%filter(t>windin)) +
  geom_point(aes(x=date, y=Cases), color="black", size=0.3) +
  geom_text(data = df_coverage,
            aes(x = min(df$date), y = Inf, label = paste0("Coverage = ", Ct_coverage, "%\nScore = ", score)),
            hjust = 0, vjust = 1.2, size = 3, color = "black") +
  facet_grid(fit~model) +
  coord_cartesian(ylim=c(0,300)) +
  ylab("Reported cases") + xlab("Date") +
  custom_theme +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  theme(legend.position="none")


# Plot parameter estimates
xlims = c(0, 0.2)

plt_param_ef = ggplot(df_efpost) +
  geom_line(aes(x=eta, y=p), color="#DE7139", lwd=1) +
  geom_vline(xintercept=0.1, color="black", linetype="dotted") +
  custom_theme +
  coord_cartesian(xlim=xlims) +
  xlab("eta") + ylab("")

plt_param_ee = ggplot(df_eepost) +
  geom_col(aes(x=k, y=p), fill="#006666") +
  custom_theme +
  xlab("k") + ylab("Posterior density") +
  scale_x_continuous(breaks=seq(1,10), labels=seq(1,10), limits=c(1,10), minor_breaks=NULL) +
  geom_vline(xintercept=7, color="black", linetype="dotted")

plt_param = plot_grid(plt_param_ee, plt_param_ef, ncol=2)



plt = plot_grid(plot_grid(plt_data, plt_param, ncol=1, rel_heights=c(1, 0.7)),
                plt_post,
                plt_predpost,
                ncol=3) 



# Zoomed-in Rt plots

plt_zoom_1 = ggplot(df %>% filter(date<=perioda_en, date>=perioda_st)) +
  geom_hline(yintercept=1, color="black") +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors, linetype=fit2), alpha=0.1) +
  geom_line(aes(x=date, y=upper, color=NameForColors, linetype=fit2)) +
  geom_line(aes(x=date, y=lower, color=NameForColors, linetype=fit2)) +
  geom_line(aes(x=date, y=mean, color=NameForColors, linetype=fit2)) +
  facet_wrap(~model, ncol=2) +
  custom_theme  +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  scale_x_date(minor_breaks="1 day") +
  theme(legend.position="none", plot.title = element_text(hjust = 0), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Period 1") +
  xlab("") + ylab("Rt")

plt_zoom_2 = ggplot(df %>% filter(date>=periodb_st, date<=periodb_en)) +
  geom_hline(yintercept=1, color="black") +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors, linetype=fit2), alpha=0.1) +
  geom_line(aes(x=date, y=upper, color=NameForColors, linetype=fit2)) +
  geom_line(aes(x=date, y=lower, color=NameForColors, linetype=fit2)) +
  geom_line(aes(x=date, y=mean, color=NameForColors, linetype=fit2)) +
  facet_wrap(~model, ncol=2) +
  custom_theme  +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  scale_x_date(minor_breaks="1 day") +
  theme(legend.position="none", plot.title = element_text(hjust = 0), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Period 2") +
  xlab("Date") + ylab("")

plt_zoom_3 = ggplot(df %>% filter(date>=periodc_st, date<=periodc_en)) +
  geom_hline(yintercept=1, color="black") +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=NameForColors, linetype=fit2), alpha=0.1) +
  geom_line(aes(x=date, y=upper, color=NameForColors, linetype=fit2)) +
  geom_line(aes(x=date, y=lower, color=NameForColors, linetype=fit2)) +
  geom_line(aes(x=date, y=mean, color=NameForColors, linetype=fit2)) +
  facet_wrap(~model, ncol=2) +
  custom_theme  +
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) +
  scale_x_date(minor_breaks="1 day") +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Period 3") +
  xlab("") + ylab("") +
  guides(color="none", fill="none", linetype=guide_legend(title=NULL, direction="horizontal"))

legend = get_legend(plt_zoom_3)

plt2 = plot_grid(plot_grid(plt_data, plt_param, ncol=1, rel_heights=c(1, 0.7)),
                 plt_post,
                 plt_predpost,
                 plt_zoom_1,
                 plt_zoom_2,
                 plt_zoom_3 + theme(legend.position="none"),
                 NULL,
                 legend,
                 NULL,
                 ncol=3,
                 rel_heights=c(0.7,0.5,0.05),
                 labels=c("A)", "B)", "C)", "D)", "E)", "F)")) 
plt2
ggsave(paste0("paper/figures/02nzcovid.png"), plt2, dpi=300, width=25, height=14, units="cm")


# Fetch numbers for table

# First window
df1 = df %>%
  filter(date >= perioda_st, date <= perioda_en) %>%
  group_by(model, fit) %>%
  summarise(lower_det = date[which(lower<1)[1]],
            mean_det = date[which(mean<1)[1]],
            upper_det = date[which(upper<1)[1]]) %>%
  ungroup() %>%
  mutate(name = paste0(model, fit)) %>%
  select(-model, -fit) %>%
  column_to_rownames(var="name") %>%
  t()

df2 = df %>%
  filter(date >= periodb_st, date <= periodb_en) %>%
  group_by(model, fit) %>%
  summarise(lower_det = date[which(lower>1)[1]],
            mean_det = date[which(mean>1)[1]],
            upper_det = date[which(upper>1)[1]]) %>%
  ungroup() %>%
  mutate(name = paste0(model, fit)) %>%
  select(-model, -fit) %>%
  column_to_rownames(var="name") %>%
  t()

df3 = df %>%
  filter(date >= periodc_st, date <= periodc_en) %>%
  group_by(model, fit) %>%
  summarise(lower_det = date[which(lower<1)[1]],
            mean_det = date[which(mean<1)[1]],
            upper_det = date[which(upper<1)[1]]) %>%
  ungroup() %>%
  mutate(name = paste0(model, fit)) %>%
  select(-model, -fit) %>%
  column_to_rownames(var="name") %>%
  t()

