## Appendix IV: alarm rate over the course of admission, by gestational age (NICU)

library(dplyr)
library(lubridate)
library(mgcv)
library(ggplot2)
library(ggthemes)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

UNIT <- "Neonatologie"
ga_levels <- c("Extremely preterm", "Very preterm", "Moderate to late preterm", "Term")

## ======== Alarm-active hours per admission, with the time since admission
subset <- data %>% semi_join(pat_char %>% filter(!is.na(GestAgeGroup)), by = c("patientID", "afdeling")) %>%
  filter(afdeling == UNIT)

subset <- subset %>%
  mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>%
  group_by(patientID, alarm_floor_hr) %>%
  summarise(alarms_per_hour = n(), .groups = 'drop') %>%
  left_join(., pat_char %>% filter(afdeling == UNIT) %>% select(patientID, GestAgeGroup, admission_start, admission_status), by = c("patientID" = "patientID")) %>%
  group_by(patientID) %>%
  arrange(alarm_floor_hr) %>%
  mutate(
    time_since_admission_hr = as.numeric(difftime(alarm_floor_hr, floor_date(admission_start, "1 hour"), units = "hours")),
    admission_period_hr = floor(time_since_admission_hr) + 1,
    admission_day = (admission_period_hr + 1) / 24,
    log_admission_day = log(admission_day)
  ) %>%
  ungroup() %>%
  filter(!is.na(GestAgeGroup)) %>%
  mutate(GestAgeGroup = factor(as.character(GestAgeGroup), levels = ga_levels)) %>%
  mutate(patientID = as.factor(patientID))

subset %>% group_by(patientID, admission_status) %>% summarise(first_hr = min(time_since_admission_hr), .groups = 'drop') %>%
  group_by(admission_status) %>% summarise(admissions = n(), min_hours = min(first_hr), median_hours = median(first_hr),
                                           q75_hours = quantile(first_hr, 0.75), max_hours = max(first_hr))
subset %>% group_by(GestAgeGroup) %>% summarise(admissions = n_distinct(patientID))

subset <- subset %>%
  group_by(GestAgeGroup) %>%
  filter(n_distinct(patientID) > 25) %>%
  ungroup() %>%
  mutate(GestAgeGroup = droplevels(GestAgeGroup))
print(table(subset$GestAgeGroup))

admission_duration_p90 <- subset %>%
  group_by(patientID, GestAgeGroup) %>%
  summarise(max_log_day = max(log_admission_day, na.rm = TRUE), .groups = 'drop') %>%
  group_by(GestAgeGroup) %>%
  summarise(p90_log_duration = quantile(max_log_day, 0.90, na.rm = TRUE), .groups = 'drop')

admission_duration_p975 <- subset %>%
  group_by(patientID, GestAgeGroup) %>%
  summarise(max_log_day = max(log_admission_day, na.rm = TRUE), .groups = 'drop') %>%
  group_by(GestAgeGroup) %>%
  summarise(p975_log_duration = quantile(max_log_day, 0.975, na.rm = TRUE), .groups = 'drop')

## ======== Spline basis size (k) per gestational-age group
group_sample_sizes <- subset %>%
  group_by(GestAgeGroup) %>%
  summarise(
    n_observations = n(),
    n_patients = n_distinct(patientID),
    .groups = 'drop'
  ) %>%
  mutate(
    spline_df = case_when(
      n_patients < 40 ~ 4,                           # few admissions: minimum for tp splines
      n_patients < 70 & n_observations < 4000 ~ 8,   # moderate admissions and hours
      TRUE ~ 15                                      # data-rich groups
    )
  )
print(group_sample_sizes)

group_levels <- levels(subset$GestAgeGroup)
k_by_group <- setNames(group_sample_sizes$spline_df, group_sample_sizes$GestAgeGroup)

## ======== Models: negative binomial GAM, log link (mgcv::bam)
overall_model <- mgcv::bam(
  alarms_per_hour ~ s(log_admission_day, bs = "tp", k = 15) +
    s(patientID, bs = "re") +
    s(patientID, log_admission_day, bs = "re"),
  family = nb(),
  data = subset,
  method = "fREML",
  discrete = TRUE
)

summary(overall_model)

## ======== Population-average predictions (random effects excluded)
predict_typical_patient <- function(fit, newdata) {
  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE,
                exclude = c("s(patientID)", "s(log_admission_day,patientID)"),
                newdata.guaranteed = TRUE)
  newdata %>%
    mutate(
      predicted = exp(pr$fit),
      conf.low = exp(pr$fit - 1.96 * pr$se.fit),
      conf.high = exp(pr$fit + 1.96 * pr$se.fit)
    )
}

placeholder_patient <- overall_model$model$patientID[1]
log_day_grid <- seq(min(subset$log_admission_day), max(subset$log_admission_day), length.out = 200)

subset %>% group_by(GestAgeGroup) %>% summarise(median_aph = median(alarms_per_hour), mean_aph = round(mean(alarms_per_hour), 1))

set.seed(20261031)
preds_groups <- NULL
kc <- mgcv::k.check(overall_model)
sm <- summary(overall_model)
model_summaries <- c("======== Overall Average (k = 15)", capture.output(sm),
                     "-- k.check", capture.output(kc), "")
model_diagnostics <- data.frame(department = "NICU", model = "Overall Average", k = 15,
                                n_alarm_active_hours = nrow(overall_model$model),
                                n_admissions = n_distinct(overall_model$model$patientID),
                                k_prime = kc["s(log_admission_day)", "k'"], edf = kc["s(log_admission_day)", "edf"],
                                k_index = kc["s(log_admission_day)", "k-index"], k_check_p = kc["s(log_admission_day)", "p-value"],
                                theta = overall_model$family$getTheta(TRUE))
for (ga in group_levels) {
  subset_ga <- subset %>% filter(GestAgeGroup == ga) %>% mutate(patientID = droplevels(patientID))
  k_value <- k_by_group[[ga]]
  group_model <- mgcv::bam(
    alarms_per_hour ~ s(log_admission_day, bs = "tp", k = k_value) +
      s(patientID, bs = "re") +
      s(patientID, log_admission_day, bs = "re"),
    family = nb(),
    data = subset_ga,
    method = "fREML",
    discrete = TRUE
  )
  nd <- data.frame(log_admission_day = log_day_grid, patientID = group_model$model$patientID[1])
  preds_groups <- bind_rows(preds_groups,
                            predict_typical_patient(group_model, nd) %>%
                              mutate(admission_day = exp(log_admission_day), GestAgeGroup = ga))
  kc <- mgcv::k.check(group_model)
  sm <- summary(group_model)
  model_summaries <- c(model_summaries, paste0("======== ", ga, " (k = ", k_value, ")"), capture.output(sm),
                       "-- k.check", capture.output(kc), "")
  model_diagnostics <- bind_rows(model_diagnostics,
                                 data.frame(department = "NICU", model = ga, k = k_value,
                                            n_alarm_active_hours = nrow(group_model$model),
                                            n_admissions = n_distinct(group_model$model$patientID),
                                            k_prime = kc["s(log_admission_day)", "k'"], edf = kc["s(log_admission_day)", "edf"],
                                            k_index = kc["s(log_admission_day)", "k-index"], k_check_p = kc["s(log_admission_day)", "p-value"],
                                            theta = group_model$family$getTheta(TRUE)))
  rm(group_model, subset_ga)
}

preds_groups_with_boundary <- preds_groups %>%
  select(admission_day, predicted, conf.low, conf.high, GestAgeGroup) %>%
  left_join(admission_duration_p90 %>% mutate(GestAgeGroup = as.character(GestAgeGroup)), by = "GestAgeGroup") %>%
  left_join(admission_duration_p975 %>% mutate(GestAgeGroup = as.character(GestAgeGroup)), by = "GestAgeGroup") %>%
  filter(log(admission_day) <= p975_log_duration) %>%
  mutate(
    p90_boundary_days = exp(p90_log_duration),
    line_typing = ifelse(admission_day > p90_boundary_days, "dotted", "solid")
  ) %>%
  select(-p90_log_duration, -p975_log_duration)

overall_p90_boundary_log <- quantile(subset %>% group_by(patientID) %>% summarise(max_log_day = max(log_admission_day, na.rm = TRUE)) %>% pull(max_log_day), 0.90, na.rm = TRUE)
overall_p975_boundary_log <- quantile(subset %>% group_by(patientID) %>% summarise(max_log_day = max(log_admission_day, na.rm = TRUE)) %>% pull(max_log_day), 0.975, na.rm = TRUE)

nd_overall <- data.frame(log_admission_day = log_day_grid, patientID = placeholder_patient)
preds_overall <- predict_typical_patient(overall_model, nd_overall) %>%
  filter(log_admission_day <= overall_p975_boundary_log) %>%
  mutate(admission_day = exp(log_admission_day),
         GestAgeGroup = "Overall Average",
         line_typing = ifelse(log_admission_day > overall_p90_boundary_log, "dotted", "solid"))

preds_combined <- bind_rows(
  preds_groups_with_boundary,
  preds_overall %>% select(admission_day, predicted, conf.low, conf.high, GestAgeGroup, line_typing)
) %>%
  mutate(GestAgeGroup = factor(GestAgeGroup, levels = c("Overall Average", group_levels)))
write.csv(preds_combined, "output//appendices//appendixIV_GA_pred.csv", row.names = FALSE)
write.csv(group_sample_sizes, "output//appendices//appendixIV_GA_k.csv", row.names = FALSE)
writeLines(model_summaries, "output//appendices//appendixIV_GA_model_summary.txt")
write.csv(model_diagnostics, "output//appendices//appendixIV_GA_diagnostics.csv", row.names = FALSE)

preds_combined <- bind_rows(preds_combined,
                            preds_combined %>% filter(line_typing == "solid") %>% group_by(GestAgeGroup) %>%
                              slice_max(admission_day, n = 1) %>% ungroup() %>% mutate(line_typing = "dotted"))

## ======== Plots
theme_Publication <- function(base_size=14, base_family="sans") {
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(fill = "white", colour = NA),
               plot.background = element_rect(fill = "white", colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line.x = element_line(colour="black"),
               axis.line.y = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
       ))
}
ga_colours <- c("Extremely preterm"        = "#E64B35",
                "Very preterm"             = "#F39B7F",
                "Moderate to late preterm" = "#4DBBD5",
                "Term"                     = "#00A087",
                "Overall Average"          = "#2c2c2c")

plot <- preds_combined %>%
  ggplot(aes(x = admission_day, y = predicted, color = GestAgeGroup, group = interaction(GestAgeGroup, line_typing))) +
  geom_ribbon(data = . %>% filter(line_typing == "solid"),
              aes(x = admission_day, ymin = conf.low, ymax = conf.high, fill = GestAgeGroup), alpha = 0.2, color = NA) +
  geom_ribbon(data = . %>% filter(line_typing == "dotted"),
              aes(x = admission_day, ymin = conf.low, ymax = conf.high), fill = "grey", color = NA, alpha = 0.2) +
  geom_line(data = . %>% filter(line_typing == "solid"), aes(x = admission_day, y = predicted), linetype = "solid", linewidth = 1.2) +
  geom_line(data = . %>% filter(line_typing == "dotted"), aes(x = admission_day, y = predicted), linetype = "dotted", linewidth = 1.2) +
  scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(xlim = c(0.1, 100), ylim = c(0, 100)) +
  theme_Publication() +
  scale_color_manual(name = "Gestational age group", values = ga_colours) +
  scale_fill_manual(name = "Gestational age group", values = ga_colours) +
  xlab("Days since admission (log scale)") +
  ylab("Predicted alarms per alarm-active hour") +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  annotation_logticks(sides = "b") +
  guides(color = guide_legend(override.aes = list(linetype = "solid")))
ggsave("output//appendices//appendixIV_GA_all_groups.png", plot, width = 40, height = 25, units = "cm", dpi = 300)

ga_headers <- c("Overall Average"          = "Overall Average (all gestational ages)",
                "Extremely preterm"        = "Extremely preterm (<28w0d)",
                "Very preterm"             = "Very preterm (28w0d - 31w6d)",
                "Moderate to late preterm" = "Moderate to late preterm (32w0d - 36w6d)",
                "Term"                     = "Term (≥37w0d)")
plot <- preds_combined %>%
  ggplot(aes(x = admission_day, y = predicted, group = interaction(GestAgeGroup, line_typing))) +
    geom_ribbon(data = . %>% filter(line_typing == "solid"),
                aes(x = admission_day, ymin = conf.low, ymax = conf.high),
                fill = "#386cb0", color = "#386cb0", alpha = 0.2) +
    geom_ribbon(data = . %>% filter(line_typing == "dotted"),
                aes(x = admission_day, ymin = conf.low, ymax = conf.high),
                fill = "grey", color = "grey", alpha = 0.2) +
    geom_line(data = . %>% filter(line_typing == "solid"),
              aes(x = admission_day, y = predicted),
              color = "#27456c", linetype = "solid", linewidth = 1.0) +
    geom_line(data = . %>% filter(line_typing == "dotted"),
              aes(x = admission_day, y = predicted), color = "grey", linetype = "dotted", linewidth = 1.0) +
    scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
    scale_y_continuous(breaks = seq(0, 100, 25)) +
    coord_cartesian(xlim = c(0.1, 100), ylim = c(0, 100)) +
    theme_Publication() +
    facet_wrap(~ GestAgeGroup, ncol = 3, labeller = as_labeller(ga_headers)) +
    xlab("Days since admission (log scale)") +
    ylab("Predicted alarms per alarm-active hour") +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
    annotation_logticks(sides = "b")
plot <- ggplotGrob(plot)
overall_strip <- which(plot$layout$name == "strip-t-1-1")
plot$grobs[[overall_strip]]$grobs[[1]]$children[[1]]$gp$fill <- "#2c2c2c"
plot$grobs[[overall_strip]]$grobs[[1]]$children[[1]]$gp$col <- "#2c2c2c"
plot$grobs[[overall_strip]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "white"
ggsave("output//appendices//appendixIV_GA_facet.png", plot, width = 40, height = 25, units = "cm", dpi = 300)
