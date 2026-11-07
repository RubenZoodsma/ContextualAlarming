## Figure 1 and Appendix IV: alarm rate over the course of admission, by reason of admission

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

## Department: set UNIT
UNIT <- "Neonatologie"
# UNIT <- "Intensive Care Kinderen"
UNIT_SHORT <- c("Neonatologie" = "NICU", "Intensive Care Kinderen" = "PICU")[[UNIT]]

## ======== Alarm-active hours per admission, with the time since admission

subset <- data %>% semi_join(pat_char %>% filter(!is.na(diag_group)), by = c("patientID", "afdeling")) %>%
  filter(afdeling == UNIT)

subset <- subset %>%
  mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>%
  group_by(patientID, alarm_floor_hr) %>%
  summarise(alarms_per_hour = n(), .groups = 'drop') %>%
  left_join(., pat_char %>% filter(afdeling == UNIT) %>% select(patientID, diag_group, admission_start, admission_status), by = c("patientID" = "patientID")) %>%
  group_by(patientID) %>%
  arrange(alarm_floor_hr) %>%
  mutate(
    time_since_admission_hr = as.numeric(difftime(alarm_floor_hr, floor_date(admission_start, "1 hour"), units = "hours")),
    admission_period_hr = floor(time_since_admission_hr) + 1,
    admission_day = (admission_period_hr + 1) / 24,
    log_admission_day = log(admission_day)
  ) %>%
  ungroup() %>%
  filter(!is.na(diag_group)) %>%
  mutate(diag_group = as.factor(diag_group)) %>%
  mutate(patientID = as.factor(patientID))

subset %>% group_by(patientID, admission_status) %>% summarise(first_hr = min(time_since_admission_hr), .groups = 'drop') %>%
  group_by(admission_status) %>% summarise(admissions = n(), min_hours = min(first_hr), median_hours = median(first_hr),
                                           q75_hours = quantile(first_hr, 0.75), max_hours = max(first_hr))

subset <- subset %>%
  group_by(diag_group) %>%
  filter(n_distinct(patientID) > 25) %>%
  ungroup() %>%
  mutate(diag_group = droplevels(diag_group))
print(table(subset$diag_group))

admission_duration_p90 <- subset %>%
  group_by(patientID, diag_group) %>%
  summarise(max_log_day = max(log_admission_day, na.rm = TRUE), .groups = 'drop') %>%
  group_by(diag_group) %>%
  summarise(p90_log_duration = quantile(max_log_day, 0.90, na.rm = TRUE), .groups = 'drop')

admission_duration_p975 <- subset %>%
  group_by(patientID, diag_group) %>%
  summarise(max_log_day = max(log_admission_day, na.rm = TRUE), .groups = 'drop') %>%
  group_by(diag_group) %>%
  summarise(p975_log_duration = quantile(max_log_day, 0.975, na.rm = TRUE), .groups = 'drop')

## ======== Spline basis size (k) per reason of admission
group_sample_sizes <- subset %>%
  group_by(diag_group) %>%
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

group_levels <- levels(subset$diag_group)
k_by_group <- setNames(group_sample_sizes$spline_df, group_sample_sizes$diag_group)

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

subset %>% group_by(diag_group) %>% summarise(median_aph = median(alarms_per_hour), mean_aph = round(mean(alarms_per_hour), 1)) %>% arrange(desc(mean_aph))
if (UNIT == "Intensive Care Kinderen") { ## Highest and lowest intensity groups for pediatric intensive care
highest_intensity_group <- 'Cardiac'
lowest_intensity_group <- 'Traumatic'
} else {    ## Highest and lowest intensity groups for Neonatal intensive care
highest_intensity_group <- 'Asphyxia'
lowest_intensity_group <- 'Convulsions (suspected)'
}

set.seed(20261031)
preds_groups <- NULL
model_summaries <- c("======== Overall Average (k = 15)", capture.output(summary(overall_model)),
                     "-- k.check", capture.output(mgcv::k.check(overall_model)), "")
for (dg in group_levels) {
  subset_dg <- subset %>% filter(diag_group == dg) %>% mutate(patientID = droplevels(patientID))
  k_value <- k_by_group[[dg]]
  group_model <- mgcv::bam(
    alarms_per_hour ~ s(log_admission_day, bs = "tp", k = k_value) +
      s(patientID, bs = "re") +
      s(patientID, log_admission_day, bs = "re"),
    family = nb(),
    data = subset_dg,
    method = "fREML",
    discrete = TRUE
  )
  nd <- data.frame(log_admission_day = log_day_grid, patientID = group_model$model$patientID[1])
  preds_groups <- bind_rows(preds_groups,
                            predict_typical_patient(group_model, nd) %>%
                              mutate(admission_day = exp(log_admission_day), line_type = dg, diag_group = dg))
  model_summaries <- c(model_summaries, paste0("======== ", dg, " (k = ", k_value, ")"), capture.output(summary(group_model)),
                       "-- k.check", capture.output(mgcv::k.check(group_model)), "")
  rm(group_model, subset_dg)
}

preds_groups_with_boundary <- preds_groups %>%
  select(admission_day, predicted, conf.low, conf.high, line_type, diag_group) %>%
  left_join(admission_duration_p90, by = "diag_group") %>%
  left_join(admission_duration_p975, by = "diag_group") %>%
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
         line_type = "Overall Average",
         diag_group = "Overall Average",
         line_typing = ifelse(log_admission_day > overall_p90_boundary_log, "dotted", "solid"))

preds_combined <- bind_rows(
  preds_groups_with_boundary %>% mutate(diag_group = as.character(diag_group)),
  preds_overall %>% select(admission_day, predicted, conf.low, conf.high, line_type, diag_group, line_typing)
) %>%
  mutate(
    line_type = factor(line_type, levels = c("Overall Average", lowest_intensity_group, highest_intensity_group)),
    diag_group = factor(diag_group, levels = c("Overall Average", group_levels))
  )
write.csv(preds_combined, paste0("output//primaryOutput//figure1_pred_", UNIT_SHORT, ".csv"), row.names = FALSE)
write.csv(group_sample_sizes, paste0("output//primaryOutput//figure1_k_", UNIT_SHORT, ".csv"), row.names = FALSE)
writeLines(model_summaries,
           paste0("output//primaryOutput//figure1_model_summary_", UNIT_SHORT, ".txt"))

preds_combined <- bind_rows(preds_combined,
                            preds_combined %>% filter(line_typing == "solid") %>% group_by(diag_group) %>%
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

plot <- preds_combined %>%
  filter(diag_group %in% c(highest_intensity_group, lowest_intensity_group, "Overall Average")) %>%
  ggplot(aes(x = admission_day, y = predicted, color = line_type, group = interaction(line_type, line_typing))) +
    geom_ribbon(data = . %>% filter(line_typing == "solid"),
                aes(x = admission_day, ymin = conf.low, ymax = conf.high),
                fill = "#386cb0", color = "#386cb0", alpha = 0.2) +
    geom_ribbon(data = . %>% filter(line_typing == "dotted"),
                aes(x = admission_day, ymin = conf.low, ymax = conf.high),
                fill = "grey", color = "grey", alpha = 0.2) +
    geom_line(data = . %>% filter(line_typing == "solid"),
              aes(x = admission_day, y = predicted),
              linetype = "solid", linewidth = 1.0) +
    geom_line(data = . %>% filter(line_typing == "dotted"),
              aes(x = admission_day, y = predicted), linetype = "dotted", linewidth = 1.0) +
    scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
    scale_y_continuous(breaks = seq(0, 100, 25)) +
    coord_cartesian(xlim = c(0.1, 100), ylim = c(0, 100)) +
    theme_Publication() +
    scale_color_manual(name = "Reason of admission", values = setNames(
      c("#2c2c2c", "#386cb0", "#d32f2f"),
      c("Overall Average", lowest_intensity_group, highest_intensity_group)
    )) +
    xlab("Days since admission (log scale)") +
    ylab("Predicted alarms per alarm-active hour") +

    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
    annotation_logticks(sides = "b") +
    guides(color = guide_legend(override.aes = list(linetype = "solid")))
ggsave(paste0("output//primaryOutput//figure1_key_", UNIT_SHORT, ".png"), plot, width = 22, height = 30, units = "cm", dpi = 300)

if (UNIT_SHORT == "PICU") {
  dept_limit <- 50
  dept_granularity <- 10
} else {
  dept_limit <- 100
  dept_granularity <- 25
}
plot <- preds_combined %>%
  ggplot(aes(x = admission_day, y = predicted, group = interaction(line_type, line_typing))) +
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
    scale_y_continuous(breaks = seq(0, dept_limit, dept_granularity)) +
    coord_cartesian(xlim = c(0.1, 100), ylim = c(0, dept_limit)) +
    theme_Publication() +
    facet_wrap(~ diag_group, ncol = 3) +
    xlab("Days since admission (log scale)") +
    ylab("Predicted alarms per alarm-active hour") +

    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
    annotation_logticks(sides = "b")

plot <- ggplotGrob(plot)
overall_strip <- which(plot$layout$name == "strip-t-1-1")
plot$grobs[[overall_strip]]$grobs[[1]]$children[[1]]$gp$fill <- "#2c2c2c"
plot$grobs[[overall_strip]]$grobs[[1]]$children[[1]]$gp$col <- "#2c2c2c"
plot$grobs[[overall_strip]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "white"
ggsave(paste0("output//appendices//appendixIV_diagnosis_", UNIT_SHORT, ".png"), plot, width = 40, height = 25, units = "cm", dpi = 300)
