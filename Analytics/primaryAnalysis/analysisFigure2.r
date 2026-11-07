## Figure 2: variation in alarm rate across contextual determinants, PICU and NICU

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

## ======== Department average and reason of admission
admitted <- pat_char %>% ungroup() %>%
  mutate(start = floor_date(pmax(admission_start, STUDY_START), "1 hour"), end = pmin(admission_end, STUDY_END)) %>% filter(end > start) %>%
  reframe(alarm_floor_hr = seq(start, end - 1, by = "1 hour"), .by = c(patientID, afdeling)) %>%
  left_join(data %>% mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>% count(patientID, afdeling, alarm_floor_hr, name = "alarms_per_hour"),
            by = c("patientID", "afdeling", "alarm_floor_hr")) %>%
  mutate(alarms_per_hour = coalesce(alarms_per_hour, 0L)) %>%
  left_join(pat_char %>% select(patientID, afdeling, diag_group), by = c("patientID", "afdeling"))

group_summary <- bind_rows(admitted, admitted %>% mutate(diag_group = "Total")) %>%
  group_by(afdeling, diag_group) %>%
  summarise(
    cases = n_distinct(patientID),
    median_active = median(alarms_per_hour[alarms_per_hour > 0]),
    .groups = 'drop')
print(group_summary, n = Inf)
admitted <- admitted %>% filter(alarms_per_hour > 0)

unit_average <- group_summary %>% filter(diag_group == "Total") %>% select(afdeling, average = median_active)
reason_rows <- group_summary %>% filter(diag_group != "Total", !is.na(diag_group), cases > 25) %>%
  group_by(afdeling) %>%
  summarise(category = "Reason of admission",
            min_label = diag_group[which.min(median_active)], minimum = min(median_active),
            max_label = diag_group[which.max(median_active)], maximum = max(median_active),
            min_group = diag_group[which.min(median_active)], max_group = diag_group[which.max(median_active)]) %>%
  left_join(unit_average, by = "afdeling")

## ======== Active alarm sources (3-hour backward window)
hours_active <- data %>% filter(!is.na(patientID)) %>% ## alarms without an admission have no alarm-active hours
  mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>%
  group_by(afdeling, patientID, alarm_floor_hr) %>%
  summarise(BP_alarm = "abp" %in% Parameter, Vasc_alarm = "vasc_pressure" %in% Parameter,
            Vent_alarm = "Ventilator" %in% Device, .groups = 'drop')
component_e <- hours_active %>% filter(BP_alarm | Vasc_alarm | Vent_alarm) %>%
  transmute(afdeling, patientID, e_s = as.numeric(alarm_floor_hr) / 3600, BP_alarm, Vasc_alarm, Vent_alarm)
target <- hours_active %>%
  transmute(afdeling, patientID, alarm_floor_hr,
            lo = as.numeric(alarm_floor_hr) / 3600 - 3, hi = as.numeric(alarm_floor_hr) / 3600 - 1) ## [h - 3, h)
monit.data <- target %>%
  left_join(target %>% inner_join(component_e, by = join_by(afdeling, patientID, lo <= e_s, hi >= e_s)) %>%
              group_by(afdeling, patientID, alarm_floor_hr) %>%
              summarise(BP = any(BP_alarm), Vasc = any(Vasc_alarm), Vent = any(Vent_alarm), .groups = 'drop'),
            by = c("afdeling", "patientID", "alarm_floor_hr")) %>%
mutate(BP = coalesce(BP, FALSE), Vasc = coalesce(Vasc, FALSE), Vent = coalesce(Vent, FALSE)) %>% ## no history -> regular
mutate(selection = case_when(
    !BP & !Vasc & !Vent ~ "regular",
    BP & !Vasc & !Vent ~ "regular + BP",
    !BP & Vasc & !Vent ~ "regular + Vasc",
    !BP & !Vasc & Vent ~ "regular + Vent",
    BP & Vasc & !Vent ~ "regular + BP + Vasc",
    BP & !Vasc & Vent ~ "regular + BP + Vent",
    !BP & Vasc & Vent ~ "regular + Vasc + Vent",
    BP & Vasc & Vent ~ "regular + BP + Vasc + Vent"
)) %>%
left_join(. ,
    data %>% filter(!is.na(patientID)) %>%
    mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")), by = c("patientID", "alarm_floor_hr", "afdeling")) %>% ## the hour's own alarms
    ungroup() %>% select(afdeling, patientID, alarm_floor_hr, selection)
rm(hours_active, component_e, target)

monit.hours <- monit.data %>%
  group_by(afdeling, selection, patientID, alarm_floor_hr) %>%
  summarise(alarms_per_hour = n(), .groups = 'drop')
monit.stats <- monit.hours %>%
  group_by(afdeling, selection) %>%
  summarise(alarm_active_hours = n(), median_alarms_per_hour = median(alarms_per_hour), .groups = 'drop')
print(monit.stats, n = Inf)
rm(monit.data)

monitoring_rows <- monit.stats %>% filter(alarm_active_hours > 1000) %>%
  group_by(afdeling) %>%
  summarise(category = "Monitoring intensity",
            min_label = "Regular", minimum = median_alarms_per_hour[selection == "regular"],
            max_label = selection[which.max(median_alarms_per_hour)], maximum = max(median_alarms_per_hour),
            min_group = "regular", max_group = selection[which.max(median_alarms_per_hour)]) %>%
  left_join(unit_average, by = "afdeling")

## ======== Gestational age (NICU only)
alarm_data_by_ga <- data %>% inner_join(pat_char %>% filter(afdeling == "Neonatologie" & !is.na(GestAgeGroup)) %>% select(patientID, afdeling, GestAgeGroup), by = c("patientID", "afdeling"))
ga.hours <- alarm_data_by_ga %>%
  mutate(alarm_start_floor_hr = floor_date(alarm_start, "1 hour")) %>%
  group_by(patientID, alarm_start_floor_hr, GestAgeGroup) %>%
  summarise(alarms_per_hour = n(), .groups = 'drop')
GA_Analysis <- ga.hours %>%
  group_by(GestAgeGroup) %>%
  summarise(n_patients = n_distinct(patientID), median_aph = median(alarms_per_hour), .groups = 'drop')
print(GA_Analysis)
rm(alarm_data_by_ga)

ga_rows <- GA_Analysis %>% filter(n_patients > 25) %>%
  summarise(afdeling = "Neonatologie", category = "Gestational age",
            min_label = as.character(GestAgeGroup[which.min(median_aph)]), minimum = min(median_aph),
            max_label = as.character(GestAgeGroup[which.max(median_aph)]), maximum = max(median_aph),
            min_group = as.character(GestAgeGroup[which.min(median_aph)]), max_group = as.character(GestAgeGroup[which.max(median_aph)])) %>%
  left_join(unit_average, by = "afdeling")

## ======== Time of day (15-minute resolution)
alarm_bins <- data %>% semi_join(pat_char, by = c("patientID", "afdeling")) %>%
  count(afdeling, patientID, bin_start = floor_date(alarm_start, "15 minutes"), name = "alarms")

windows <- bind_rows(
  alarm_bins %>% mutate(window_start = bin_start),
  alarm_bins %>% mutate(window_start = bin_start - minutes(15)),
  alarm_bins %>% mutate(window_start = bin_start - minutes(30)),
  alarm_bins %>% mutate(window_start = bin_start - minutes(45))) %>%
  filter(window_start >= STUDY_START & window_start + hours(1) <= STUDY_END) %>%
  group_by(afdeling, patientID, window_start) %>%
  summarise(alarms_per_hour = sum(alarms), .groups = 'drop')
rm(alarm_bins)

rate_over_day <- windows %>%
  mutate(time_of_day = hour(window_start) + minute(window_start) / 60) %>%
  group_by(afdeling, time_of_day) %>%
  summarise(median_alarms_per_hour = median(alarms_per_hour), .groups = 'drop') %>%
  mutate(window_mid = (time_of_day + 0.5) %% 24,
         window_label = paste0(sprintf("%02d:%02d", floor(time_of_day), round((time_of_day %% 1) * 60)), "-",
                               sprintf("%02d:%02d", floor((time_of_day + 1) %% 24), round((time_of_day %% 1) * 60))))

time_rows <- rate_over_day %>%
  group_by(afdeling) %>%
  summarise(category = "Time of day",
            min_bin = window_mid[which.min(median_alarms_per_hour)], min_label = window_label[which.min(median_alarms_per_hour)], minimum = min(median_alarms_per_hour),
            max_bin = window_mid[which.max(median_alarms_per_hour)], max_label = window_label[which.max(median_alarms_per_hour)], maximum = max(median_alarms_per_hour),
            min_group = as.character(time_of_day[which.min(median_alarms_per_hour)]), max_group = as.character(time_of_day[which.max(median_alarms_per_hour)])) %>%
  left_join(unit_average, by = "afdeling")

## ======== Figure values
values <- bind_rows(monitoring_rows, time_rows, reason_rows, ga_rows) %>%
  mutate(department = c("Intensive Care Kinderen" = "PICU", "Neonatologie" = "NICU")[afdeling],
         minus_pct = round((average - minimum) / average * 100),
         plus_pct = round((maximum - average) / average * 100)) %>%
  select(category, department, min_label, minimum, average, max_label, maximum, minus_pct, plus_pct, min_bin, max_bin)
print(values, width = Inf)
write.csv(values, "output//primaryOutput//figure2_summary_values.csv", row.names = FALSE)

## ======== Cluster bootstrap (admissions resampled within department): 95% CIs
B <- 10000   ## resamples for the percentile 95% CIs

boot_counts <- bind_rows(
  admitted %>% count(afdeling, patientID, alarms_per_hour) %>% mutate(category = "Average", group = "Total"),
  admitted %>% filter(!is.na(diag_group)) %>% count(afdeling, patientID, group = diag_group, alarms_per_hour) %>% mutate(category = "Reason of admission"),
  monit.hours %>% count(afdeling, patientID, group = selection, alarms_per_hour) %>% mutate(category = "Monitoring intensity"),
  ga.hours %>% mutate(afdeling = "Neonatologie", group = as.character(GestAgeGroup)) %>% count(afdeling, patientID, group, alarms_per_hour) %>% mutate(category = "Gestational age"),
  windows %>% mutate(group = as.character(hour(window_start) + minute(window_start) / 60)) %>% count(afdeling, patientID, group, alarms_per_hour) %>% mutate(category = "Time of day"))
rm(admitted, monit.hours, ga.hours, windows)

eligible_groups <- bind_rows(
  tibble(category = "Average", afdeling = unit_average$afdeling, group = "Total"),
  group_summary %>% filter(diag_group != "Total", !is.na(diag_group), cases > 25) %>% transmute(category = "Reason of admission", afdeling, group = diag_group),
  monit.stats %>% filter(alarm_active_hours > 1000) %>% transmute(category = "Monitoring intensity", afdeling, group = selection),
  GA_Analysis %>% filter(n_patients > 25) %>% transmute(category = "Gestational age", afdeling = "Neonatologie", group = as.character(GestAgeGroup)),
  rate_over_day %>% transmute(category = "Time of day", afdeling, group = as.character(time_of_day)))
boot_counts <- boot_counts %>% semi_join(eligible_groups, by = c("category", "afdeling", "group"))

observed_groups <- bind_rows(monitoring_rows, time_rows, reason_rows, ga_rows) %>% select(category, afdeling, min_group, max_group)

clusters <- boot_counts %>% distinct(afdeling, patientID) %>% mutate(cluster = row_number())
boot_counts <- boot_counts %>%
  left_join(clusters, by = c("afdeling", "patientID")) %>%
  arrange(category, afdeling, group, alarms_per_hour) %>%
  mutate(cv = consecutive_id(category, afdeling, group, alarms_per_hour))
cv_table <- boot_counts %>% distinct(cv, category, afdeling, group, alarms_per_hour) %>% arrange(cv)

set.seed(20261029)
boot_results <- NULL
for (b in 0:B) {
  if (b == 0) {
    w <- rep(1L, nrow(clusters))
  } else {
    w <- tabulate(clusters %>% group_by(afdeling) %>% slice_sample(prop = 1, replace = TRUE) %>% pull(cluster),
                  nbins = nrow(clusters))
  }
  boot_medians <- cv_table %>%
    mutate(n = rowsum(boot_counts$n * w[boot_counts$cluster], boot_counts$cv)[, 1]) %>%
    group_by(category, afdeling, group) %>%
    mutate(cum = cumsum(n), N = sum(n)) %>%
    summarise(median = (alarms_per_hour[which(cum >= floor((N[1] + 1) / 2))[1]] +
                        alarms_per_hour[which(cum >= floor(N[1] / 2) + 1)[1]]) / 2, .groups = 'drop')
  boot_results <- bind_rows(boot_results,
    boot_medians %>% filter(category != "Average") %>%
      left_join(observed_groups, by = c("category", "afdeling")) %>%
      group_by(category, afdeling) %>%
      summarise(minimum = if (category[1] == "Monitoring intensity") median[group == "regular"][1] else min(median),
                maximum = max(median),
                min_group_median = median[group == min_group[1]][1],
                max_group_median = median[group == max_group[1]][1], .groups = 'drop') %>%
      left_join(boot_medians %>% filter(category == "Average") %>% select(afdeling, average = median), by = "afdeling") %>%
      mutate(minus_pct = (average - minimum) / average * 100,
             plus_pct = (maximum - average) / average * 100,
             minus_pct_fixed = (average - min_group_median) / average * 100,
             plus_pct_fixed = (max_group_median - average) / average * 100,
             b = b))
}

boot_results %>% filter(b == 0) %>%
  mutate(department = c("Intensive Care Kinderen" = "PICU", "Neonatologie" = "NICU")[afdeling]) %>%
  left_join(values, by = c("category", "department"), suffix = c("_b0", "")) %>%
  transmute(category, department, d_min = minimum_b0 - minimum, d_max = maximum_b0 - maximum, d_avg = average_b0 - average,
            d_min_group = min_group_median - minimum, d_max_group = max_group_median - maximum)

boot_ci <- boot_results %>% filter(b > 0) %>%
  mutate(department = c("Intensive Care Kinderen" = "PICU", "Neonatologie" = "NICU")[afdeling]) %>%
  group_by(category, department) %>%
  summarise(across(c(average, minimum, maximum, min_group_median, max_group_median, minus_pct, plus_pct, minus_pct_fixed, plus_pct_fixed),
                   list(low = ~ quantile(.x, 0.025, na.rm = TRUE), high = ~ quantile(.x, 0.975, na.rm = TRUE))),
            n_missing = sum(is.na(min_group_median) | is.na(max_group_median)), .groups = 'drop')
print(boot_ci, width = Inf)

values_ci <- values %>% left_join(boot_ci, by = c("category", "department"))
write.csv(values_ci, "output//primaryOutput//figure2_summary_values_ci.csv", row.names = FALSE)

## ======== Plot
categories <- c("Monitoring intensity", "Reason of admission", "Time of day", "Gestational age")
category_display <- c("Monitoring intensity" = "Active alarm sources", "Reason of admission" = "Reason of admission",
                      "Time of day" = "Time of day", "Gestational age" = "Gestational age")

fig <- values_ci %>% mutate(
  lo_name = case_when(category == "Monitoring intensity" ~ "Standard alarm sources",
                      category == "Time of day" ~ sub("-", "–", min_label), TRUE ~ min_label),
  hi_name = case_when(category == "Monitoring intensity" ~ c("regular" = "Standard alarm sources", "regular + Vent" = "+ Ventilator", "regular + BP" = "+ ABP", "regular + Vasc" = "+ Vasc",
                                                             "regular + BP + Vent" = "+ Ventilator + ABP", "regular + Vasc + Vent" = "+ Ventilator + Vasc",
                                                             "regular + BP + Vasc" = "+ ABP + Vasc", "regular + BP + Vasc + Vent" = "+ Ventilator + ABP + Vasc")[max_label],

                      category == "Time of day" ~ sub("-", "–", max_label), TRUE ~ max_label),
  i = match(category, categories),
  y_position = 9.7 - (i - 1) * 1.1 + ifelse(department == "PICU", 0.2, -0.2),
  lo_pct = sprintf("−%d%%", minus_pct), hi_pct = sprintf("+%d%%", plus_pct),
  lo_x = pmin(minimum, min_group_median_low) - 0.8, hi_x = pmax(maximum, max_group_median_high) + 0.8)
ci <- bind_rows(fig %>% transmute(department, y_position, low = min_group_median_low, high = min_group_median_high),
                fig %>% transmute(department, y_position, low = max_group_median_low, high = max_group_median_high))
cap <- 0.07

axis_y <- 5.55
x_ticks <- seq(0, 60, 10)
unit_colours <- c("PICU" = "#2F4858", "NICU" = "#E66101")

ly <- 7.75
legend_layers <- list(
  annotate("rect", xmin = 51.6, xmax = 71.6, ymin = 6.9, ymax = 8.42, fill = "white", colour = "grey55", linewidth = 0.3),
  annotate("text", x = 52.2, y = 8.3, label = "Legend", hjust = 0, vjust = 1, size = 2.8, fontface = "bold"),
  annotate("segment", x = 55, xend = 65, y = ly, yend = ly, colour = "grey82", linewidth = 4.5),
  annotate("segment", x = c(54.2, 64.2), xend = c(55.8, 65.8), y = ly, yend = ly, linewidth = 0.45),
  annotate("segment", x = c(54.2, 55.8, 64.2, 65.8), xend = c(54.2, 55.8, 64.2, 65.8), y = ly - 0.05, yend = ly + 0.05, linewidth = 0.45),
  annotate("segment", x = 60, xend = 60, y = ly - 0.09, yend = ly + 0.09, linewidth = 1.2),
  annotate("text", x = 60, y = ly + 0.24, label = "Department median", size = 2.5, fontface = "bold"),
  annotate("text", x = 59.6, y = ly + 0.1, label = "−%", hjust = 1, vjust = 0, size = 2.3),
  annotate("text", x = 60.4, y = ly + 0.1, label = "+%", hjust = 0, vjust = 0, size = 2.3),
  annotate("text", x = 55, y = ly - 0.15, label = "Lowest group", vjust = 1, size = 2.4, fontface = "italic"),
  annotate("text", x = 65, y = ly - 0.15, label = "Highest group", vjust = 1, size = 2.4, fontface = "italic"),
  annotate("segment", x = 52.4, xend = 54, y = 7.15, yend = 7.15, linewidth = 0.45),
  annotate("segment", x = c(52.4, 54), xend = c(52.4, 54), y = 7.1, yend = 7.2, linewidth = 0.45),
  annotate("text", x = 54.6, y = 7.15, label = "95% CI of group median", hjust = 0, size = 2.4),
  annotate("segment", x = 62.6, xend = 64.2, y = 7.15, yend = 7.15, colour = unit_colours["PICU"], alpha = 0.65, linewidth = 3.5),
  annotate("text", x = 64.8, y = 7.15, label = "PICU", hjust = 0, size = 2.4),
  annotate("segment", x = 67.4, xend = 69, y = 7.15, yend = 7.15, colour = unit_colours["NICU"], alpha = 0.65, linewidth = 3.5),
  annotate("text", x = 69.6, y = 7.15, label = "NICU", hjust = 0, size = 2.4))

plot <- ggplot() +
  annotate("segment", x = x_ticks, xend = x_ticks, y = axis_y, yend = 10.2, colour = "grey92", linewidth = 0.3) +
  geom_segment(data = fig, aes(x = minimum, xend = maximum, y = y_position, yend = y_position), colour = "white", linewidth = 7) +
  geom_segment(data = fig, aes(x = minimum, xend = maximum, y = y_position, yend = y_position, colour = department), linewidth = 7, alpha = 0.65, show.legend = FALSE) +
  geom_segment(data = ci, aes(x = low, xend = high, y = y_position, yend = y_position), linewidth = 0.6) +
  geom_segment(data = ci, aes(x = low, xend = low, y = y_position - cap, yend = y_position + cap), linewidth = 0.6) +
  geom_segment(data = ci, aes(x = high, xend = high, y = y_position - cap, yend = y_position + cap), linewidth = 0.6) +
  geom_segment(data = fig, aes(x = average, xend = average, y = y_position - 0.13, yend = y_position + 0.13), linewidth = 1.6) +
  geom_text(data = fig, aes(x = lo_x, y = y_position, label = lo_name), hjust = 1, size = 3.2, fontface = "italic") +
  geom_text(data = fig, aes(x = ifelse(maximum > 50, maximum, hi_x), y = ifelse(maximum > 50, y_position - 0.17, y_position),
                            label = hi_name, hjust = ifelse(maximum > 50, 1, 0), vjust = ifelse(maximum > 50, 1, 0.5)),
            size = 3.2, fontface = "italic") +
  geom_text(data = fig, aes(x = average - 0.5, y = y_position + ifelse(department == "PICU", 0.15, -0.15), label = lo_pct,
                            vjust = ifelse(department == "PICU", 0, 1)), hjust = 1, size = 2.8) +
  geom_text(data = fig, aes(x = average + 0.5, y = y_position + ifelse(department == "PICU", 0.15, -0.15), label = hi_pct,
                            vjust = ifelse(department == "PICU", 0, 1)), hjust = 0, size = 2.8) +
  annotate("text", x = 25, y = 9.7 - 3 * 1.1 + 0.2, label = "-- Gestational age not applicable for pediatric ICU --",
           size = 3, fontface = "italic", colour = "grey45") +
  annotate("text", x = -16, y = 9.7 - (seq_along(categories) - 1) * 1.1, label = category_display[categories],
           hjust = 0, fontface = "bold.italic", size = 3.8) +
  annotate("segment", x = 0, xend = 65, y = axis_y, yend = axis_y, linewidth = 0.5) +
  annotate("segment", x = x_ticks, xend = x_ticks, y = axis_y, yend = axis_y - 0.08, linewidth = 0.5) +
  annotate("text", x = x_ticks, y = axis_y - 0.15, label = x_ticks, vjust = 1, size = 3.4) +
  annotate("text", x = 32.5, y = axis_y - 0.5, label = "Median alarm rate (alarms per alarm-active hour)", fontface = "bold", size = 3.6) +
  legend_layers +
  scale_colour_manual(values = unit_colours) +
  coord_cartesian(xlim = c(-16, 72), ylim = c(axis_y - 0.75, 10.2), expand = FALSE, clip = "off") +
  theme_void(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = unit(c(6, 6, 6, 6), "mm"))
ggsave("output//primaryOutput//figure2_summary.png", plot, width = 36, height = 15, units = "cm", dpi = 300)
