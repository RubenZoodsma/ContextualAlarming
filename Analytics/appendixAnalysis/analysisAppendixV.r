## Appendix V: alarm burden (A) and alarm composition (B) over the course of the day, PICU and NICU

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

study_days <- seq(as.Date(STUDY_START, tz = "Europe/Amsterdam"), as.Date(STUDY_END, tz = "Europe/Amsterdam") - 1, by = "day")
data_tod <- data %>%
  transmute(afdeling, bed_naam, patientID, Alarm.group,
            day = as_date(alarm_start),
            time_of_day = hour(alarm_start) + floor(minute(alarm_start) / 15) * 0.25)

## ======== A: alarm burden over the day, whole department and per unit
bed_map <- read.csv("..//..//Data//bed_identifiers.csv", sep = ";", fileEncoding = "UTF-8-BOM") %>%
  filter(Unit %in% c("unit_1", "unit_2", "unit_3", "unit_4")) %>%
  select(bed_naam = BedName, Unit, bay_department = Department)
data_tod <- data_tod %>% left_join(bed_map, by = "bed_naam")

unit_letters <- data_tod %>% filter(!is.na(Unit)) %>% count(bay_department, Unit, name = "alarms") %>%
  arrange(bay_department, Unit) %>% group_by(bay_department) %>% mutate(bay = paste("Unit", LETTERS[row_number()])) %>% ungroup()
print(unit_letters)
data_tod <- data_tod %>% left_join(unit_letters %>% select(bay_department, Unit, bay), by = c("bay_department", "Unit"))
data_tod %>% count(afdeling, bay_department, name = "alarms") %>% print()

ward_counts <- bind_rows(
  data_tod %>% count(afdeling, day, time_of_day, name = "alarms") %>% mutate(Unit = "Overall"),
  data_tod %>% filter(!is.na(bay)) %>% count(afdeling = bay_department, Unit = bay, day, time_of_day, name = "alarms"))
ward_bins <- ward_counts %>%
  complete(nesting(afdeling, Unit), day = study_days, time_of_day = seq(0, 23.75, 0.25), fill = list(alarms = 0)) %>%
  semi_join(ward_counts %>% distinct(afdeling, Unit, day), by = c("afdeling", "Unit", "day")) %>%   ## days in use only
  mutate(alarms_per_hour = alarms * 4)
ward_bins %>% distinct(afdeling, Unit, day) %>% count(afdeling, Unit, name = "days_in_use") %>% mutate(study_days = length(study_days))
rm(ward_counts)

ward_over_day <- ward_bins %>%
  group_by(afdeling, Unit, time_of_day) %>%
  summarise(
    days = n(),
    mean_alarms_per_hour = mean(alarms_per_hour),
    median_alarms_per_hour = median(alarms_per_hour),
    q25_alarms_per_hour = quantile(alarms_per_hour, 0.25),
    q75_alarms_per_hour = quantile(alarms_per_hour, 0.75),
    .groups = 'drop') %>%
  mutate(bin_mid = time_of_day + 0.125)
ward_reference <- ward_bins %>% group_by(afdeling, Unit) %>% summarise(mean_alarms_per_hour = round(mean(alarms_per_hour)), .groups = 'drop')
print(ward_reference)
ward_over_day %>% group_by(afdeling, Unit) %>%
  summarise(lowest = round(min(mean_alarms_per_hour)), lowest_bin_start = time_of_day[which.min(mean_alarms_per_hour)],
            highest = round(max(mean_alarms_per_hour)), highest_bin_start = time_of_day[which.max(mean_alarms_per_hour)])
write.csv(ward_over_day, "output//appendices//appendixV_A_ward_over_day.csv", row.names = FALSE)
rm(ward_bins)

## ======== B: alarm composition over the day
composition_reference <- data_tod %>% group_by(afdeling) %>%
  summarise(phys_reference = round(sum(Alarm.group == "Physiological") / n() * 100, 1),
            tech_reference = round(sum(Alarm.group == "Technical") / n() * 100, 1))
print(composition_reference)

composition_over_day <- data_tod %>%
  group_by(afdeling, day, time_of_day) %>%
  summarise(alarms = n(), phys = sum(Alarm.group == "Physiological"), .groups = 'drop') %>%
  mutate(phys_share = phys / alarms * 100) %>%
  group_by(afdeling, time_of_day) %>%
  summarise(
    alarms = sum(alarms),
    phys_percent = sum(phys) / sum(alarms) * 100,
    q25_daily_phys_percent = quantile(phys_share, 0.25),
    q75_daily_phys_percent = quantile(phys_share, 0.75),
    .groups = 'drop') %>%
  left_join(composition_reference, by = "afdeling") %>%
  mutate(phys_difference = (phys_percent - phys_reference) / phys_reference * 100,
         q25_difference = (q25_daily_phys_percent - phys_reference) / phys_reference * 100,
         q75_difference = (q75_daily_phys_percent - phys_reference) / phys_reference * 100,
         bin_mid = time_of_day + 0.125)
composition_over_day %>% group_by(afdeling) %>%
  summarise(lowest_phys = round(min(phys_percent), 1), lowest_bin_start = time_of_day[which.min(phys_percent)],
            highest_phys = round(max(phys_percent), 1), highest_bin_start = time_of_day[which.max(phys_percent)])
write.csv(composition_over_day, "output//appendices//appendixV_B_composition_over_day.csv", row.names = FALSE)

## ======== Alarm rate per alarm-active hour over the day
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
  summarise(
    windows = n(),
    mean_alarms_per_hour = round(mean(alarms_per_hour), 1),
    median_alarms_per_hour = median(alarms_per_hour),
    q25_alarms_per_hour = quantile(alarms_per_hour, 0.25),
    q75_alarms_per_hour = quantile(alarms_per_hour, 0.75),
    .groups = 'drop') %>%
  mutate(window_mid = (time_of_day + 0.5) %% 24,
         window_label = paste0(sprintf("%02d:%02d", floor(time_of_day), round((time_of_day %% 1) * 60)), "-",
                               sprintf("%02d:%02d", floor((time_of_day + 1) %% 24), round((time_of_day %% 1) * 60))))
windows %>% group_by(afdeling) %>% summarise(windows = n(), median_all_windows = median(alarms_per_hour))
rm(windows)
write.csv(rate_over_day, "output//appendices//appendixV_rate_per_patient_over_day.csv", row.names = FALSE)

admitted <- pat_char %>% ungroup() %>%
  mutate(start = floor_date(pmax(admission_start, STUDY_START), "1 hour"), end = pmin(admission_end, STUDY_END)) %>% filter(end > start) %>%
  reframe(alarm_floor_hr = seq(start, end - 1, by = "1 hour"), .by = c(patientID, afdeling)) %>%
  left_join(data %>% mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>% count(patientID, afdeling, alarm_floor_hr, name = "alarms_per_hour"),
            by = c("patientID", "afdeling", "alarm_floor_hr")) %>%
  mutate(alarms_per_hour = coalesce(alarms_per_hour, 0L))
unit_average <- admitted %>% group_by(afdeling) %>% summarise(average = median(alarms_per_hour[alarms_per_hour > 0]))
rm(admitted)

extremes <- rate_over_day %>%
  left_join(unit_average, by = "afdeling") %>%
  group_by(afdeling, average) %>%
  summarise(
    lowest = min(median_alarms_per_hour),
    lowest_windows = paste(window_label[median_alarms_per_hour == min(median_alarms_per_hour)], collapse = ", "),
    highest = max(median_alarms_per_hour),
    highest_windows = paste(window_label[median_alarms_per_hour == max(median_alarms_per_hour)], collapse = ", "),
    .groups = 'drop') %>%
  mutate(lowest_pct = round((lowest - average) / average * 100),
         highest_pct = round((highest - average) / average * 100),
         highest_ratio = round(highest / average, 2),
         highest_vs_lowest_ratio = round(highest / lowest, 2))
print(extremes, width = Inf)
write.csv(extremes, "output//appendices//appendixV_extremes.csv", row.names = FALSE)

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
unit_panels <- c("Intensive Care Kinderen" = "Pediatric ICU", "Neonatologie" = "Neonatal ICU")
unit_short <- c("Intensive Care Kinderen" = "PICU", "Neonatologie" = "NICU")
shift_labels <- data.frame(x = c(3.5, 11, 19),
                           label = c("Night", "Day", "Evening"))

unit_columns <- c("Overall", sort(unique(unit_letters$bay)))
unit_headers <- setNames(c("Overall: all units combined", sort(unique(unit_letters$bay))), unit_columns)
empty_panels <- expand.grid(department = unit_panels, column = unit_columns, stringsAsFactors = FALSE) %>%
  anti_join(ward_over_day %>% distinct(department = unit_panels[afdeling], column = Unit), by = c("department", "column")) %>%
  left_join(unit_letters %>% count(department = unit_panels[bay_department], name = "units"), by = "department") %>%
  mutate(department = factor(department, levels = unit_panels), column = factor(column, levels = unit_columns),
         note = paste0("Not applicable:\nthe ", department, " has ", units, " units"))
plot <- ward_over_day %>%
  mutate(department = factor(unit_panels[afdeling], levels = unit_panels), column = factor(Unit, levels = unit_columns)) %>%
  ggplot(aes(x = bin_mid)) +
  geom_ribbon(aes(ymin = q25_alarms_per_hour, ymax = q75_alarms_per_hour), fill = "#386cb0", alpha = 0.12) +
  geom_line(aes(y = q25_alarms_per_hour), color = "#386cb0", linewidth = 0.5) +
  geom_line(aes(y = q75_alarms_per_hour), color = "#386cb0", linewidth = 0.5) +
  geom_line(aes(y = mean_alarms_per_hour), color = "grey20", linewidth = 0.7) +
  geom_vline(xintercept = c(7, 15, 23), linetype = "dashed", color = "grey50") +
  geom_hline(data = ward_reference %>% mutate(department = factor(unit_panels[afdeling], levels = unit_panels), column = factor(Unit, levels = unit_columns)),
             aes(yintercept = mean_alarms_per_hour), linetype = "dashed", color = "black") +
  geom_label(data = ward_reference %>% mutate(department = factor(unit_panels[afdeling], levels = unit_panels), column = factor(Unit, levels = unit_columns)),
             aes(x = 0.3, y = mean_alarms_per_hour, label = paste0("Mean = ", mean_alarms_per_hour)),
             hjust = 0, vjust = -0.3, size = 2.6, linewidth = 0, fill = "white", label.r = unit(0, "pt")) +
  geom_text(data = shift_labels, aes(x = x, y = Inf, label = label), vjust = 1.5, fontface = "bold.italic", size = 2.8) +
  geom_blank(data = data.frame(column = factor(rep(unit_columns, each = 2), levels = unit_columns),
                               y = rep(c(0, 1000), length(unit_columns)) * ifelse(rep(unit_columns, each = 2) == "Overall", 1, 0.45)),
             aes(x = 12, y = y)) +
  geom_segment(data = data.frame(column = factor("Overall", levels = unit_columns)),
               x = -Inf, xend = Inf, y = -Inf, yend = -Inf, linewidth = 1.6, color = "grey30", inherit.aes = FALSE) +
  geom_rect(data = empty_panels, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey92", inherit.aes = FALSE) +
  geom_label(data = empty_panels, aes(x = 12, y = 225, label = note), size = 3, color = "grey45", fill = "grey97",
             linewidth = 0.3, label.r = unit(2, "pt"), lineheight = 0.95, inherit.aes = FALSE) +
  facet_grid(column ~ department, scales = "free_y", labeller = labeller(column = as_labeller(unit_headers))) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  scale_y_continuous(breaks = function(limits) if (limits[2] > 600) seq(0, 1000, 200) else seq(0, 450, 50),
                     expand = expansion(mult = 0.02)) +
  theme_Publication(base_size = 11) +
  theme(panel.spacing.x = unit(0.6, "cm"),
        panel.spacing.y = unit(c(1.8, rep(0.4, max(0, length(unit_columns) - 2))), "cm"),
        plot.title = element_text(hjust = 0)) +
  xlab("Time of day") +
  ylab("Alarms per hour")
plot <- ggplotGrob(plot)
overall_strip <- which(plot$layout$name == "strip-r-1")
plot$grobs[[overall_strip]]$grobs[[1]]$children[[1]]$gp$fill <- "#2c2c2c"
plot$grobs[[overall_strip]]$grobs[[1]]$children[[1]]$gp$col <- "#2c2c2c"
plot$grobs[[overall_strip]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "white"
ggsave("output//appendices//appendixV_A_burden_over_day.png", plot, width = 30, height = 30, units = "cm", dpi = 300)

plot <- composition_over_day %>%
  mutate(panel = factor(unit_panels[afdeling], levels = unit_panels)) %>%
  ggplot(aes(x = bin_mid)) +
  geom_ribbon(aes(ymin = q25_difference, ymax = q75_difference), fill = "#386cb0", alpha = 0.12) +
  geom_line(aes(y = q25_difference), color = "#386cb0", linewidth = 0.6) +
  geom_line(aes(y = q75_difference), color = "#386cb0", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.8) +
  geom_line(aes(y = phys_difference), color = "grey20", linewidth = 0.8) +
  geom_vline(xintercept = c(7, 15, 23), linetype = "dashed", color = "grey50") +
  geom_label(data = composition_reference %>% mutate(panel = factor(unit_panels[afdeling], levels = unit_panels)),
             aes(x = 11, y = Inf, label = paste0(unit_short[afdeling], ": ", phys_reference, "% of all alarms physiological")),
             vjust = 3.4, size = 3, linewidth = 0, fill = "white", label.r = unit(0, "pt")) +
  geom_label(data = composition_reference %>% mutate(panel = factor(unit_panels[afdeling], levels = unit_panels)),
             aes(x = 11, y = -Inf, label = paste0(unit_short[afdeling], ": ", tech_reference, "% of all alarms technical")),
             vjust = -1.6, size = 3, linewidth = 0, fill = "white", label.r = unit(0, "pt")) +
  geom_text(data = shift_labels, aes(x = x, y = Inf, label = label), vjust = 1.5, fontface = "bold.italic", size = 3.5) +
  geom_segment(data = data.frame(panel = factor(unit_panels[["Intensive Care Kinderen"]], levels = unit_panels), y = c(4, -4), yend = c(20, -20)),
               aes(x = 26.3, xend = 26.3, y = y, yend = yend), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.6, color = "grey40") +
  geom_text(data = data.frame(panel = factor(unit_panels[["Intensive Care Kinderen"]], levels = unit_panels), y = c(24, -24), hjust = c(0, 1),
                              label = c("Relatively more\nphysiological alarms", "Relatively more\ntechnical alarms")),
            aes(x = 26.3, y = y, label = label, hjust = hjust), angle = 90, vjust = 0.5, size = 2.6, fontface = "bold", lineheight = 0.9) +
  facet_wrap(~ panel, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 27)) +
  scale_y_continuous(breaks = seq(-60, 80, 10)) +
  coord_cartesian(xlim = c(0, 24), ylim = c(-60, 80), clip = "off") +
  theme_Publication(base_size = 12) +
  theme(panel.spacing = unit(2.2, "cm"), plot.title = element_text(hjust = 0)) +
  xlab("Time of day") +
  ylab("Relative difference (%) from the unit's\nshare of physiological alarms")
ggsave("output//appendices//appendixV_B_composition_over_day.png", plot, width = 34, height = 14, units = "cm", dpi = 300)
