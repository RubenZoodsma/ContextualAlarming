## Table 2: alarm count and duration by device and priority, PICU and NICU

library(dplyr)
library(tidyr)
library(readr)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

data %>% group_by(afdeling) %>% summarise(alarms = n(),
  missing_device = sum(is.na(Device)),
  missing_type = sum(is.na(Alarm.group)),
  missing_duration = sum(is.na(duur_sec)),
  first_alarm = min(alarm_start), last_alarm = max(alarm_start))

## ======== Alarm count, n (%)
device_count <- data %>% group_by(afdeling, Device, alarm_prio_bron, Alarm.group) %>%
  summarise(count = n(), .groups = 'drop')
device_count <- bind_rows(device_count,
                          device_count %>% group_by(afdeling, Device, Alarm.group) %>% summarise(count = sum(count), .groups = 'drop') %>% mutate(alarm_prio_bron = "Total")) %>%
  mutate(alarm_prio_bron = factor(alarm_prio_bron, levels = c("Laag", "Medium", "Hoog", "Total")))
device_count %>% filter(alarm_prio_bron == "Total") %>% group_by(afdeling) %>% summarise(count = sum(count))

device_rows <- device_count %>% group_by(afdeling, Device, alarm_prio_bron) %>% summarise(n = sum(count), .groups = 'drop') %>%
  group_by(afdeling, Device) %>% mutate(pct = round(n / n[alarm_prio_bron == "Total"] * 100, 1)) %>%
  group_by(afdeling) %>% mutate(pct = ifelse(alarm_prio_bron == "Total", round(n / sum(n[alarm_prio_bron == "Total"]) * 100, 1), pct)) %>%
  ungroup() %>% arrange(afdeling, Device, alarm_prio_bron)
print(device_rows, n = Inf, width = Inf)

type_rows <- device_count %>% group_by(afdeling, Device, alarm_prio_bron) %>%
  mutate(pct = round(count / sum(count) * 100, 1)) %>%
  ungroup() %>% arrange(afdeling, Device, alarm_prio_bron, Alarm.group)
print(type_rows, n = Inf, width = Inf)

total_rows <- device_count %>% group_by(afdeling, alarm_prio_bron) %>% summarise(n = sum(count), .groups = 'drop') %>%
  group_by(afdeling) %>% mutate(pct = round(n / n[alarm_prio_bron == "Total"] * 100, 1)) %>% ungroup()
print(total_rows, n = Inf, width = Inf)
total_type_rows <- device_count %>% group_by(afdeling, alarm_prio_bron, Alarm.group) %>% summarise(n = sum(count), .groups = 'drop') %>%
  group_by(afdeling, alarm_prio_bron) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>% ungroup()
print(total_type_rows, n = Inf, width = Inf)

alarm_count <- bind_rows(device_rows %>% mutate(Alarm.group = "All"), type_rows %>% rename(n = count),
                    total_rows %>% mutate(Device = "Total", Alarm.group = "All"), total_type_rows %>% mutate(Device = "Total")) %>%
            select(afdeling, Device, Alarm.group, alarm_prio_bron, n, pct)
write.csv(alarm_count, "output//primaryOutput//table2_alarm_count.csv", row.names = FALSE)

## ======== Alarm duration, s, median [IQR]
duration <- bind_rows(
  data %>% group_by(afdeling, alarm_prio_bron, Alarm.group) %>% summarise(alarms = n(), median_duration = median(duur_sec, na.rm=TRUE), q25_duration = quantile(duur_sec, 0.25, na.rm=TRUE), q75_duration = quantile(duur_sec, 0.75, na.rm=TRUE), .groups = 'drop'),
  data %>% group_by(afdeling, alarm_prio_bron) %>% summarise(alarms = n(), median_duration = median(duur_sec, na.rm=TRUE), q25_duration = quantile(duur_sec, 0.25, na.rm=TRUE), q75_duration = quantile(duur_sec, 0.75, na.rm=TRUE), .groups = 'drop') %>% mutate(Alarm.group = "Average"),
  data %>% group_by(afdeling, Alarm.group) %>% summarise(alarms = n(), median_duration = median(duur_sec, na.rm=TRUE), q25_duration = quantile(duur_sec, 0.25, na.rm=TRUE), q75_duration = quantile(duur_sec, 0.75, na.rm=TRUE), .groups = 'drop') %>% mutate(alarm_prio_bron = "Total"),
  data %>% group_by(afdeling) %>% summarise(alarms = n(), median_duration = median(duur_sec, na.rm=TRUE), q25_duration = quantile(duur_sec, 0.25, na.rm=TRUE), q75_duration = quantile(duur_sec, 0.75, na.rm=TRUE), .groups = 'drop') %>% mutate(alarm_prio_bron = "Total", Alarm.group = "Average")) %>%
  mutate(alarm_prio_bron = factor(alarm_prio_bron, levels = c("Laag", "Medium", "Hoog", "Total"))) %>%
  arrange(afdeling, Alarm.group == "Average", Alarm.group, alarm_prio_bron)
print(duration, n = Inf, width = Inf)
write.csv(duration, "output//primaryOutput//table2_duration.csv", row.names = FALSE)

## ======== Table 2 layout
alarm_count %>% filter(Alarm.group == "unknown", alarm_prio_bron == "Total") %>% select(afdeling, Device, n)

count_layout <- tibble(
  Device      = c("Monitor", "Monitor", "Monitor", NA, "Ventilator", "Ventilator", "Ventilator", "Ventilator", NA,
                  "Incubator", "Incubator", "Incubator", NA, "Thermoregulator", "Thermoregulator", "Thermoregulator", NA, "Total", "Total", "Total", "Total"),
  Alarm.group = c("All", "Physiological", "Technical", NA, "All", "Physiological", "Technical", "unknown", NA,
                  "All", "Physiological", "Technical", NA, "All", "Physiological", "Technical", NA, "All", "Physiological", "Technical", "unknown"),
  label       = c("Physiological Monitor", "Physiological", "Technical", "", "Mech. Ventilator", "Physiological", "Technical", "Unknown", "",
                  "Incubator", "Physiological", "Technical", "", "Thermoregulator", "Physiological", "Technical", "", "Total", "Physiological", "Technical", "Unknown")) %>%
  mutate(row = row_number())

count_rows <- count_layout %>%
  cross_join(tibble(afdeling = c("Intensive Care Kinderen", "Neonatologie"))) %>%
  cross_join(tibble(alarm_prio_bron = c("Laag", "Medium", "Hoog", "Total"))) %>%
  left_join(alarm_count %>% mutate(alarm_prio_bron = as.character(alarm_prio_bron)), by = c("afdeling", "Device", "Alarm.group", "alarm_prio_bron")) %>%
  mutate(unit = if_else(afdeling == "Intensive Care Kinderen", "picu", "nicu"),
         k = match(alarm_prio_bron, c("Laag", "Medium", "Hoog", "Total")),
         n_txt = case_when(is.na(Device) ~ "", is.na(n) | n == 0 ~ "-", TRUE ~ as.character(n)),
         pct_txt = case_when(is.na(Device) | is.na(n) | n == 0 ~ "", TRUE ~ sprintf("%.1f", pct)))
count_rows <- bind_rows(count_rows %>% transmute(row, label, col = paste0(unit, "_c", 2 * k - 1), value = n_txt),
                        count_rows %>% transmute(row, label, col = paste0(unit, "_c", 2 * k), value = pct_txt)) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  arrange(row) %>% mutate(section = "Alarm count")

duration_rows <- tibble(Alarm.group = c("Physiological", "Technical", "Average")) %>% mutate(label = Alarm.group, row = row_number()) %>%
  cross_join(tibble(afdeling = c("Intensive Care Kinderen", "Neonatologie"))) %>%
  cross_join(tibble(alarm_prio_bron = c("Laag", "Medium", "Hoog", "Total"))) %>%
  left_join(duration %>% mutate(alarm_prio_bron = as.character(alarm_prio_bron)), by = c("afdeling", "Alarm.group", "alarm_prio_bron")) %>%
  transmute(row, label,
            col = paste0(if_else(afdeling == "Intensive Care Kinderen", "picu", "nicu"), "_c", 2 * match(alarm_prio_bron, c("Laag", "Medium", "Hoog", "Total")) - 1),
            value = if_else(is.na(median_duration), "-", sprintf("%.0f [%.0f - %.0f]", median_duration, q25_duration, q75_duration))) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  arrange(row) %>% mutate(section = "Alarm duration, s")

table2 <- bind_rows(count_rows, tibble(section = "Alarm duration, s", label = c("", "", "")), duration_rows) %>%
  select(section, label, all_of(paste0("picu_c", 1:8)), all_of(paste0("nicu_c", 1:8)))
print(table2, n = Inf, width = Inf)

write_excel_csv(table2, "output//submissionTables//table2.csv", na = "")
