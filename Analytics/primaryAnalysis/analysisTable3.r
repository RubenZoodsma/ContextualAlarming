## Table 3: alarm rate by active alarm sources, by gestational age (NICU) and per nurse per shift, PICU and NICU

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(hms)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

## ======== Alarm rate by active alarm sources (3-hour backward window)

hours_active <- data %>% filter(!is.na(patientID)) %>%
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
              summarise(BP = any(BP_alarm),
                        Vasc = any(Vasc_alarm),
                        Vent = any(Vent_alarm),
                            .groups = 'drop'),
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
    ungroup() %>% select(afdeling, patientID, alarm_floor_hr, selection, Alarm.group, duur_sec)
rm(hours_active, component_e, target)

monit.stats <- bind_rows(monit.data, monit.data %>% mutate(selection = "Average")) %>%
  group_by(afdeling, selection, patientID, alarm_floor_hr) %>%
  summarise(
    alarms_per_hour = n(),
    phys_alarms = sum(Alarm.group == "Physiological", na.rm = TRUE),
    tech_alarms = sum(Alarm.group == "Technical", na.rm = TRUE),
    .groups = 'drop') %>%
  group_by(afdeling, selection) %>%
  summarise(
    alarm_active_hours = n(),
    mean_alarms_per_hour = round(mean(alarms_per_hour, na.rm = TRUE), 1),
    median_alarms_per_hour = median(alarms_per_hour, na.rm = TRUE),
    q25_alarms_per_hour = quantile(alarms_per_hour, 0.25, na.rm = TRUE),
    q75_alarms_per_hour = quantile(alarms_per_hour, 0.75, na.rm = TRUE),
    total_alarms = sum(alarms_per_hour, na.rm = TRUE),
    phys_percentage = round(sum(phys_alarms, na.rm = TRUE) / total_alarms * 100, 0),
    tech_percentage = round(sum(tech_alarms, na.rm = TRUE) / total_alarms * 100, 0),
    .groups = 'drop') %>%
  left_join(bind_rows(monit.data, monit.data %>% mutate(selection = "Average")) %>%
    group_by(afdeling, selection) %>%
    summarise(
      median_duration = median(duur_sec, na.rm = TRUE),
      q25_duration = quantile(duur_sec, 0.25, na.rm = TRUE),
      q75_duration = quantile(duur_sec, 0.75, na.rm = TRUE),
      .groups = 'drop'),
    by = c("afdeling", "selection")) %>%
  arrange(afdeling, factor(selection, levels = c("regular", "regular + Vent", "regular + BP", "regular + Vasc", "regular + BP + Vent",
                                                 "regular + Vasc + Vent", "regular + BP + Vasc", "regular + BP + Vasc + Vent", "Average")))
print(monit.stats, n = Inf, width = Inf)

monit.stats %>% group_by(afdeling) %>% summarise(levels_add_up = sum(alarm_active_hours[selection != "Average"]) == alarm_active_hours[selection == "Average"])
write.csv(monit.stats, "output//primaryOutput//table3_monitoring.csv", row.names = FALSE)
rm(monit.data)

## ======== Alarm rate by gestational age (NICU only)

alarm_data_by_ga <- data %>% inner_join(pat_char %>% filter(afdeling == "Neonatologie" & !is.na(GestAgeGroup)) %>% select(patientID, afdeling, GestAgeGroup), by = c("patientID", "afdeling"))

pat_char %>% ungroup() %>% filter(afdeling == "Neonatologie") %>% summarise(admissions = n(), missing_ga = sum(is.na(GestAgeGroup)))
GA_Analysis <- bind_rows(alarm_data_by_ga %>% mutate(GestAgeGroup = as.character(GestAgeGroup)), alarm_data_by_ga %>% mutate(GestAgeGroup = "Average")) %>%
  mutate(alarm_start_floor_hr = floor_date(alarm_start, "1 hour")) %>%
  group_by(patientID, alarm_start_floor_hr, GestAgeGroup) %>%
  summarise(alarms_per_hour = n(), .groups = 'drop') %>%
  group_by(GestAgeGroup) %>%
  summarise(
    n_patients = n_distinct(patientID),
    mean_aph = round(mean(alarms_per_hour, na.rm = TRUE), 1),
    median_aph = median(alarms_per_hour, na.rm = TRUE),
    q25_aph = quantile(alarms_per_hour, 0.25, na.rm = TRUE),
    q75_aph = quantile(alarms_per_hour, 0.75, na.rm = TRUE),
    .groups = 'drop') %>%
  left_join(
    bind_rows(alarm_data_by_ga %>% mutate(GestAgeGroup = as.character(GestAgeGroup)), alarm_data_by_ga %>% mutate(GestAgeGroup = "Average")) %>%
      group_by(GestAgeGroup) %>%
      summarise(
        phys_percent = round(sum(Alarm.group == "Physiological", na.rm = TRUE) / n() * 100, 0),
        tech_percent = round(sum(Alarm.group == "Technical", na.rm = TRUE) / n() * 100, 0),
        median_duration = median(duur_sec, na.rm = TRUE),
        q25_duration = quantile(duur_sec, 0.25, na.rm = TRUE),
        q75_duration = quantile(duur_sec, 0.75, na.rm = TRUE),
        .groups = 'drop'),
    by = "GestAgeGroup") %>%
  arrange(GestAgeGroup == "Average", factor(GestAgeGroup, levels = levels(pat_char$GestAgeGroup)))
print(GA_Analysis, width = Inf)
write.csv(GA_Analysis, "output//primaryOutput//table3_ga.csv", row.names = FALSE)
rm(alarm_data_by_ga)

## ======== Alarm burden to the individual nurse
shift_days <- seq(as.Date(STUDY_START, tz = "Europe/Amsterdam"), as.Date(STUDY_END, tz = "Europe/Amsterdam"), by = "day")

shifts <- bind_rows(
  tibble(day = shift_days, shift = "Day",     shift_start = as.POSIXct(paste(shift_days, "07:00:00"), tz = "Europe/Amsterdam"),     shift_end = as.POSIXct(paste(shift_days, "15:00:00"), tz = "Europe/Amsterdam")),
  tibble(day = shift_days, shift = "Evening", shift_start = as.POSIXct(paste(shift_days, "15:00:00"), tz = "Europe/Amsterdam"),     shift_end = as.POSIXct(paste(shift_days, "23:00:00"), tz = "Europe/Amsterdam")),
  tibble(day = shift_days, shift = "Night",   shift_start = as.POSIXct(paste(shift_days - 1, "23:00:00"), tz = "Europe/Amsterdam"), shift_end = as.POSIXct(paste(shift_days, "07:00:00"), tz = "Europe/Amsterdam"))) %>%
  mutate(shift_start = pmax(shift_start, STUDY_START), shift_end = pmin(shift_end, STUDY_END)) %>% filter(shift_end > shift_start) %>%
  mutate(shift_hours = as.numeric(difftime(shift_end, shift_start, units = "hours")))

census <- bed_intervals %>% semi_join(pat_char, by = c("patientID", "afdeling")) %>%
  inner_join(shifts, by = join_by(overlaps(bed_start, bed_end, shift_start, shift_end)), relationship = "many-to-many") %>%
  mutate(present_h = as.numeric(difftime(pmin(bed_end, shift_end), pmax(bed_start, shift_start), units = "hours"))) %>%
  filter(present_h > 0) %>%
  group_by(afdeling, day, shift, shift_hours) %>% summarise(census = n_distinct(bedID), .groups = 'drop')

nurses <- census %>%
  mutate(providers = case_when(
           shift == "Day" ~ ceiling(census / 1.5),
           shift == "Evening" ~ ceiling(census / 1.75),
           shift == "Night" ~ ceiling(census / 2.0)))

alarms_shift <- data %>% filter(!is.na(patientID)) %>%
  inner_join(shifts, by = join_by(between(alarm_start, shift_start, shift_end, bounds = "[)"))) %>%
  semi_join(nurses, by = c("afdeling", "day", "shift")) %>%
  select(afdeling, day, shift, Alarm.group, duur_sec)

shift_load <- nurses %>%
  left_join(alarms_shift %>% count(afdeling, day, shift, name = "total_alarms"), by = c("afdeling", "day", "shift")) %>%
  mutate(total_alarms = coalesce(total_alarms, 0L),
         alarms_per_provider_per_hour = total_alarms / (providers * shift_hours))

HCP_load <- bind_rows(shift_load, shift_load %>% mutate(shift = "Average")) %>%
  group_by(afdeling, shift) %>%
  summarise(
    shifts = n(),
    mean_alarms_per_provider_per_hour = round(mean(alarms_per_provider_per_hour, na.rm = TRUE), 1),
    median_alarms_per_provider_per_hour = round(median(alarms_per_provider_per_hour, na.rm = TRUE), 1),
    q25_alarms_per_provider_per_hour = round(quantile(alarms_per_provider_per_hour, 0.25, na.rm = TRUE), 1),
    q75_alarms_per_provider_per_hour = round(quantile(alarms_per_provider_per_hour, 0.75, na.rm = TRUE), 1),
    .groups = 'drop') %>%
  left_join(bind_rows(alarms_shift, alarms_shift %>% mutate(shift = "Average")) %>%
    group_by(afdeling, shift) %>%
    summarise(
      phys_percent = round(sum(Alarm.group == "Physiological", na.rm = TRUE) / n() * 100, 0),
      tech_percent = round(sum(Alarm.group == "Technical", na.rm = TRUE) / n() * 100, 0),
      median_duration = median(duur_sec, na.rm = TRUE),
      q25_duration = quantile(duur_sec, 0.25, na.rm = TRUE),
      q75_duration = quantile(duur_sec, 0.75, na.rm = TRUE),
      .groups = 'drop'),
    by = c("afdeling", "shift")) %>%
  arrange(afdeling, factor(shift, levels = c("Day", "Evening", "Night", "Average")))
print(HCP_load, width = Inf)
write.csv(HCP_load, "output//primaryOutput//table3_nurse.csv", row.names = FALSE)

## ======== Table 3 layout

source_layout <- tibble(
  selection = c("regular", NA, "regular + Vent", "regular + BP", "regular + Vasc", NA,
                "regular + BP + Vent", "regular + Vasc + Vent", "regular + BP + Vasc", NA, "regular + BP + Vasc + Vent", "Average"),
  label     = c("Standard monitoring only", "", "+ Ventilator", "+ ABP", "+ Vasc", "",
                "+ Ventilator + ABP", "+ Ventilator + Vasc", "+ ABP + Vasc", "", "+ Ventilator + ABP + Vasc", "Total")) %>%
  mutate(row = row_number())
source_rows <- source_layout %>%
  cross_join(tibble(afdeling = c("Intensive Care Kinderen", "Neonatologie"))) %>%
  left_join(monit.stats, by = c("afdeling", "selection")) %>%
  mutate(unit = if_else(afdeling == "Intensive Care Kinderen", "picu", "nicu"),
         c1 = paste0(alarm_active_hours, if_else(alarm_active_hours <= 1000, " †", "")),
         c2 = sprintf("%.1f", mean_alarms_per_hour),
         c3 = sprintf("%.0f", median_alarms_per_hour),
         c4 = sprintf("[%.0f - %.0f]", q25_alarms_per_hour, q75_alarms_per_hour),
         c5 = sprintf("%.0f", phys_percentage),
         c6 = sprintf("%.0f", tech_percentage),
         c7 = sprintf("%.0f [%.0f - %.0f]", median_duration, q25_duration, q75_duration)) %>%
  mutate(across(c1:c7, ~ case_when(is.na(selection) ~ "", is.na(alarm_active_hours) ~ "-", TRUE ~ .x))) %>%
  select(row, label, unit, c1:c7) %>%
  pivot_wider(names_from = unit, values_from = c1:c7, names_glue = "{unit}_{.value}") %>%
  arrange(row) %>% mutate(section = "Active alarm sources", label2 = "")

ga_rows <- tibble(
  GestAgeGroup = c("Extremely preterm", "Very preterm", "Moderate to late preterm", "Term", "Average"),
  label        = c("Extremely preterm", "Very preterm", "Moderate to late preterm", "Term", "Total"),
  label2       = c("< 28+0 weeks", "28+0 - 31+6 weeks", "32+0 - 36+6 weeks", "≥ 37+0 weeks", "")) %>%
  left_join(GA_Analysis, by = "GestAgeGroup") %>%
  mutate(nicu_c1 = paste0(n_patients, if_else(GestAgeGroup != "Average" & n_patients <= 25, " †", "")),
         nicu_c2 = sprintf("%.1f", mean_aph),
         nicu_c3 = sprintf("%.0f", median_aph),
         nicu_c4 = sprintf("[%.0f - %.0f]", q25_aph, q75_aph),
         nicu_c5 = sprintf("%.0f", phys_percent),
         nicu_c6 = sprintf("%.0f", tech_percent),
         nicu_c7 = sprintf("%.0f [%.0f - %.0f]", median_duration, q25_duration, q75_duration)) %>%
  mutate(across(nicu_c1:nicu_c7, ~ if_else(is.na(n_patients), "-", .x))) %>%
  select(label, label2, nicu_c1:nicu_c7) %>% mutate(section = "Gestational age")

shift_rows <- tibble(
  shift  = c("Day", "Evening", "Night", "Average"),
  label  = c("Day, 07h - 15h", "Evening, 15h - 23h", "Night, 23h - 07h", "All shifts"),
  label2 = c("1 : 1.5", "1 : 1.75", "1 : 2", "")) %>%
  mutate(row = row_number()) %>%
  cross_join(tibble(afdeling = c("Intensive Care Kinderen", "Neonatologie"))) %>%
  left_join(HCP_load, by = c("afdeling", "shift")) %>%
  mutate(unit = if_else(afdeling == "Intensive Care Kinderen", "picu", "nicu"),
         c1 = as.character(shifts),
         c2 = sprintf("%.1f", mean_alarms_per_provider_per_hour),
         c3 = sprintf("%.1f", median_alarms_per_provider_per_hour),
         c4 = sprintf("[%.1f - %.1f]", q25_alarms_per_provider_per_hour, q75_alarms_per_provider_per_hour),
         c5 = sprintf("%.0f", phys_percent),
         c6 = sprintf("%.0f", tech_percent),
         c7 = sprintf("%.0f [%.0f - %.0f]", median_duration, q25_duration, q75_duration)) %>%
  select(row, label, label2, unit, c1:c7) %>%
  pivot_wider(names_from = unit, values_from = c1:c7, names_glue = "{unit}_{.value}") %>%
  arrange(row) %>% mutate(section = "Alarm burden to the individual nurse")

table3 <- bind_rows(source_rows,
                    tibble(section = "Gestational age", label = c("", "", ""), label2 = ""), ga_rows,
                    tibble(section = "Alarm burden to the individual nurse", label = c("", "", ""), label2 = ""), shift_rows) %>%
  select(section, label, label2, all_of(paste0("picu_c", 1:7)), all_of(paste0("nicu_c", 1:7)))
print(table3, n = Inf, width = Inf)

write_excel_csv(table3, "output//submissionTables//table3.csv", na = "")
