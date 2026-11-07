## Table 1: baseline characteristics, PICU and NICU

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")

## ======== Denominators

pat_char <- pat_char %>% mutate(
  admission_length_study = as.numeric(difftime(pmin(admission_end, STUDY_END),
                           pmax(admission_start, STUDY_START), units = "hours")))

bed_hours <- bed_intervals %>% semi_join(pat_char, by = c("patientID", "afdeling")) %>%
  mutate(bed_hours = pmax(0, as.numeric(difftime(pmin(bed_end, STUDY_END),
                    pmax(bed_start, STUDY_START), units = "hours")))) %>%
  group_by(afdeling) %>% summarise(bed_hours = sum(bed_hours))

admitted <- pat_char %>% ungroup() %>%
  mutate(start = floor_date(pmax(admission_start, STUDY_START), "1 hour"), end = pmin(admission_end, STUDY_END)) %>% filter(end > start) %>%
  reframe(alarm_floor_hr = seq(start, end - 1, by = "1 hour"), .by = c(patientID, afdeling)) %>%
  left_join(data %>% mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>% count(patientID, afdeling, alarm_floor_hr, name = "alarms_per_hour"),
            by = c("patientID", "afdeling", "alarm_floor_hr")) %>%
  mutate(alarms_per_hour = coalesce(alarms_per_hour, 0L)) %>%
  left_join(pat_char %>% select(patientID, afdeling, diag_group, admission_length_study), by = c("patientID", "afdeling"))

data %>% filter(!is.na(patientID)) %>%
  distinct(afdeling, patientID, alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>%
  anti_join(admitted %>% filter(alarms_per_hour > 0), by = c("afdeling", "patientID", "alarm_floor_hr")) %>%
  count(afdeling)

## ======== Demographics
demo <- pat_char %>% group_by(afdeling) %>% summarise(
  unique_cases = n_distinct(patientID),
  completes = sum(admission_status == "complete_case"),
  unique_pt = n_distinct(hix.id),
  los_hours = sum(admission_length_study, na.rm=TRUE)) %>%
  left_join(bed_hours, by = "afdeling") %>%
  mutate(bed_pct = round(bed_hours / los_hours * 100, 1))
print(demo)

los <- pat_char %>% filter(admission_status == "complete_case") %>% group_by(afdeling) %>%
  summarise(
  median_los = median(admission_length_study, na.rm=TRUE),
  q25_los = quantile(admission_length_study, 0.25, na.rm=TRUE),
  q75_los = quantile(admission_length_study, 0.75, na.rm=TRUE)
)
print(los)

alarms <- data %>% mutate(alarms_all = n(),
  unique_alarms_all = n_distinct(paste(conditie, alarm_prio_bron, device_naam))) %>% group_by(afdeling) %>%
  summarise(total_alarms = n(),
            pct_share = round(total_alarms / alarms_all[1] * 100, 1),
            unique_alarms = n_distinct(paste(conditie, alarm_prio_bron, device_naam)),
            unique_share = round(unique_alarms / unique_alarms_all[1] * 100, 1))
print(alarms)

age <- pat_char %>% ungroup() %>% filter(afdeling=="Intensive Care Kinderen") %>% group_by(afdeling) %>%
  summarise(n_missing = sum(is.na(age_at_admission_years)),
            median_age = median(age_at_admission_years, na.rm=TRUE),
            age_q25 = quantile(age_at_admission_years, 0.25, na.rm=TRUE),
            age_q75 = quantile(age_at_admission_years, 0.75, na.rm=TRUE))
print(age)

ga_summary <- pat_char %>% ungroup() %>% filter(afdeling=="Neonatologie") %>% group_by(afdeling) %>%
  mutate(ga = gestational_age_weeks + coalesce(gestational_age_days, 0) / 7) %>%
  summarise(n_missing = sum(is.na(ga)),
            median_ga = median(ga, na.rm=TRUE),
            ga_q25 = quantile(ga, 0.25, na.rm=TRUE),
            ga_q75 = quantile(ga, 0.75, na.rm=TRUE)) %>%
  mutate(across(c(median_ga, ga_q25, ga_q75), ~ sprintf("%dw%dd", round(.x * 7) %/% 7, round(.x * 7) %% 7), .names = "{.col}_wd"))
print(ga_summary)

## ======== Reason of admission

group_summary <- bind_rows(admitted, admitted %>% mutate(diag_group = "Total")) %>%
  group_by(afdeling, diag_group) %>%
  summarise(
    cases = n_distinct(patientID),
    admitted_hours = n(),
    alarm_active_hours = sum(alarms_per_hour > 0),
    alarm_active_pct = round(alarm_active_hours / admitted_hours * 100),
    median_admitted = median(alarms_per_hour), q25_admitted = quantile(alarms_per_hour, 0.25), q75_admitted = quantile(alarms_per_hour, 0.75),
    median_active = median(alarms_per_hour[alarms_per_hour > 0]), q25_active = quantile(alarms_per_hour[alarms_per_hour > 0], 0.25), q75_active = quantile(alarms_per_hour[alarms_per_hour > 0], 0.75),
    .groups = 'drop') %>%
  left_join(bind_rows(pat_char, pat_char %>% mutate(diag_group = "Total")) %>% group_by(afdeling, diag_group) %>% summarise(los_hours = sum(admission_length_study, na.rm=TRUE), .groups = 'drop'),
            by = c("afdeling", "diag_group")) %>%
  mutate(alarm_active_pct_los = round(alarm_active_hours / los_hours * 100),
         small_group = diag_group != "Total" & cases <= 25) %>%
  arrange(afdeling, diag_group == "Total", diag_group %in% "Other indications", desc(median_active))
print(group_summary, n = Inf, width = Inf)
write.csv(group_summary, "output//primaryOutput//table1_reason_of_admission.csv", row.names = FALSE)

## ======== Shift nurse occupation
shift_days <- seq(as.Date(STUDY_START, tz = "Europe/Amsterdam"), as.Date(STUDY_END, tz = "Europe/Amsterdam"), by = "day")

shifts <- bind_rows(
  tibble(day = shift_days, shift = "Day",     shift_start = as.POSIXct(paste(shift_days, "07:00:00"), tz = "Europe/Amsterdam"),     shift_end = as.POSIXct(paste(shift_days, "15:00:00"), tz = "Europe/Amsterdam")),
  tibble(day = shift_days, shift = "Evening", shift_start = as.POSIXct(paste(shift_days, "15:00:00"), tz = "Europe/Amsterdam"),     shift_end = as.POSIXct(paste(shift_days, "23:00:00"), tz = "Europe/Amsterdam")),
  tibble(day = shift_days, shift = "Night",   shift_start = as.POSIXct(paste(shift_days - 1, "23:00:00"), tz = "Europe/Amsterdam"), shift_end = as.POSIXct(paste(shift_days, "07:00:00"), tz = "Europe/Amsterdam"))) %>%
  mutate(shift_start = pmax(shift_start, STUDY_START), shift_end = pmin(shift_end, STUDY_END)) %>% filter(shift_end > shift_start) %>%
  mutate(shift_hours = as.numeric(difftime(shift_end, shift_start, units = "hours")))

uniqueBed_Shift <- bed_intervals %>% semi_join(pat_char, by = c("patientID", "afdeling")) %>%
  inner_join(shifts, by = join_by(overlaps(bed_start, bed_end, shift_start, shift_end)), relationship = "many-to-many") %>%
  mutate(present_h = as.numeric(difftime(pmin(bed_end, shift_end), pmax(bed_start, shift_start), units = "hours"))) %>%
  filter(present_h > 0) %>%
  group_by(afdeling, day, shift) %>% summarise(uniqueBed_Shift = n_distinct(bedID), .groups = 'drop')

nurses <- uniqueBed_Shift %>%
  mutate(shift = factor(shift, levels = c("Day", "Evening", "Night")),
         providers = case_when(
           shift == "Day" ~ ceiling(uniqueBed_Shift / 1.5),
           shift == "Evening" ~ ceiling(uniqueBed_Shift / 1.75),
           shift == "Night" ~ ceiling(uniqueBed_Shift / 2.0)))
nurse_summary <- bind_rows(nurses, nurses %>% mutate(shift = "Average")) %>%
  group_by(afdeling, shift) %>%
  summarise(mean = round(mean(providers), 1), median = median(providers),
            q25 = quantile(providers, 0.25), q75 = quantile(providers, 0.75), .groups = 'drop') %>%
  arrange(afdeling, match(shift, c("Day", "Evening", "Night", "Average")))
print(nurse_summary)

## ======== Table 1 layout
demographics <- demo %>% left_join(los, by = "afdeling") %>% left_join(alarms, by = "afdeling") %>%
  left_join(group_summary %>% filter(diag_group == "Total") %>% select(afdeling, admitted_hours, alarm_active_hours, alarm_active_pct), by = "afdeling") %>%
  left_join(age, by = "afdeling") %>% left_join(ga_summary, by = "afdeling") %>%
  transmute(afdeling,
    `Cases, n` = as.character(unique_cases),
    `Unique patients, n` = as.character(unique_pt),
    `Admitted clock hours, n` = as.character(admitted_hours),
    `Length-of-stay, hours (median [IQR])` = sprintf("%.1f [%.1f - %.1f]", median_los, q25_los, q75_los),
    `Total alarm-active hours, n [%]` = sprintf("%d [%.0f]", alarm_active_hours, alarm_active_pct),
    `Total alarms, n [%]` = sprintf("%d [%.0f]", total_alarms, pct_share),
    `Unique alarm labels, n` = as.character(unique_alarms),
    `Age at admission, years, median [IQR]` = if_else(is.na(median_age), "-", sprintf("%.1f [%.1f - %.1f]", median_age, age_q25, age_q75)),
    `Gestational age, median [IQR]` = if_else(is.na(median_ga), "-", sprintf("%s [%s - %s]", median_ga_wd, ga_q25_wd, ga_q75_wd))) %>%
  pivot_longer(-afdeling, names_to = "picu_label") %>%
  pivot_wider(names_from = afdeling, values_from = value) %>%
  transmute(section = "Demographics", picu_label, picu_c1 = `Intensive Care Kinderen`, nicu_label = NA_character_, nicu_c1 = Neonatologie)

reason <- group_summary %>% filter(!is.na(diag_group)) %>%
  transmute(afdeling,
    label = paste0(diag_group, if_else(small_group, " †", "")),
    c1 = as.character(cases),
    c2 = sprintf("%d [%.0f]", admitted_hours, alarm_active_pct),
    c3 = sprintf("%.0f [%.0f - %.0f]", median_active, q25_active, q75_active),
    c4 = sprintf("%.0f [%.0f - %.0f]", median_admitted, q25_admitted, q75_admitted))
reason <- bind_cols(reason %>% filter(afdeling == "Intensive Care Kinderen") %>% select(-afdeling) %>% rename_with(~ paste0("picu_", .x)),
                    reason %>% filter(afdeling == "Neonatologie") %>% select(-afdeling) %>% rename_with(~ paste0("nicu_", .x))) %>%
  mutate(section = "Reason of admission", .before = 1)

shift <- nurse_summary %>%
  transmute(afdeling,
    label = case_when(shift == "Day" ~ "Day, 07h - 15h", shift == "Evening" ~ "Evening, 15h - 23h", shift == "Night" ~ "Night, 23h - 07h", TRUE ~ "Average"),
    c1 = case_when(shift == "Day" ~ "1 : 1.5", shift == "Evening" ~ "1 : 1.75", shift == "Night" ~ "1 : 2", TRUE ~ ""),
    c2 = sprintf("%.1f", mean),
    c3 = sprintf("%.0f", median),
    c4 = sprintf("[%.0f - %.0f]", q25, q75))
shift <- bind_cols(shift %>% filter(afdeling == "Intensive Care Kinderen") %>% select(-afdeling) %>% rename_with(~ paste0("picu_", .x)),
                   shift %>% filter(afdeling == "Neonatologie") %>% select(-afdeling) %>% rename_with(~ paste0("nicu_", .x))) %>%
  mutate(section = "Shift nurse occupation", .before = 1)

table1 <- bind_rows(demographics, reason, shift) %>%
  select(section, picu_label, picu_c1, picu_c2, picu_c3, picu_c4, nicu_label, nicu_c1, nicu_c2, nicu_c3, nicu_c4)
print(table1, n = Inf, width = Inf)

write_excel_csv(table1, "output//submissionTables//table1.csv", na = "")
