## Appendix III: top-quartile admissions by alarm volume and by alarm rate, PICU and NICU

library(dplyr)
library(tidyr)
library(lubridate)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

## ======== Admitted clock hours
pat_char <- pat_char %>% mutate(
  admission_length_study = as.numeric(difftime(pmin(admission_end, STUDY_END),
                           pmax(admission_start, STUDY_START), units = "hours")))

admitted <- pat_char %>% ungroup() %>%
  mutate(start = floor_date(pmax(admission_start, STUDY_START), "1 hour"), end = pmin(admission_end, STUDY_END)) %>% filter(end > start) %>%
  reframe(alarm_floor_hr = seq(start, end - 1, by = "1 hour"), .by = c(patientID, afdeling)) %>%
  left_join(data %>% mutate(alarm_floor_hr = floor_date(alarm_start, "1 hour")) %>% count(patientID, afdeling, alarm_floor_hr, name = "alarms_per_hour"),
            by = c("patientID", "afdeling", "alarm_floor_hr")) %>%
  mutate(alarms_per_hour = coalesce(alarms_per_hour, 0L))

## ======== Volume and rate per admission
admissions <- pat_char %>% ungroup() %>%
  select(afdeling, patientID, hix.id, admission_status, admission_length_study, diag_group,
         age_at_admission_years, gestational_age_weeks, gestational_age_days, GestAgeGroup) %>%
  left_join(data %>% group_by(afdeling, patientID) %>% summarise(total_alarms = n(), .groups = 'drop'), by = c("afdeling", "patientID")) %>%
  left_join(admitted %>% group_by(afdeling, patientID) %>%
              summarise(alarm_active_hours = sum(alarms_per_hour > 0), median_rate = median(alarms_per_hour[alarms_per_hour > 0]), .groups = 'drop'),
            by = c("afdeling", "patientID")) %>%
  mutate(total_alarms = coalesce(total_alarms, 0L), alarm_active_hours = coalesce(alarm_active_hours, 0L)) %>%
  group_by(afdeling) %>%
  mutate(top_volume = total_alarms >= quantile(total_alarms, 0.75),
         top_rate = coalesce(median_rate >= quantile(median_rate, 0.75, na.rm=TRUE), FALSE)) %>%
  ungroup()

admissions %>% group_by(afdeling) %>%
  summarise(admissions = n(), no_alarms = sum(total_alarms == 0), no_rate = sum(is.na(median_rate)),
            volume_threshold = quantile(total_alarms, 0.75), rate_threshold = quantile(median_rate, 0.75, na.rm=TRUE),
            alarm_responsibility_pct = sum(total_alarms[top_volume]) / sum(total_alarms) * 100)
data %>% count(afdeling, name = "alarms_all") %>%
  left_join(data %>% semi_join(pat_char, by = c("patientID", "afdeling")) %>% count(afdeling, name = "alarms_admissions"), by = "afdeling") %>%
  left_join(admitted %>% group_by(afdeling) %>% summarise(alarms_admitted_hours = sum(alarms_per_hour)), by = "afdeling")

admissions_cat <- bind_rows(admissions %>% mutate(category = "Volume", top_quartile = top_volume),
                            admissions %>% mutate(category = "Rate", top_quartile = top_rate))

alarms_cat <- data %>% select(afdeling, patientID, Device, Alarm.group, conditie) %>%
  inner_join(admissions_cat %>% select(afdeling, patientID, category, top_quartile), by = c("afdeling", "patientID"), relationship = "many-to-many")

small_groups <- admitted %>% left_join(pat_char %>% select(patientID, afdeling, diag_group), by = c("patientID", "afdeling")) %>%
  group_by(afdeling, diag_group) %>% summarise(cases = n_distinct(patientID), .groups = 'drop') %>%
  mutate(small_group = cases <= 25)

## ======== Characteristics: top quartile vs remaining admissions
ga_levels <- c("Extremely preterm", "Very preterm", "Moderate to late preterm", "Term")
characteristics <- bind_rows(
  admissions_cat %>% group_by(afdeling, category, top_quartile) %>% summarise(n = n(), .groups = 'drop') %>%
    group_by(afdeling, category) %>% mutate(value = sprintf("%s (%.1f)", format(n, big.mark = ",", trim = TRUE), n / sum(n) * 100)) %>% ungroup() %>%
    mutate(row = "Admissions, n (% of unit)", ord = 1),
  admissions_cat %>% group_by(afdeling, category, top_quartile) %>% summarise(n = sum(total_alarms), .groups = 'drop') %>%
    group_by(afdeling, category) %>% mutate(value = sprintf("%s (%.1f)", format(n, big.mark = ",", trim = TRUE), n / sum(n) * 100)) %>% ungroup() %>%
    mutate(row = "Alarms, n (% of unit)", ord = 2),
  admissions_cat %>% group_by(afdeling, category, top_quartile) %>%
    summarise(value = sprintf("%s [%s-%s]", formatC(median(alarm_active_hours), format = "f", digits = 0, big.mark = ","),
                              formatC(quantile(alarm_active_hours, 0.25), format = "f", digits = 0, big.mark = ","),
                              formatC(quantile(alarm_active_hours, 0.75), format = "f", digits = 0, big.mark = ",")), .groups = 'drop') %>%
    mutate(row = "Alarm-active hours per admission, median [IQR]", ord = 3),
  admissions_cat %>% group_by(afdeling, category, top_quartile) %>%
    summarise(value = sprintf("%.1f [%.1f-%.1f]", median(median_rate, na.rm=TRUE), quantile(median_rate, 0.25, na.rm=TRUE),
                              quantile(median_rate, 0.75, na.rm=TRUE)), .groups = 'drop') %>%
    mutate(row = "Alarms per alarm-active hour, median [IQR]", ord = 4),
  admissions_cat %>% filter(admission_status == "complete_case") %>% group_by(afdeling, category, top_quartile) %>%
    summarise(value = sprintf("%.0f [%.0f-%.0f]", median(admission_length_study, na.rm=TRUE), quantile(admission_length_study, 0.25, na.rm=TRUE),
                              quantile(admission_length_study, 0.75, na.rm=TRUE)), .groups = 'drop') %>%
    mutate(row = "Length-of-stay, hours, median [IQR]", ord = 5),
  admissions_cat %>% filter(afdeling == "Intensive Care Kinderen") %>% group_by(afdeling, category, top_quartile) %>%
    summarise(value = sprintf("%.1f [%.1f-%.1f]", median(age_at_admission_years, na.rm=TRUE), quantile(age_at_admission_years, 0.25, na.rm=TRUE),
                              quantile(age_at_admission_years, 0.75, na.rm=TRUE)), .groups = 'drop') %>%
    mutate(row = "Age, years, median [IQR]", ord = 6),
  admissions_cat %>% filter(afdeling == "Neonatologie") %>%
    mutate(ga = gestational_age_weeks + coalesce(gestational_age_days, 0) / 7) %>% group_by(afdeling, category, top_quartile) %>%
    summarise(median_ga = median(ga, na.rm=TRUE), ga_q25 = quantile(ga, 0.25, na.rm=TRUE), ga_q75 = quantile(ga, 0.75, na.rm=TRUE), .groups = 'drop') %>%
    mutate(across(c(median_ga, ga_q25, ga_q75), ~ sprintf("%dw%dd", round(.x * 7) %/% 7, round(.x * 7) %% 7), .names = "{.col}_wd"),
           value = sprintf("%s [%s-%s]", median_ga_wd, ga_q25_wd, ga_q75_wd), row = "Gestational age, median [IQR]", ord = 7),
  admissions_cat %>% filter(afdeling == "Neonatologie") %>%
    mutate(level = coalesce(as.character(GestAgeGroup), "Missing")) %>% count(afdeling, category, top_quartile, level) %>%
    complete(nesting(afdeling, level), nesting(category, top_quartile), fill = list(n = 0L)) %>%
    group_by(afdeling, category, top_quartile) %>% mutate(value = sprintf("%s (%.1f)", format(n, big.mark = ",", trim = TRUE), n / sum(n) * 100)) %>% ungroup() %>%
    mutate(row = paste("Gestational age group:", level), ord = 8 + coalesce(match(level, ga_levels), 5L) / 10),
  admissions_cat %>% count(afdeling, category, top_quartile, diag_group) %>%
    complete(nesting(afdeling, diag_group), nesting(category, top_quartile), fill = list(n = 0L)) %>%
    group_by(afdeling, category, top_quartile) %>% mutate(value = sprintf("%s (%.1f)", format(n, big.mark = ",", trim = TRUE), n / sum(n) * 100)) %>% ungroup() %>%
    left_join(small_groups %>% select(afdeling, diag_group, small_group), by = c("afdeling", "diag_group")) %>%
    mutate(row = paste0("Reason of admission: ", diag_group, ifelse(coalesce(small_group, TRUE), " †", "")), ord = 9),
  alarms_cat %>% group_by(afdeling, category, top_quartile) %>%
    summarise(Physiological = round(sum(Alarm.group == "Physiological", na.rm = TRUE) / n() * 100, 0),
              Technical = round(sum(Alarm.group == "Technical", na.rm = TRUE) / n() * 100, 0), .groups = 'drop') %>%
    pivot_longer(c(Physiological, Technical), names_to = "level", values_to = "pct") %>%
    mutate(value = sprintf("%.0f", pct), row = paste("Alarm nature, %:", level), ord = 10),
  alarms_cat %>% count(afdeling, category, top_quartile, Device) %>%
    complete(nesting(afdeling, Device), nesting(category, top_quartile), fill = list(n = 0L)) %>%
    group_by(afdeling, category, top_quartile) %>% mutate(value = sprintf("%.1f", n / sum(n) * 100)) %>% ungroup() %>%
    mutate(row = paste("Device, %:", Device), ord = 11)) %>%
  mutate(column = paste(c("Intensive Care Kinderen" = "PICU", "Neonatologie" = "NICU")[afdeling], category,
                        ifelse(top_quartile, "top", "remaining"), sep = "_")) %>%
  select(ord, row, column, value) %>%
  pivot_wider(names_from = column, values_from = value) %>%
  select(row, PICU_Volume_top, PICU_Volume_remaining, PICU_Rate_top, PICU_Rate_remaining,
         NICU_Volume_top, NICU_Volume_remaining, NICU_Rate_top, NICU_Rate_remaining, ord) %>%
  arrange(ord, row) %>% select(-ord)
print(characteristics, n = Inf, width = Inf)
write.csv(characteristics, "output//supplOutput//appendixIII_A_characteristics.csv", row.names = FALSE, na = "")

## ======== 25 most common alarm labels per alarm nature, top quartile by volume
labels <- alarms_cat %>% filter(category == "Volume", Alarm.group %in% c("Physiological", "Technical")) %>%
  count(afdeling, top_quartile, Alarm.group, Device, conditie) %>%
  group_by(afdeling, top_quartile, Alarm.group) %>%
  mutate(pct = round(n / sum(n) * 100, 1), rank = rank(-n, ties.method = "first")) %>% ungroup()

top25 <- labels %>% filter(top_quartile, rank <= 25) %>%
  select(afdeling, Alarm.group, rank_top = rank, Device, conditie, n_top = n, pct_top = pct) %>%
  left_join(labels %>% filter(!top_quartile) %>% select(afdeling, Alarm.group, Device, conditie, n_remaining = n, pct_remaining = pct, rank_remaining = rank),
            by = c("afdeling", "Alarm.group", "Device", "conditie")) %>%
  mutate(n_remaining = coalesce(n_remaining, 0L), pct_remaining = coalesce(pct_remaining, 0),
         pct_diff = round(pct_top - pct_remaining, 1)) %>%
  arrange(afdeling, Alarm.group, rank_top)
print(top25, n = Inf, width = Inf)
write.csv(top25, "output//supplOutput//appendixIII_B_top25_labels.csv", row.names = FALSE, na = "")
