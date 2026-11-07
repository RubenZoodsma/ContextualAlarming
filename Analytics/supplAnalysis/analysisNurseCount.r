## Supplementary analysis: bed-based vs admission-based census per shift, PICU and NICU

library(dplyr)
library(lubridate)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")
data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

## ======== Shifts
shift_days <- seq(as.Date(STUDY_START, tz = "Europe/Amsterdam"), as.Date(STUDY_END, tz = "Europe/Amsterdam"), by = "day")

shifts <- bind_rows(
  tibble(day = shift_days, shift = "Day",     shift_start = as.POSIXct(paste(shift_days, "07:00:00"), tz = "Europe/Amsterdam"),     shift_end = as.POSIXct(paste(shift_days, "15:00:00"), tz = "Europe/Amsterdam")),
  tibble(day = shift_days, shift = "Evening", shift_start = as.POSIXct(paste(shift_days, "15:00:00"), tz = "Europe/Amsterdam"),     shift_end = as.POSIXct(paste(shift_days, "23:00:00"), tz = "Europe/Amsterdam")),
  tibble(day = shift_days, shift = "Night",   shift_start = as.POSIXct(paste(shift_days - 1, "23:00:00"), tz = "Europe/Amsterdam"), shift_end = as.POSIXct(paste(shift_days, "07:00:00"), tz = "Europe/Amsterdam"))) %>%
  mutate(shift_start = pmax(shift_start, STUDY_START), shift_end = pmin(shift_end, STUDY_END)) %>% filter(shift_end > shift_start) %>%
  mutate(shift_hours = as.numeric(difftime(shift_end, shift_start, units = "hours")))

## ======== Census per shift: occupied beds and unique admissions
census <- bed_intervals %>% semi_join(pat_char, by = c("patientID", "afdeling")) %>%
  inner_join(shifts, by = join_by(overlaps(bed_start, bed_end, shift_start, shift_end)), relationship = "many-to-many") %>%
  mutate(present_h = as.numeric(difftime(pmin(bed_end, shift_end), pmax(bed_start, shift_start), units = "hours"))) %>%
  filter(present_h > 0) %>%
  group_by(afdeling, day, shift, shift_hours) %>% summarise(census = n_distinct(bedID), patients = n_distinct(patientID), .groups = 'drop') %>%
  mutate(ratio = case_when(shift == "Day" ~ 1.5, shift == "Evening" ~ 1.75, shift == "Night" ~ 2.0),
         providers_beds = ceiling(census / ratio),
         providers_patients = ceiling(patients / ratio))

## ======== Alarm burden per nurse for both counts
alarm_bins <- data %>% filter(!is.na(patientID)) %>%
  count(afdeling, bin_start = floor_date(alarm_start, "15 minutes"), name = "alarms")

bins <- census %>% select(afdeling, day, shift) %>%
  left_join(shifts %>% select(day, shift, shift_start, shift_end), by = c("day", "shift")) %>%
  reframe(bin_start = seq(shift_start, shift_end - minutes(15), by = "15 min"), .by = c(afdeling, day, shift)) %>%
  left_join(alarm_bins, by = c("afdeling", "bin_start")) %>%
  mutate(alarms = coalesce(alarms, 0L))
rm(alarm_bins)

shift_load <- bins %>% group_by(afdeling, day, shift) %>% summarise(total_alarms = sum(alarms), .groups = 'drop') %>%
  right_join(census, by = c("afdeling", "day", "shift")) %>%
  mutate(total_alarms = coalesce(total_alarms, 0L),
         burden_beds = total_alarms / (providers_beds * shift_hours),
         burden_patients = total_alarms / (providers_patients * shift_hours))
rm(bins)

## ======== Comparison per department and shift
census_comparison <- bind_rows(shift_load, shift_load %>% mutate(shift = "All shifts")) %>%
  group_by(afdeling, shift) %>%
  summarise(
    shifts = n(),
    equal_pct = 100 * mean(patients == census),
    median_beds = median(burden_beds), q25_beds = quantile(burden_beds, 0.25), q75_beds = quantile(burden_beds, 0.75),
    median_patients = median(burden_patients), q25_patients = quantile(burden_patients, 0.25), q75_patients = quantile(burden_patients, 0.75),
    .groups = 'drop') %>%
  mutate(median_change_pct = 100 * (median_patients / median_beds - 1),
         q25_change_pct = 100 * (q25_patients / q25_beds - 1),
         q75_change_pct = 100 * (q75_patients / q75_beds - 1)) %>%
  arrange(afdeling, factor(shift, levels = c("Day", "Evening", "Night", "All shifts")))
print(census_comparison, n = Inf, width = Inf)
