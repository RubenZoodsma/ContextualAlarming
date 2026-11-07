## Appendix II: alarm specifications per unique alarm label, PICU and NICU

library(dplyr)
library(readxl)
library(openxlsx)
options(scipen=999)

load('output//DATA_FINAL.rData')
STUDY_START <- as.POSIXct("2023-08-01 00:00:00", tz = "Europe/Amsterdam")
STUDY_END   <- as.POSIXct("2025-08-01 00:00:00", tz = "Europe/Amsterdam")

data <- data %>% filter(alarm_start >= STUDY_START & alarm_start < STUDY_END)

priority_english <- c("Laag" = "Low", "Medium" = "Medium", "Hoog" = "High")
setdiff(unique(data$alarm_prio_bron), names(priority_english))
specifications <- data %>% group_by(conditie, alarm_prio_bron, device_naam, afdeling,
 Device, Alarm.group, Parameter) %>%
  summarise(Quant_count = n(),
  Quant_duration_median = median(duur_sec, na.rm=TRUE),
  Quant_duration_q25 = quantile(duur_sec, 0.25, na.rm=TRUE),
  Quant_duration_q75 = quantile(duur_sec, 0.75, na.rm=TRUE), .groups = 'drop') %>%
  rename(Alarm_label = conditie, Alarm_priority = alarm_prio_bron, Alarm_device = device_naam, Alarm_department = afdeling,
         Specs_device = Device, Specs_group = Alarm.group, Specs_parameter = Parameter) %>%
  mutate(Alarm_department = case_when(
    Alarm_department == "Neonatologie" ~ "Neonatal intensive care",
    Alarm_department == "Intensive Care Kinderen" ~ "Pediatric intensive care",
    TRUE ~ "Other"),
    Alarm_priority = unname(priority_english[Alarm_priority])) %>%
  select(Alarm_label, Alarm_department, Alarm_priority, Alarm_device,
         Specs_device, Specs_group, Specs_parameter,
         Quant_count, Quant_duration_median, Quant_duration_q25, Quant_duration_q75) %>%
  filter(Quant_count > 0)

specifications %>% count(Alarm_label, Alarm_department, Alarm_priority, Alarm_device) %>% filter(n > 1)
specifications %>% group_by(Alarm_department) %>% summarise(messages = n(), alarms = sum(Quant_count))
specifications %>% count(Alarm_department, name = "unique_messages")
specifications %>% distinct(Alarm_label) %>% nrow()

write.csv(specifications, "output//appendices//appendixII_alarm_specs.csv", row.names = FALSE)

submitted <- read_excel("output//appendices//App_1_AlarmSpecs_revised.xlsx", sheet = "Alarm specifications")

## ======== Write to Excel
wb <- loadWorkbook("output//appendices//App_1_AlarmSpecs_revised.xlsx")
removeTable(wb, sheet = "Alarm specifications", table = "tabel1")
deleteData(wb, sheet = "Alarm specifications", cols = 1:ncol(submitted), rows = 1:(nrow(submitted) + 1), gridExpand = TRUE)
writeDataTable(wb, sheet = "Alarm specifications", x = specifications, startRow = 1, startCol = 1, colNames = TRUE,
               tableName = "tabel1", tableStyle = "TableStyleLight1", withFilter = TRUE, keepNA = FALSE)
saveWorkbook(wb, "output//appendices//appendixII_alarm_specs.xlsx", overwrite = TRUE)
