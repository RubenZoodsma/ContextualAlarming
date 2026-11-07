## Appendix VI: alarm burden per nurse over the course of the day, by nurse:patient ratio, PICU and NICU

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

## ======== Shifts and nurses per ratio scenario
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

scenario_levels <- c("1:1", "1:1.5", "1:2", "Institutional (1:1.5 / 1:1.75 / 1:2)")
nurses <- bind_rows(
  census %>% mutate(scenario = "1:1",   providers = ceiling(census / 1)),
  census %>% mutate(scenario = "1:1.5", providers = ceiling(census / 1.5)),
  census %>% mutate(scenario = "1:2",   providers = ceiling(census / 2)),
  census %>% mutate(scenario = "Institutional (1:1.5 / 1:1.75 / 1:2)",
                    providers = case_when(
                      shift == "Day" ~ ceiling(census / 1.5),
                      shift == "Evening" ~ ceiling(census / 1.75),
                      shift == "Night" ~ ceiling(census / 2.0))))

## ======== Alarms per 15-minute bin of every shift
alarm_bins <- data %>% filter(!is.na(patientID)) %>%
  count(afdeling, bin_start = floor_date(alarm_start, "15 minutes"), name = "alarms")

bins <- census %>% select(afdeling, day, shift) %>%
  left_join(shifts %>% select(day, shift, shift_start, shift_end), by = c("day", "shift")) %>%
  reframe(bin_start = seq(shift_start, shift_end - minutes(15), by = "15 min"), .by = c(afdeling, day, shift)) %>%
  left_join(alarm_bins, by = c("afdeling", "bin_start")) %>%
  mutate(alarms = coalesce(alarms, 0L),
         time_of_day = hour(bin_start) + minute(bin_start) / 60)
rm(alarm_bins)

shift_load <- bins %>% group_by(afdeling, day, shift) %>% summarise(total_alarms = sum(alarms), .groups = 'drop') %>%
  right_join(nurses, by = c("afdeling", "day", "shift")) %>%
  mutate(total_alarms = coalesce(total_alarms, 0L),
         alarms_per_provider_per_hour = total_alarms / (providers * shift_hours))
burden_by_ratio <- bind_rows(shift_load, shift_load %>% mutate(shift = "Average")) %>%
  group_by(afdeling, scenario, shift) %>%
  summarise(
    shifts = n(),
    median_providers = median(providers),
    mean_alarms_per_provider_per_hour = round(mean(alarms_per_provider_per_hour, na.rm = TRUE), 1),
    median_alarms_per_provider_per_hour = round(median(alarms_per_provider_per_hour, na.rm = TRUE), 1),
    q25_alarms_per_provider_per_hour = round(quantile(alarms_per_provider_per_hour, 0.25, na.rm = TRUE), 1),
    q75_alarms_per_provider_per_hour = round(quantile(alarms_per_provider_per_hour, 0.75, na.rm = TRUE), 1),
    .groups = 'drop') %>%
  arrange(afdeling, factor(shift, levels = c("Day", "Evening", "Night", "Average")), factor(scenario, levels = scenario_levels))
print(burden_by_ratio, n = Inf, width = Inf)
write.csv(burden_by_ratio, "output//appendices//appendixVI_burden_by_ratio.csv", row.names = FALSE)

## ======== Burden per nurse over the day
burden_over_day <- bins %>%
  left_join(nurses %>% select(afdeling, day, shift, scenario, providers), by = c("afdeling", "day", "shift"), relationship = "many-to-many") %>%
  mutate(alarms_per_provider_per_hour = alarms * 4 / providers) %>%
  group_by(afdeling, scenario, time_of_day) %>%
  summarise(
    bins = n(),
    mean_alarms_per_provider_per_hour = mean(alarms_per_provider_per_hour),
    median_alarms_per_provider_per_hour = median(alarms_per_provider_per_hour),
    q25_alarms_per_provider_per_hour = quantile(alarms_per_provider_per_hour, 0.25),
    q75_alarms_per_provider_per_hour = quantile(alarms_per_provider_per_hour, 0.75),
    .groups = 'drop') %>%
  mutate(bin_mid = time_of_day + 0.125)
write.csv(burden_over_day, "output//appendices//appendixVI_burden_over_day.csv", row.names = FALSE)
rm(bins)

shift_steps <- burden_over_day %>%
  filter(time_of_day %in% c(6.75, 7, 14.75, 15, 22.75, 23)) %>%
  mutate(shift_change = case_when(time_of_day %in% c(6.75, 7) ~ "07:00", time_of_day %in% c(14.75, 15) ~ "15:00", .default = "23:00"),
         side = if_else(time_of_day %in% c(7, 15, 23), "after", "before")) %>%
  select(afdeling, scenario, shift_change, side, median_alarms_per_provider_per_hour) %>%
  pivot_wider(names_from = side, values_from = median_alarms_per_provider_per_hour) %>%
  mutate(step_pct = round((after - before) / before * 100)) %>%
  arrange(afdeling, shift_change, factor(scenario, levels = scenario_levels))
print(shift_steps, n = Inf, width = Inf)
write.csv(shift_steps, "output//appendices//appendixVI_shift_steps.csv", row.names = FALSE)

## ======== Plot
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
shift_labels <- data.frame(x = c(3.5, 11, 19),
                           label = c("Night", "Day", "Evening"))
scenario_rows <- c("1:1" = "Uniform 1:1", "1:1.5" = "Uniform 1:1.5", "1:2" = "Uniform 1:2",
                   "Institutional (1:1.5 / 1:1.75 / 1:2)" = "Institutional\n(1:1.5 / 1:1.75 / 1:2)")
row_levels <- unname(scenario_rows[c(4, 3, 2, 1)])
y_max <- max(burden_over_day$q75_alarms_per_provider_per_hour)
scenario_reference <- burden_by_ratio %>% filter(shift == "Average") %>%
  mutate(department = factor(unit_panels[afdeling], levels = unit_panels), row = factor(scenario_rows[scenario], levels = row_levels))
plot_data <- burden_over_day %>%
  mutate(department = factor(unit_panels[afdeling], levels = unit_panels), row = factor(scenario_rows[scenario], levels = row_levels))

plot_top <- plot_data %>% filter(row == row_levels[1]) %>%
  ggplot(aes(x = bin_mid)) +
  geom_ribbon(aes(ymin = q25_alarms_per_provider_per_hour, ymax = q75_alarms_per_provider_per_hour), fill = "#386cb0", alpha = 0.12) +
  geom_line(aes(y = q25_alarms_per_provider_per_hour), color = "#386cb0", linewidth = 0.5) +
  geom_line(aes(y = q75_alarms_per_provider_per_hour), color = "#386cb0", linewidth = 0.5) +
  geom_line(aes(y = median_alarms_per_provider_per_hour), color = "grey20", linewidth = 0.7) +
  geom_vline(xintercept = c(7, 15, 23), linetype = "dashed", color = "grey50") +
  geom_hline(data = scenario_reference %>% filter(row == row_levels[1]), aes(yintercept = median_alarms_per_provider_per_hour), linetype = "dashed", color = "black") +
  geom_label(data = scenario_reference %>% filter(row == row_levels[1]), aes(x = 0.3, y = median_alarms_per_provider_per_hour,
                                            label = paste0("All shifts: median = ", median_alarms_per_provider_per_hour)),
             hjust = 0, vjust = -0.3, size = 2.6, linewidth = 0, fill = "white", label.r = unit(0, "pt")) +
  geom_text(data = shift_labels, aes(x = x, y = Inf, label = label), vjust = 1.5, fontface = "bold.italic", size = 2.8) +
  facet_grid(row ~ department) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(0, 200, 20), limits = c(0, y_max), expand = expansion(mult = c(0.02, 0.12))) +
  theme_Publication(base_size = 11) +
  theme(panel.spacing.x = unit(0.6, "cm"), plot.margin = unit(c(10, 5, 12, 5), "mm")) +
  xlab("Time of day") +
  ylab("Alarms per nurse per hour")

plot_bottom <- plot_data %>% filter(row != row_levels[1]) %>%
  ggplot(aes(x = bin_mid)) +
  geom_ribbon(aes(ymin = q25_alarms_per_provider_per_hour, ymax = q75_alarms_per_provider_per_hour), fill = "#386cb0", alpha = 0.12) +
  geom_line(aes(y = q25_alarms_per_provider_per_hour), color = "#386cb0", linewidth = 0.5) +
  geom_line(aes(y = q75_alarms_per_provider_per_hour), color = "#386cb0", linewidth = 0.5) +
  geom_line(aes(y = median_alarms_per_provider_per_hour), color = "grey20", linewidth = 0.7) +
  geom_vline(xintercept = c(7, 15, 23), linetype = "dashed", color = "grey50") +
  geom_hline(data = scenario_reference %>% filter(row != row_levels[1]), aes(yintercept = median_alarms_per_provider_per_hour), linetype = "dashed", color = "black") +
  geom_label(data = scenario_reference %>% filter(row != row_levels[1]), aes(x = 0.3, y = median_alarms_per_provider_per_hour,
                                            label = paste0("All shifts: median = ", median_alarms_per_provider_per_hour)),
             hjust = 0, vjust = -0.3, size = 2.6, linewidth = 0, fill = "white", label.r = unit(0, "pt")) +
  geom_text(data = shift_labels, aes(x = x, y = Inf, label = label), vjust = 1.5, fontface = "bold.italic", size = 2.8) +
  facet_grid(row ~ department) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(0, 200, 20), limits = c(0, y_max), expand = expansion(mult = c(0.02, 0.12))) +
  theme_Publication(base_size = 11) +
  theme(panel.spacing.x = unit(0.6, "cm"), panel.spacing.y = unit(0.4, "cm"), plot.margin = unit(c(0, 5, 5, 5), "mm")) +
  xlab("Time of day") +
  ylab("Alarms per nurse per hour")

plot_top <- ggplotGrob(plot_top)
institutional_strip <- which(plot_top$layout$name == "strip-r-1")
plot_top$grobs[[institutional_strip]]$grobs[[1]]$children[[1]]$gp$fill <- "#2c2c2c"
plot_top$grobs[[institutional_strip]]$grobs[[1]]$children[[1]]$gp$col <- "#2c2c2c"
plot_top$grobs[[institutional_strip]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- "white"
plot_bottom <- ggplotGrob(plot_bottom)
plot <- rbind(plot_top, plot_bottom, size = "max")
plot$heights[plot$layout$t[grepl("^panel", plot$layout$name)]] <- unit(1, "null")
ggsave("output//appendices//appendixVI_burden_over_day.png", plot, width = 30, height = 34, units = "cm", dpi = 300)
