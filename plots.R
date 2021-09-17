# # # # # # # # # # # # # # # # # # # # 
# Script to produce occupancy figures
# 
# David J Price
# 19 July 2021
# 

## Inputs:
# scenario.folder name: folder in current working directory containing the output 
#                       from the clinical pathways model
# quantium_scenario_target_dates.csv: file containing dates that vaccination targets 
#                                     are anticipated to be reached under different 
#                                     vaccination strategies
# dim_scenario.csv: list containing details of vaccination strategies considered 
#                   by Quantium

## Outputs:
# Creates a folder for each model output file (if it doesn't already exist)
# Saves figures for each of: all infections, symptomatic infections, deaths, ward occupancy, ICU occupancy, absenteeism 


require(tidyverse)
require(lubridate)
require(R.matlab)
source("./functions.R")


# Point to folder with output from clinical pathway model
scenario.folder <- "./scenarios/"

# Specify age groups
age.groups <- paste(seq(0,80, by=5), seq(0,80, by=5)+5, sep = "_")
age.groups[17] <- "80+"

# Percentage of each age group in workforce
perc.workforce <- c(0, 0, 0, 2 * 0.27006, 0.83684, 0.83684, 0.862148, 0.862148, 0.863779, 
                    0.863779, 0.799347, 0.799347, 0.46397, 0.46397, 0.115252, 0.115252, 0.02064) # Provided by Treasury

pc.workforce <- data.frame(age.groups, "pc"=perc.workforce)

# Specify plotting colours and transparencies
infections.colour <- "steelblue3"
symp.cases.colour <- "#E0D03A"
death.colour <- "#FC8D62"

ed.excess.colour <- "grey40"

ward.colour <- "#BA53C9"
icu.colour <- "#66c2a5"

outer.alpha <- 0.15
inner.alpha <- 0.4
line.alpha <- 0.8


## Get dates to plot on top ----
the.dates <- read_csv("quantium_scenario_target_dates.csv", 
                      col_types = list(
                        `50%` = col_date(format = "%d/%m/%y"),
                        `60%` = col_date(format = "%d/%m/%y"),
                        `70%` = col_date(format = "%d/%m/%y"),
                        `80%` = col_date(format = "%d/%m/%y")
                      ))

## Get strategy details ----
lookup.strategy <- read_csv("dim_scenario.csv") %>% 
  mutate(priority_order = case_when(
    priority_order == "Oldest to youngest" ~ "Oldest",
    priority_order == "Youngest to oldest (40+ first then 16+)" ~ "Middle",
    priority_order == "Random" ~ "AllAdults",
    priority_order == "Phased" ~ "Phased"
  ))


# Specify capacities (for horizontal lines) and plot limits
ward.capacity <- 25756
icu.capacity <- 1848

inf.upper.y <- 122500
death.upper.y <- 180
ward.upper.y <- max(ward.capacity,10000)
icu.upper.y <- max(icu.capacity, 3200)
deaths.upper.y <- 200

max.day <- 180
my.day.breaks <- seq(0, max.day, by = 30)

# Get all files in current folder
filenames <- list.files(path = scenario.folder)


for (filename in filenames){
  print(filename)
  dat <- readMat(paste0(scenario.folder,filename))
  
  folder.name <- gsub(pattern = ".mat", replacement = "", x = gsub(pattern = "MOCsimulations_", replacement = "", x = filename))
  dir.create(folder.name, showWarnings = F)
  
  #the.dates
  
  tmp <- gsub(pattern = "([0-9]+).*$", replacement = "\\1", x = filename)
  coverage <- as.numeric(gsub(pattern = "MOCsimulations_coverage", replacement = "", x = tmp))
  
  tmp <- gsub(pattern = "MOCsimulations_coverage[0-9]+_rollout", replacement = "", x = filename)
  rollout <- gsub(pattern = "_VOCDelta_optimalTTIQ_[0-9]+.mat", replacement = "", x = tmp)
  rollout <- gsub(pattern = "_VOCDelta_partial_lowphsmTTIQ_[0-9]+.mat", replacement = "", x = rollout)
  rollout <- gsub(pattern = "_VOCDelta_partial_LMLTTIQ_[0-9]+.mat", replacement = "", x = rollout)
  rollout <- gsub(pattern = "_VOCDelta_partialTTIQ_[0-9]+.mat", replacement = "", x = rollout) # add in removing number from filename
  
  # lookup strategy
  
  az.dose.gap <- "12 weeks"
  az.age.cutoff <- "60"
  
  vacc.start.date <- ymd("2020-11-30")
  
  the.scenario <- the.dates %>% 
    slice(lookup.strategy %>% 
            filter(priority_order == rollout, 
                   az_dose_gap == az.dose.gap, 
                   az_age_cutoff == az.age.cutoff) %>% 
            select(scenario) %>% pull)
  

  vlines <- the.scenario %>% pivot_longer(cols = `50%`:`80%`) %>% 
    mutate(cov = as.numeric(gsub(pattern = "%",replacement = "",x = name))) %>% 
    filter(cov > coverage) %>% 
    mutate(value = as.numeric(value - vacc.start.date)) %>% 
    select(name, value)
  
  
  # Combine all vax/unvax from each scenario here
  df.az1 <- convert.mat.no.age(dat$AZ1.data, "AZ1")
  df.az2 <- convert.mat.no.age(dat$AZ2.data, "AZ2")
  df.M1 <- convert.mat.no.age(dat$M1.data, "M1")
  df.M2 <- convert.mat.no.age(dat$M2.data, "M2")
  df.P1 <- convert.mat.no.age(dat$P1.data, "P1")
  df.P2 <- convert.mat.no.age(dat$P2.data, "P2")
  df.unvax <- convert.mat.no.age(dat$unvax.data, "unvax")
  
  
  df <- df.az1 %>% select(-group) %>% 
    left_join(., df.az2 %>% select(-group), by = c("var","sim","day")) %>%
    left_join(., df.M1 %>% select(-group), by = c("var","sim","day")) %>% 
    left_join(., df.M2 %>% select(-group), by = c("var","sim","day")) %>% 
    left_join(., df.P1 %>% select(-group), by = c("var","sim","day")) %>% 
    left_join(., df.P2 %>% select(-group), by = c("var","sim","day")) %>% 
    left_join(., df.unvax %>% select(-group), by = c("var","sim","day")) %>% 
    mutate(value = value.x + value.y + value.x.x + value.y.y + value.x.x.x + value.y.y.y + value) %>% 
    select(-value.x, -value.y, -value.x.x, -value.y.y, -value.x.x.x, -value.y.y.y)
  
  # Add asymptomatic and symptomatic together for 'total cases'

  # Specify day of seeding
  if (coverage == 70){
    start.date <- 342 + 1
  } else if (coverage == 80){
    start.date <- 358 + 1
  }
  
  # Align vertical lines for vaccination targets
  vlines <- vlines %>% 
    mutate(value = value - start.date)
  
  
  df <- df %>% 
    filter(day >= 0, day <= max.day) %>% 
    mutate(day = day )
  
  
  
  # inf.traj.plot <- df %>% 
  #   filter(var == "infections") %>% 
  #   ggplot() +
  #   aes(x = day, y = value, group = sim) +
  #   
  #   geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
  #   annotate(geom = "text", x = vlines$value+1, y = 
  #              (df %>% summarise(m = max(value)) %>% pull())*0.95, 
  #            label = vlines$name, colour = "grey70",
  #            hjust = "left") +
  #   
  #   geom_line(alpha = 0.25) +
  #   scale_x_continuous("Day", breaks = my.day.breaks) + 
  #   scale_y_continuous("Daily New Infections") + 
  #   cowplot::theme_cowplot() +
  #   theme(text = element_text(size = 16))
  # 
  # ggsave(filename = paste0(folder.name,"/",folder.name, "_infections_trajectory.png"), plot = inf.traj.plot, height = 4, width = 6, units="in", bg = "white")
  
  
  # Generate quantiles for plotting
  plot.dat <- df %>% 
    group_by(var, day) %>% 
    summarise(x=list(enframe(
      quantile(value, probs=c(0.05, 
                              #0.1, 0.15, 0.2, 
                              0.25, 0.75, 
                              #0.8, 0.85, 0.9, 
                              0.95)),"quantiles", "value"))) %>%
    unnest(x) %>% mutate(quantiles = as.numeric(gsub(pattern = "%",replacement = "",x = quantiles))) %>%
    mutate(quantiles = paste0("q",quantiles)) %>% pivot_wider(names_from = quantiles, values_from=value)
  
  
  
  ## Incidence ----
  inf.plot <- plot.dat %>% 
    filter(var == "infections") %>% 
    ggplot() + 
    aes(x = day) +
    geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
    annotate(geom = "text", x = vlines$value+1, 
             # y = (plot.dat %>% filter(var == "infections") %>% summarise(m = max(q95)) %>% pull()), 
             y = inf.upper.y*0.95,
             label = vlines$name, colour = "grey70",
             hjust = "left") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha, fill = infections.colour) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha, fill = infections.colour) +
    geom_line(aes(y = q5), colour = infections.colour, alpha = line.alpha) +
    geom_line(aes(y = q95), colour = infections.colour, alpha = line.alpha) +
    scale_x_continuous("Day", breaks = my.day.breaks) + 
    scale_y_continuous("Daily New Infections") + 
    coord_cartesian(ylim=c(0, inf.upper.y)) +
    cowplot::theme_cowplot() +
    theme(text = element_text(size = 16))
  
  ggsave(filename = paste0(folder.name,"/",folder.name, "_infections.png"), plot = inf.plot, height = 4, width = 6, units="in", bg = "white")
  
  
  ## Incidence ----
  symp.plot <- plot.dat %>% 
    filter(var == "incidence") %>% 
    ggplot() + 
    aes(x = day) +
    geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
    annotate(geom = "text", x = vlines$value+1, y = 
               (plot.dat %>% filter(var == "incidence") %>% summarise(m = max(q95)) %>% pull()), 
             label = vlines$name, colour = "grey70",
             hjust = "left") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha, fill = symp.cases.colour) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha, fill = symp.cases.colour) +
    geom_line(aes(y = q5), colour = symp.cases.colour, alpha = line.alpha) +
    geom_line(aes(y = q95), colour = symp.cases.colour, alpha = line.alpha) +
    scale_x_continuous("Day", breaks = my.day.breaks) + 
    scale_y_continuous("Daily New Symptomatic Cases") + 
    cowplot::theme_cowplot() +
    theme(text = element_text(size = 16))
  
  ggsave(filename = paste0(folder.name,"/",folder.name, "_symptomatic_cases.png"), plot = symp.plot, height = 4, width = 6, units="in", bg = "white")
  
  ## Deaths ----
  death.plot <- plot.dat %>% 
    filter(var == "new.deaths") %>% 
    ggplot() + 
    aes(x = day) +
    geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
    annotate(geom = "text", x = vlines$value+1, y = 
               (plot.dat %>% filter(var == "new.deaths") %>% summarise(m = max(q95)) %>% pull()), 
             label = vlines$name, colour = "grey70",
             hjust = "left") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha, fill = death.colour) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha, fill = death.colour) +
    geom_line(aes(y = q5), colour = death.colour, alpha = line.alpha) +
    geom_line(aes(y = q95), colour = death.colour, alpha = line.alpha) +
    scale_x_continuous("Day", breaks = my.day.breaks) + 
    scale_y_continuous("Daily Deaths", limits = my_lims) +
    coord_cartesian(ylim = c(0, death.upper.y)) +
    cowplot::theme_cowplot() +
    theme()
  
  ggsave(filename = paste0(folder.name,"/",folder.name, "_deaths.png"), plot = death.plot, height = 4, width = 6, units="in", bg = "white")
  


  # Ward and ICU ccupancy plots
  occ.plot.dat <- df %>% 
    filter(var %in% c("total.beds.occupied", "ICU.occupancy")) %>% 
    pivot_wider(names_from = var, values_from = value) %>% 
    mutate(ward.occupancy = total.beds.occupied - ICU.occupancy) %>% 
    pivot_longer(cols = ICU.occupancy:ward.occupancy) %>% 
    group_by(name, day) %>% 
    summarise(x=list(enframe(
      quantile(value, probs=c(0.05, 
                              #0.1, 0.15, 0.2, 
                              0.25, 0.75, 
                              #0.8, 0.85, 0.9, 
                              0.95)),"quantiles", "value"))) %>%
    unnest(x) %>% mutate(quantiles = as.numeric(gsub(pattern = "%",replacement = "",x = quantiles))) %>%
    mutate(quantiles = paste0("q",quantiles)) %>% pivot_wider(names_from = quantiles, values_from=value)
  
  
  ward.plot <- occ.plot.dat %>% 
    filter(name == "ward.occupancy") %>% 
    ggplot() + 
    aes(x = day) +
    geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
    annotate(geom = "text", x = vlines$value+1, y = ward.upper.y*0.95,
             label = vlines$name, colour = "grey70",
             hjust = "left") +
    geom_hline(yintercept = ward.capacity, lty = 2, colour = "grey30") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha, fill = ward.colour) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha, fill = ward.colour) +
    geom_line(aes(y = q5), colour = ward.colour, alpha = line.alpha) +
    geom_line(aes(y = q95), colour = ward.colour, alpha = line.alpha) +
    scale_x_continuous("Day", breaks = my.day.breaks) + 
    scale_y_continuous("Occupied Ward Beds") + 
    coord_cartesian(ylim = c(0,ward.upper.y)) +
    cowplot::theme_cowplot() +
    theme()
  
  ggsave(filename = paste0(folder.name,"/",folder.name, "_ward_occupancy.png"), plot = ward.plot, height = 4, width = 6, units="in", bg = "white")
  
  
  icu.plot <- occ.plot.dat %>% 
    filter(name == "ICU.occupancy") %>% 
    ggplot() + 
    aes(x = day) +
    geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
    annotate(geom = "text", x = vlines$value+1, y = icu.upper.y*0.95,
             label = vlines$name, colour = "grey70",
             hjust = "left") +
    geom_hline(yintercept = icu.capacity, lty = 2, colour = "grey30") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha, fill = icu.colour) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha, fill = icu.colour) +
    geom_line(aes(y = q5), colour = icu.colour, alpha = line.alpha) +
    geom_line(aes(y = q95), colour = icu.colour, alpha = line.alpha) +
    scale_x_continuous("Day", breaks = my.day.breaks) + 
    scale_y_continuous("Occupied ICU Beds") + 
    coord_cartesian(ylim = c(0,icu.upper.y)) +
    cowplot::theme_cowplot() +
    theme()
  
  ggsave(filename = paste0(folder.name,"/",folder.name, "_icu_occupancy.png"), plot = icu.plot, height = 4, width = 6, units="in", bg = "white")
  
  
  ## Absenteeism ----
  # Assumes symptomatic individuals miss 10 days of work
  
  # Convert MATLAB array but retain age groups
  df.az1 <- convert.mat(dat$AZ1.data, "AZ1")
  df.az2 <- convert.mat(dat$AZ2.data, "AZ2")
  df.M1 <- convert.mat(dat$M1.data, "M1")
  df.M2 <- convert.mat(dat$M2.data, "M2")
  df.P1 <- convert.mat(dat$P1.data, "P1")
  df.P2 <- convert.mat(dat$P2.data, "P2")
  df.unvax <- convert.mat(dat$unvax.data, "unvax")

  df <- df.az1 %>% select(-group) %>%
    left_join(., df.az2 %>% select(-group), by = c("var","sim","day","age.groups")) %>%
    left_join(., df.M1 %>% select(-group), by = c("var","sim","day","age.groups")) %>%
    left_join(., df.M2 %>% select(-group), by = c("var","sim","day","age.groups")) %>%
    left_join(., df.P1 %>% select(-group), by = c("var","sim","day","age.groups")) %>%
    left_join(., df.P2 %>% select(-group), by = c("var","sim","day","age.groups")) %>%
    left_join(., df.unvax %>% select(-group), by = c("var","sim","day","age.groups")) %>%
    mutate(value = value.x + value.y + value.x.x + value.y.y + value.x.x.x + value.y.y.y + value) %>%
    select(-value.x, -value.y, -value.x.x, -value.y.y, -value.x.x.x, -value.y.y.y)
  
  
  df <- df %>%
    filter(var == "incidence") %>%
    left_join(., pc.workforce, by = "age.groups") %>%
    mutate(value = ceiling(value * pc)) %>%
    select(-pc) %>%
    group_by(var, sim, day) %>%
    summarise(value = sum(value)) %>%
    filter(day >= 0, day <= max.day)

  # Rescale symptomatic infections by percentage of each age group in the workforce, to subsequently count
  # number of days absent from work

  abs.dat <- df %>%
    group_by(sim) %>%
    mutate(absent = slider::slide_index_dbl(
      .x = value,
      .i = day,
      .f = sum,
      .before = 10,
      .after = 0)) %>%
    select(-value, -var) %>%
    group_by(day) %>%
    summarise(x=list(enframe(
      quantile(absent, probs=c(0.05,
                               #0.1, 0.15, 0.2,
                               0.25, 0.75,
                               #0.8, 0.85, 0.9,
                               0.95)),"quantiles", "value"))) %>%
    unnest(x) %>% mutate(quantiles = as.numeric(gsub(pattern = "%",replacement = "",x = quantiles))) %>%
    mutate(quantiles = paste0("q",quantiles)) %>% pivot_wider(names_from = quantiles, values_from=value)

  abs.plot <- abs.dat %>%
    ggplot() +
    aes(x = day) +
    geom_vline(data = vlines, aes(xintercept = value), colour = "grey70") +
    annotate(geom = "text", x = vlines$value+1, y =
               (abs.dat %>% summarise(m = max(q95)) %>% pull()),
             label = vlines$name, colour = "grey70",
             hjust = "left") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = inner.alpha, fill = symp.cases.colour) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = outer.alpha, fill = symp.cases.colour) +
    geom_line(aes(y = q5), colour = symp.cases.colour, alpha = line.alpha) +
    geom_line(aes(y = q95), colour = symp.cases.colour, alpha = line.alpha) +
    scale_x_continuous("Day", breaks = my.day.breaks) +
    scale_y_continuous("Absenteeism") +
    cowplot::theme_cowplot() +
    theme()

  ggsave(filename = paste0(folder.name,"/",folder.name, "_absenteeism.png"), plot = abs.plot, height = 4, width = 6, units="in", bg = "white")
}
