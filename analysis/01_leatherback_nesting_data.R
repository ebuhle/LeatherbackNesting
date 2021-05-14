#################################################################
# LEATHERBACK SEA TURTLE NESTING DATA
#################################################################


#================================================================
# SETUP
#================================================================

library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(zoo)
library(matrixStats)
library(here)


#================================================================
# DATA
#================================================================

# Import raw nest data and rename variables
nest_raw <- read_excel(here("data", "Historical Leatherback Data_updated1.18.2020.xlsx"),
                       sheet = "2007-2020") %>% as.data.frame() %>% 
  rename(name = `Turtle Name`, ID = `Turtle ID`, date_first = `Date of 1st Encounter`,
         year = Year, date_encounter = `Date of Encounter`, encounter = `Encounter ID`,
         ccl_min = CCL_min, ccl_max = CCL_max, ccw = CCW,
         fat_depth = `Fat Depth`, time_encounter = `Time Spotted Turtle`, 
         behav_encounter = `Behavior when first spotted`, crawl = `Crawl ID_am`,
         date_survey = `Survey Date_am`, beach = Beach, zone = Zone, crawl_type = `Crawl Type`,
         nest = `Nest ID`, dist_hwl = `Distance to High Water Line`,
         dist_dune = `Distance to Toe of Dune`, lat = Latitude, lon = Longitude,
         date_emergence = `Date of Hatchling Emergence`, incubation = `Incubation Time`,
         hatched = `# Hatched Eggshells`, unhatched = `# Uhatched/Whole Eggs`,
         live = `# Live in Nest`, dead = `# Dead in Nest`, live_pipped = `# Live Pipped`,
         dead_pipped = `# Dead Pipped`, clutch = `Total Clutch Size`, 
         hatch_rate = `Hatch Success`, emergence_rate = `Emergence Success`, 
         fate = `Fate Code`, notes = Notes) %>%
  mutate(doy_encounter = yday(date_encounter), .after = date_encounter)

# Import weather data
weather_raw <- read.csv(here("data","juno_weather.csv"), header = TRUE) %>% select(-1) %>% 
  rename(ppt = precip) %>% mutate(date = as_date(date, format = "%m/%d/%Y")) %>% 
  arrange(date)

# Interpolate outlier dp_avg and p_avg on two dates
# Ignore 26 NAs in SST (outside nesting season)
weather <- weather_raw %>% mutate(year = year(date), doy = yday(date), .after = date) %>%
  mutate(dp_avg = replace(dp_avg, p_avg == 0 & row_number() != n(), NA), # avoid trouble b/c
         p_avg = ifelse(p_avg == 0 & row_number() != n(), NA, p_avg),    # of terminal NAs
         dp_avg = c(zoo::na.approx(dp_avg), NA), p_avg = c(zoo::na.approx(p_avg), NA)) %>%
  as.data.frame()

# Filter by window of encounter dates
# Average weather variables by year
weather_agg <- weather %>% group_by(year) %>% 
  filter(between(doy, min(nest_raw$doy_encounter, na.rm = TRUE), max(nest_raw$doy_encounter, na.rm = TRUE))) %>%
  summarize(across(humid_max:sst, mean), .groups = "keep") %>% 
  rename_with(~ gsub("avg", "yavg", .x), ends_with("avg")) %>% 
  rename(ppt_yavg = ppt, sst_yavg = sst) %>% select(year, ends_with("yavg")) %>% 
  as.data.frame()

# Encounter data
# Identify first encounter of each turtle each year she was sighted
# Merge weather into encounter data
# Create centered and scaled predictors
turtle <- nest_raw %>% filter(!is.na(ID)) %>% select(-notes) %>% group_by(name, year) %>% 
  mutate(first_of_year = rank(doy_encounter) == 1, .after = doy_encounter) %>% ungroup() %>% 
  mutate(year_ctr = scale(year, scale = FALSE), fyear = factor(year), .after = year) %>% 
  left_join(select(weather, c(date, ends_with("avg"), ppt, sst)), by = c("date_encounter" = "date")) %>% 
  left_join(weather_agg) %>% arrange(name, year, date_encounter) %>% 
  mutate(across(ends_with("yavg"), scale, .names = "{gsub('yavg', 'std', .col)}")) %>% 
  as.data.frame()

# Subset of encounter data that recorded max CCL
# Average measurements within a year (variation is miniscule)
# Compute specific growth rate: SGR = 100 * log(final / initial) / t
# Turtle was only encountered once IFF is.na(size$sgr)
# Center and scale initial max CCL based on non-missing growth rates
size <- turtle %>% filter(!is.na(ccl_max)) %>% group_by(name, ID, year) %>% 
  summarize(fyear = unique(fyear), ccl_max = mean(ccl_max), 
            neophyte = ifelse(any(date_encounter == date_first), "neophyte", "remigrant")) %>% 
  group_by(name) %>% arrange(year) %>% 
  mutate(dyear = c(NA, diff(year)), ccl_max0 = lag(ccl_max),
         lgr = (ccl_max - ccl_max0) / dyear,
         sgr = 100 * (log(ccl_max) - log(ccl_max0)) / dyear) %>% 
  ungroup() %>% 
  mutate(ccl_max0_std = (ccl_max0 - mean(ccl_max0[!is.na(sgr)])) / sd(ccl_max0[!is.na(sgr)]),
         year_ctr = scale(year, scale = FALSE)) %>% 
  select(name, ID, year, year_ctr, fyear, neophyte, dyear, 
         ccl_max, ccl_max0, ccl_max0_std, lgr, sgr) %>% 
  as.data.frame()

# Frequencies of neophytes vs. repeat breeders
neophyte <- turtle %>% group_by(name, year) %>% 
  summarize(neophyte = ifelse(any(date_encounter == date_first), "neophyte", "remigrant")) %>% 
  group_by(year, neophyte) %>% summarize(count = n()) %>% ungroup() %>% 
  pivot_wider(id_cols = year, names_from = neophyte, values_from = count, values_fill = 0) %>% 
  mutate(fyear = factor(year), year_ctr = scale(year, scale = FALSE), .after = year) %>% 
  mutate(count = neophyte + remigrant) %>% 
  cbind(with(., Hmisc::binconf(x = neophyte, n = count))) %>%
  rename(p_neophyte = PointEst) %>% as.data.frame(row.names = 1:nrow(.))
  

  

