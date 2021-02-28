#################################################################
# LEATHERBACK SEA TURTLE NESTING ANALYSES
#################################################################

### NOTES
# (1) One turtle (ID 1587, encountered only once) has name "2-May". Error?
# (2) One encounter (664, HAEDI) has date 1/5/1900


#================================================================
# SETUP
#================================================================

if(.Platform$OS.type == "windows") options(device = windows)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
library(dplyr)
library(readxl)
library(lubridate)
library(matrixStats)
library(rstan)
library(shinystan)
library(here)

#================================================================
# DATA
#================================================================

# Import raw nest data and rename variables
# Remove (for now) one encounter whose date is recorded in Excel as 1/4/1900
nest_raw <- read_excel(here("data", "Historical Leatherback Data_updated1.18.2020.xlsx"),
                       sheet = "2007-2020") %>% as.data.frame() %>% 
  rename(name = `Turtle Name`, ID = `Turtle ID`, date_1st = `Date of 1st Encounter`,
         year = Year, date_encounter = `Date of Encounter`, encounter = `Encounter ID`,
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
  mutate(across(starts_with("date"), as_date)) %>% 
  mutate(doy_encounter = yday(date_encounter), .after = date_encounter) %>% 
  filter(encounter != 664) %>% 
  arrange(name, year, date_encounter)

# Remove surveys without an ID'd female (for now)
nest <- filter(nest_raw, !is.na(ID)) %>% select(-notes)


#================================================================
# MODELING
#
# Fit multivariate random walk state-space (MRWSS) models
# to different response variables
#================================================================

# DOY of encounter
mrwss_doy <- stan(file = here("analysis","MRWSS.stan"),
                  data = list(N = nrow(nest), year = as.numeric(factor(nest$year)),
                              turtle = as.numeric(factor(nest$name)), 
                              y = nest$doy_encounter),
                  pars = c("mu","sigma_omega","sigma_nu","tau","theta","x"),
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)

print(mrwss_doy, pars = c("theta","x"), include = FALSE, probs = c(0.025, 0.5, 0.975))





#================================================================
# FIGURES
#================================================================

#----------------------------------------------------------------
# Exploratory plots
#----------------------------------------------------------------

# Time series of encounter DOY, all data
nest %>% ggplot(aes(x = year, y = as_date(format(date_encounter, "%m-%d"), format = "%m-%d"))) +
  geom_jitter(width = 0.2, pch = 16, alpha = 0.5) + 
  scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("DOY") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor.x = element_blank())

# Time series of encounter DOY, female averages
nest %>% group_by(year, name) %>% summarize(mean_date_encounter = mean(date_encounter), .groups = "drop") %>% 
  ggplot(aes(x = year, y = as_date(format(mean_date_encounter, "%m-%d"), format = "%m-%d"), group = name)) +
  geom_line(alpha = 0.5) + scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("DOY") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor.x = element_blank())

# Time series of body size, all data
nest %>% ggplot(aes(x = year, y = CCL_max)) + 
  geom_jitter(width = 0.2, pch = 16, alpha = 0.5) + 
  scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor.x = element_blank())

# Time series of body size, female averages
nest %>% group_by(year, name) %>% summarize(mean_CCL_max = mean(CCL_max), .groups = "drop") %>% 
  ggplot(aes(x = year, y = mean_CCL_max, group = name)) +
  geom_line(alpha = 0.5) + scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor.x = element_blank())





