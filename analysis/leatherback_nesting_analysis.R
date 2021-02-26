#################################################################
# LEATHERBACK SEA TURTLE NESTING ANALYSES
#################################################################

### NOTES
# (1) One turtle (ID 1587, encountered only once) has name "2-May". Error?


#================================================================
# SETUP
#================================================================

if(.Platform$OS.type == "windows") options(device = windows)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
library(dplyr)
library(readxl)
library(matrixStats)
library(rstan)
library(shinystan)
library(here)

#================================================================
# DATA
#================================================================

# Import raw nest data and rename variables
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
  arrange(name, year, date_encounter)

# Remove surveys without an ID'd female (for now)
nest <- filter(nest_raw, !is.na(ID)) %>% select(-notes)





