#################################################################
# LEATHERBACK SEA TURTLE NESTING ANALYSES
#################################################################

### NOTES
# (1) One turtle (ID 1587, encountered only once) has name "2-May". Error?
# (2) One encounter (664, HAEDI) has date 1/5/1900


#================================================================
# SETUP
#================================================================

library(dplyr)
library(readxl)
library(lubridate)
library(matrixStats)
library(rstan)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(here)
if(.Platform$OS.type == "windows") options(device = windows)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

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
  mutate(date = as_date(date, format = "%m/%d/%Y"))

# Filter by window of encounter dates
weather <- weather_raw %>% mutate(year = year(date), doy = yday(date), .after = date) %>% 
  filter(doy >= min(nest_raw$doy_encounter, na.rm = TRUE) & 
           doy <= max(nest_raw$doy_encounter, na.rm = TRUE)) %>% as.data.frame()

# Average weather variables by year
weather_agg <- weather %>% group_by(year) %>% 
  summarize(across(humid_max:sst, mean), .groups = "keep") %>% 
  rename_with(~ gsub("avg", "yavg", .x), ends_with("avg")) %>% 
  select(year, ends_with("yavg")) %>% as.data.frame()

# Remove surveys without an ID'd female (for now)
# Identify first encounter of each turtle each year she was sighted
# Shift year so origin is first year in data
# Merge weather into nest data
nest <- nest_raw %>% filter(!is.na(ID)) %>% select(-notes) %>% group_by(name, year) %>% 
  mutate(first_of_year = rank(doy_encounter) == 1, .after = doy_encounter) %>% ungroup() %>% 
  mutate(year0 = year - min(year)) %>% arrange(name, year, date_encounter) %>% 
  left_join(select(weather, c(date, ends_with("avg"))), by = c("date_encounter" = "date")) %>% 
  left_join(weather_agg) %>% as.data.frame()


#================================================================
# MODELING
#
# Fit hierarchical models to different response variables
# using rstanarm
#================================================================

#----------------------------------------------------------------
# Breeding phenology: encounter date
#----------------------------------------------------------------

# DOY of encounter
# turtle-level and year-level hierarchical intercepts
# linear change over time (year)
lmer_doy <- stan_lmer(doy_encounter ~ year0 + (1 | year0) + (1 | name), data = nest, 
                      chains = getOption("mc.cores"), iter = 5000, warmup = 1000)

print(lmer_doy)
summary(lmer_doy)





# DOY of *first* encounter
# turtle-level and year-level hierarchical intercepts
# linear change over time (year)
# => very similar to lmer_doy but fewer obs/turtle so more uncertainty
lmer_doy1st <- stan_lmer(doy_encounter ~ year0 + (1 | year0) + (1 | name), 
                      data = nest, subset = first_of_year,
                      chains = getOption("mc.cores"), iter = 5000, warmup = 1000)

print(lmer_doy1st)
summary(lmer_doy1st)



#----------------------------------------------------------------
# Body size
#----------------------------------------------------------------

# Maximum curved carapace length
# turtle-level and year-level hierarchical intercepts
# linear change over time (year)
lmer_ccl <- stan_lmer(ccl_max ~ year + (1 | year) + (1 | name), data = nest, 
                      chains = getOption("mc.cores"), iter = 5000, warmup = 1000)

print(lmer_ccl)
summary(lmer_ccl)



#----------------------------------------------------------------
# Diagnostic plots 
#----------------------------------------------------------------

mod <- lmer_doy1st
yrep <- posterior_predict(mod)
indx <- sample(nrow(yrep), 500)

# PPD marginal density
ppc_dens_overlay(mod$y, yrep[indx,]) + ggtitle(deparse(mod$glmod$formula))

# PPD marginal density grouped by year (takes a while)
ppc_dens_overlay_grouped(mod$y, yrep[indx[1:100],], group = mod$glmod$fr$year) +
  ggtitle(deparse(mod$glmod$formula))

# PPD scatterplots grouped by year
ppc_scatter_avg_grouped(mod$y, yrep, group = mod$glmod$fr$year) +
  geom_abline(intercept = 0, slope = 1) + ggtitle(deparse(mod$glmod$formula))

# Normal QQ plot of year-level random effects
ranef(mod)$year %>% rename(intercept = `(Intercept)`) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + ggtitle(deparse(mod$glmod$formula))

# Normal QQ plot of turtle-level random effects
ranef(mod)$name %>% rename(intercept = `(Intercept)`) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + ggtitle(deparse(mod$glmod$formula))



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
nest %>% ggplot(aes(x = year, y = ccl_max)) + 
  geom_jitter(width = 0.2, pch = 16, alpha = 0.5) + 
  scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor.x = element_blank())

# Time series of body size, female averages
nest %>% group_by(year, name) %>% summarize(mean_ccl_max = mean(ccl_max), .groups = "drop") %>% 
  ggplot(aes(x = year, y = mean_ccl_max, group = name)) +
  geom_line(alpha = 0.5) + scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor.x = element_blank())





################
## DEPRECATED ##
################

#================================================================
# MODELING
#
# Fit hierarchical random walk state-space (HRWSS) models
# to different response variables
#================================================================

# DOY of encounter
hrwss_doy <- stan(file = here("analysis","HRWSS.stan"),
                  data = list(N = nrow(nest), year = as.numeric(factor(nest$year)),
                              turtle = as.numeric(factor(nest$name)), 
                              y = nest$doy_encounter),
                  pars = c("mu","sigma_alpha","sigma_theta","sigma","alpha","theta","y_hat","LL"),
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000,
                control = list(adapt_delta = 0.95))

print(hrwss_doy, pars = c("theta","y_hat","LL"), include = FALSE, probs = c(0.025, 0.5, 0.975))


# Max curved carapace length
hrwss_ccl <- stan(file = here("analysis","HRWSS.stan"),
                  data = list(N = nrow(ccl), year = as.numeric(factor(ccl$year)),
                              turtle = as.numeric(factor(ccl$name)), y = ccl$CCL_max),
                  pars = c("mu","sigma_alpha","sigma_theta","sigma","alpha","theta","y_hat","LL"),
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000,
                  control = list(adapt_delta = 0.95))

print(hrwss_ccl, pars = c("theta","y_hat","LL"), include = FALSE, probs = c(0.025, 0.5, 0.975))





