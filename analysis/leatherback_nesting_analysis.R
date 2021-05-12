#################################################################
# LEATHERBACK SEA TURTLE NESTING ANALYSES
#################################################################


#================================================================
# SETUP
#================================================================

library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(matrixStats)
library(rstan)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(ggridges)
library(here)
if(.Platform$OS.type == "windows") options(device = windows)
theme_set(theme_bw(base_size = 14))
bayesplot_theme_update(panel.grid = element_blank(), plot.title = element_text(size = 10))
color_scheme_set("gray")
source(here("analysis", "GeomFlatViolin.R"))
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
  rename(ppt = precip) %>% mutate(date = as_date(date, format = "%m/%d/%Y"))

# Change outlier dp_avg and p_avg on two dates to NA
weather <- weather_raw %>% mutate(year = year(date), doy = yday(date), .after = date) %>%
  mutate(dp_avg = replace(dp_avg, p_avg == 0, NA), p_avg = ifelse(p_avg == 0, NA, p_avg)) %>% 
  as.data.frame()

# Filter by window of encounter dates
# Average weather variables by year
weather_agg <- weather %>% group_by(year) %>% 
  filter(between(doy, min(nest_raw$doy_encounter, na.rm = TRUE), max(nest_raw$doy_encounter, na.rm = TRUE))) %>%
  summarize(across(humid_max:sst, mean), .groups = "keep") %>% 
  rename_with(~ gsub("avg", "yavg", .x), ends_with("avg")) %>% 
  rename(ppt_yavg = ppt, sst_yavg = sst) %>% select(year, ends_with("yavg")) %>% 
  as.data.frame()

# Remove surveys without an ID'd female (for now)
# Identify first encounter of each turtle each year she was sighted
# Merge weather into nest data
# Create centered and scaled predictors
nest <- nest_raw %>% filter(!is.na(ID)) %>% select(-notes) %>% group_by(name, year) %>% 
  mutate(first_of_year = rank(doy_encounter) == 1, .after = doy_encounter) %>% ungroup() %>% 
  mutate(year_ctr = scale(year, scale = FALSE), fyear = factor(year), .after = year) %>% 
  left_join(select(weather, c(date, ends_with("avg"), ppt, sst)), by = c("date_encounter" = "date")) %>% 
  left_join(weather_agg) %>% arrange(name, year, date_encounter) %>% 
  mutate(across(ends_with("yavg"), scale, .names = "{gsub('yavg', 'std', .col)}")) %>% 
  as.data.frame()


#================================================================
# MODELING
#
# Fit hierarchical models to different response variables
# using rstanarm
#================================================================

#----------------------------------------------------------------
# Breeding phenology: encounter date
#----------------------------------------------------------------

# turtle-level and year-level hierarchical intercepts
lmer_doy0 <- stan_lmer(doy_encounter ~ (1 | name) + (1 | fyear), data = nest, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy0)
summary(lmer_doy0, pars = "alpha", regex_pars = "igma", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy0, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# linear trend
lmer_doy1 <- stan_lmer(doy_encounter ~ year_ctr + (1 | name) + (1 | fyear), data = nest, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy1)
summary(lmer_doy1, pars = c("alpha","beta"), regex_pars = "igma", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy1, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# linear trend
# all weather variables (average annual anomalies)
lmer_doy2 <- stan_lmer(doy_encounter ~ year_ctr + humid_std + ws_std + p_std + dp_std + 
                         ppt_std + sst_std + (1 | name) + (1 | fyear), data = nest, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy2)
summary(lmer_doy2, pars = c("alpha","beta"), regex_pars = "igma", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy2, pars = "varying")

# # DOY of a female's *first* encounter in a given year
# # turtle-level and year-level hierarchical intercepts
# # linear change over time (year)
# # => very similar to lmer_doy but fewer obs/turtle so more uncertainty
# lmer_doy1st1 <- stan_lmer(doy_encounter ~ year_ctr + (1 | name) + (1 | fyear),
#                          data = nest, subset = first_of_year,
#                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
# print(lmer_doy1st1)


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

mod_name <- "lmer_doy2"
mod <- get(mod_name)
yrep <- posterior_predict(mod)
indx <- sample(nrow(yrep), 100)

# PPD marginal density
ppc_dens_overlay(mod$y, yrep[indx,]) + ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))
ggsave(filename=here("analysis", "results", paste0(mod_name, "_ppc_dens_overlay.png")),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

# PPD marginal density grouped by year (takes a while)
ppc_dens_overlay_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$year) +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPD scatterplots grouped by year
ppc_scatter_avg_grouped(mod$y, yrep, group = mod$glmod$fr$year) +
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# Normal QQ plot of year-level random effects
ranef(mod)$fyear %>% rename(intercept = `(Intercept)`) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# Normal QQ plot of turtle-level random effects
ranef(mod)$name %>% rename(intercept = `(Intercept)`) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))
 


#================================================================
# FIGURES
#================================================================

#----------------------------------------------------------------
# Exploratory plots of weather data
#----------------------------------------------------------------

# Annual climatologies of each variable
dev.new(width = 12, height = 5)

weather %>% group_by(date) %>% 
  summarize(across(c(ends_with("avg"), ppt, sst), mean, na.rm = TRUE), .groups = "drop") %>% 
  pivot_longer(-date, names_to = "variable") %>% 
  mutate(variable = recode(variable, humid_avg = "Humidity~('%')", ws_avg = "Wind~speed~(mph)",
                           p_avg = "Pressure~'(in Hg)'", t_avg = "Temperature~(degree * C)", 
                           dp_avg = "Dew~point~(degree * C)", ppt = "Precipitation~(cm)", 
                           sst = "SST~(degree * C)")) %>% 
  mutate(doy = yday(date), .after = date) %>% group_by(variable, doy) %>% 
  summarize(lb = quantile(value, 0.05, na.rm = TRUE), med = median(value, na.rm = TRUE), 
            ub = quantile(value, 0.95, na.rm = TRUE)) %>% 
  ggplot(aes(x = as_date(as_date(doy), format = "%m-%d"), y = med)) +
  annotate("rect", xmin = min(as_date(nest$doy_encounter), na.rm = TRUE),
           xmax = max(as_date(nest$doy_encounter), na.rm = TRUE),
           ymin = -Inf, ymax = Inf, fill = "steelblue4", alpha = 0.5) +
  geom_line(lwd = 0.7, col = "gray40") +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "gray40", alpha = 0.7) +
  scale_x_date(date_minor_breaks = "1 month", date_labels = "%b") + 
  xlab("Date") + ylab("Daily average") +
  facet_wrap(vars(variable), ncol = 4, scales = "free_y", labeller = label_parsed) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = NA))

ggsave(filename=here("analysis", "results", "climatology.png"),
       width=12, height=5, units="in", dpi=300, type="cairo-png")

#----------------------------------------------------------------
# Breeding phenology: encounter date
#----------------------------------------------------------------

# Time series of encounter DOY, all data
ppd <- as.data.frame(t(yrep[indx,])) %>% mutate(fyear = nest$fyear) %>% 
  pivot_longer(-fyear, names_to = "iter", values_to = "doy_yrep") %>% 
  mutate(date_encounter = as_date(doy_yrep), "%m-%d")

dev.new(width = 7, height = 5)

nest %>% 
  ggplot(aes(x = fyear, y = as_date(format(date_encounter, "%m-%d"), format = "%m-%d"))) +
  geom_flat_violin(aes(fill = alpha("steelblue4", 0.5)), scale = "width", width = 0.8, trim = FALSE, col = NA) +
  geom_flat_violin(data = ppd, aes(fill = "lightgray"), scale = "width", width = -0.8, col = NA) +
  geom_jitter(aes(group = 1), width = 0.2, pch = 16, size = 1, col = "steelblue4", alpha = 0.5) + 
  scale_y_date(date_breaks = "1 month", date_labels = "%b", expand = expansion(add = -30)) +
  scale_fill_identity(name = "", guide = "legend", labels = c("observed", "predicted")) +
  xlab("Year") + ylab("Encounter DOY") + 
  theme(legend.position = "top", legend.box.margin = margin(b = -12), panel.grid = element_blank())

ggsave(filename=here("analysis", "results", "doy_all_encounters.png"),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

# Time series of encounter DOY, female averages
dev.new(width = 7, height = 5)

nest %>% group_by(year, name) %>% summarize(mean_date_encounter = mean(date_encounter), .groups = "drop") %>% 
  ggplot(aes(x = year, y = as_date(format(mean_date_encounter, "%m-%d"), format = "%m-%d"), group = name)) +
  geom_line(col = "steelblue4", alpha = 0.5) + scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Encounter DOY") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

ggsave(filename=here("analysis", "results", "doy_female_avg.png"),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

# Joyplot of predicted DOY for each female in an average year
# Only include females encountered in >= 3 years (relatively informative for intercept)
mod_name <- "lmer_doy2"
mod <- get(mod_name)

dat <- nest %>% group_by(name) %>% 
  select(name, year_ctr, fyear, doy_encounter, humid_std:sst_std) %>% 
  summarize(n_year = n_distinct(year_ctr), year_ctr = 0, fyear = fyear[1], doy_encounter = 0, 
            across(ends_with("std"), function(x) 0), .groups = "keep") %>% 
  filter(n_year > 3) %>% as.data.frame()

pred <- posterior_linpred(mod, newdata = dat, re.form = ~ (1 | name))

dev.new(height = 10, width = 7)

cbind(dat, t(pred)) %>% 
  pivot_longer(-(name:sst_std), names_to = "iter", values_to = "doy_pred") %>% 
  mutate(name = reorder(name, doy_pred, mean)) %>% 
  ggplot(aes(x = as_date(as_date(doy_pred), format = "%m-%d"), 
             y = name, height = stat(density))) +
  geom_density_ridges(col = alpha("steelblue4", 0.7), fill = alpha("steelblue4", 0.4), 
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_x_date(date_breaks = "2 weeks", date_minor_breaks = "1 week", date_labels = "%b %d",
               limits = function(x) c(ceiling_date(x[1], "week"), x[2] - 3)) +
xlab("Encounter DOY") + 
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11, margin = margin(r = -10)),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.line.x = element_line(size = 0.5), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename=here("analysis", "results", "doy_female_joyplot.png"),
       width=7*0.9, height=10*0.9, units="in", dpi=300, type="cairo-png")


#----------------------------------------------------------------
# Body size
#----------------------------------------------------------------

# Time series of body size, all data
nest %>% 
  ggplot(aes(x = year, y = ccl_max)) + 
  geom_jitter(width = 0.2, pch = 16, alpha = 0.5) + 
  scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") + theme_bw(base_size = 14) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

# Time series of body size, female averages
nest %>% group_by(year, name) %>% summarize(mean_ccl_max = mean(ccl_max), .groups = "drop") %>% 
  ggplot(aes(x = year, y = mean_ccl_max, group = name)) +
  geom_line(alpha = 0.5) + scale_x_continuous(breaks = sort(unique(nest$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())





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





