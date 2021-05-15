#################################################################
# LEATHERBACK SEA TURTLE NESTING ANALYSES
#################################################################


#================================================================
# SETUP
#================================================================

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

# Load and manipulate data
source(here("analysis", "01_leatherback_nesting_data.R"))

# Load saved stanreg and loo objects
if(file.exists(here("analysis", "results", "leatherback_nesting_models.RData")))
  load(here("analysis", "results", "leatherback_nesting_models.RData"))


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
lmer_doy0 <- stan_lmer(doy_encounter ~ (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy0)
summary(lmer_doy0, pars = "alpha", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy0, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
lmer_doy1 <- stan_lmer(doy_encounter ~ year_ctr + (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy1)
summary(lmer_doy1, pars = c("alpha","beta"), probs = c(0.025, 0.5, 0.975))
summary(lmer_doy1, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
# all weather variables (average annual anomalies)
# except dew point (correlation with temperature is 0.91)
lmer_doy2 <- stan_lmer(doy_encounter ~ year_ctr + humid_std + ws_std + p_std + t_std + 
                         ppt_std + sst_std + (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy2)
summary(lmer_doy2, pars = c("alpha","beta"), probs = c(0.025, 0.5, 0.975))
summary(lmer_doy2, pars = "varying", regex_pars = "igma")


#----------------------------------------------------------------
# Somatic growth: 
# change in max curved carapace length / year 
# for turtles encountered multiple times
#----------------------------------------------------------------

# turtle-level and year-level hierarchical intercepts
lmer_sgr0 <- stan_lmer(sgr ~ (1 | name) + (1 | fyear), data = size, 
                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr0, 3)
summary(lmer_sgr0, pars = "alpha", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr0, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# initial size effect
lmer_sgr1 <- stan_lmer(sgr ~ ccl_max0_std + (1 | name) + (1 | fyear), data = size, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr1, 3)
summary(lmer_sgr1, pars = "alpha", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr1, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# initial size effect with turtle-varying slopes 
# (diagonal random effects covariance matrix)
lmer_sgr2 <- stan_lmer(sgr ~ ccl_max0_std + (ccl_max0_std || name) + (1 | fyear), data = size, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr2, 3)
summary(lmer_sgr2, pars = "alpha", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr2, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# initial size effect with turtle-varying slopes 
# linear time trend
lmer_sgr3 <- stan_lmer(sgr ~ ccl_max0_std + year_ctr + (ccl_max0_std || name) + (1 | fyear), data = size, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr3, 3)
summary(lmer_sgr3, pars = "alpha", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr3, pars = "varying", regex_pars = "igma")


#----------------------------------------------------------------
# Proportion of neophytes
#----------------------------------------------------------------

# intercept only
glm_neo0 <- stan_glm(cbind(neophyte, remigrant) ~ 1, 
                     data = neophyte, family = binomial, 
                     chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glm_neo0, 2)
summary(glm_neo0, probs = c(0.025, 0.5, 0.975), digits = 2)

# linear trend
glm_neo1 <- stan_glm(cbind(neophyte, remigrant) ~ year_ctr, 
                     data = neophyte, family = binomial, 
                     chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glm_neo1, 2)
summary(glm_neo1, probs = c(0.025, 0.5, 0.975), digits = 2)

## Model selection using loo
loo_neo <- lapply(stanreg_list(glm_neo0, glm_neo1), loo)
loo_neo
print(loo_compare(loo_neo)[,], 3)


#----------------------------------------------------------------
# Hatchling emergence success
# emerged = hatched - live in nest - dead in nest
#----------------------------------------------------------------

#---------------#
# All nest data #
#---------------#

# year-level hierarchical intercept
glmer_anest0 <- stan_glmer(emergence_rate ~ (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_anest0, 2)
summary(glmer_anest0, pars = "alpha", probs = c(0.025, 0.5, 0.975))
summary(glmer_anest0, pars = "varying", regex_pars = "igma")

# year-level hierarchical intercept
# linear trend
glmer_anest1 <- stan_glmer(emergence_rate ~ year_ctr + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_anest1, 2)
summary(glmer_anest1, pars = c("alpha", "beta"), probs = c(0.025, 0.5, 0.975), digits = 3)
summary(glmer_anest1, pars = "varying", regex_pars = "igma")

# year-level hierarchical intercept
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
glmer_anest2 <- stan_glmer(emergence_rate ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                             (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_anest2, 3)
summary(glmer_anest2, pars = c("alpha", "beta"), probs = c(0.025, 0.5, 0.975), digits = 3)
summary(glmer_anest2, pars = "varying", regex_pars = "igma")

#--------------------------------------------#
# Only nests where a female was encountered  #
#--------------------------------------------#

# turtle-level and year-level hierarchical intercepts
glmer_enest0 <- stan_glmer(emergence_rate ~ (1 | name) + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch, 
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_enest0, 2)
summary(glmer_enest0, pars = "alpha", probs = c(0.025, 0.5, 0.975))
summary(glmer_enest0, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
glmer_enest1 <- stan_glmer(emergence_rate ~ year_ctr + (1 | name) + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_enest1, 2)
summary(glmer_enest1, pars = c("alpha", "beta"), probs = c(0.025, 0.5, 0.975), digits = 3)
summary(glmer_enest1, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
glmer_enest2 <- stan_glmer(emergence_rate ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                             (1 | name) + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_enest2, 3)
summary(glmer_enest2, pars = c("alpha", "beta"), probs = c(0.025, 0.5, 0.975), digits = 3)
summary(glmer_enest2, pars = "varying", regex_pars = "igma")

#-----------------------------------------------------------#
# All nests - unknown females drawn from hyperdistribution  #
#-----------------------------------------------------------#

# turtle-level and year-level hierarchical intercepts
glmer_hnest0 <- stan_glmer(emergence_rate ~ (1 | hname) + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch, 
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_hnest0, 2)
summary(glmer_hnest0, pars = "alpha", probs = c(0.025, 0.5, 0.975))
summary(glmer_hnest0, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
glmer_hnest1 <- stan_glmer(emergence_rate ~ year_ctr + (1 | hname) + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_hnest1, 2)
summary(glmer_hnest1, pars = c("alpha", "beta"), probs = c(0.025, 0.5, 0.975), digits = 3)
summary(glmer_hnest1, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
glmer_hnest2 <- stan_glmer(emergence_rate ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                             (1 | hname) + (1 | fyear), 
                           data = nest, family = binomial, weights = clutch,
                           chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glmer_hnest2, 3)
summary(glmer_hnest2, pars = c("alpha", "beta"), probs = c(0.025, 0.5, 0.975), digits = 3)
summary(glmer_hnest2, pars = "varying", regex_pars = "igma")



#================================================================
# Diagnostic plots 
#================================================================

mod_name <- "glmer_hnest2"
mod <- get(mod_name)
yrep <- posterior_predict(mod)
indx <- sample(nrow(yrep), 100)

#----------------------------------------------------------------
# Continuous responses
#----------------------------------------------------------------

# PPD marginal density
ppc_dens_overlay(mod$y, yrep[indx,]) + ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))
ggsave(filename=here("analysis", "results", paste0(mod_name, "_ppc_dens_overlay.png")),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

# PPD marginal density grouped by year (takes a while)
ppc_dens_overlay_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$year) +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPD scatterplots grouped by year
ppc_scatter_avg_grouped(mod$y, yrep, group = mod$glmod$fr$fyear) +
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

#----------------------------------------------------------------
# Binomial responses
#----------------------------------------------------------------

# PPD marginal histogram for binomial responses
ppc_rootogram(mod$y[,1], yrep[1:5,]) +   
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPD marginal histogram for binomial responses
ppc_hist(mod$y[,1]/rowSums(mod$y), sweep(yrep[1:5,], 2, rowSums(mod$y), "/"), binwidth = 0.05) + 
  scale_x_continuous(limits = c(0,1)) +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPD empirical CDF overlay for binomial responses
ppc_ecdf_overlay(mod$y[,1]/rowSums(mod$y), sweep(yrep[indx,], 2, rowSums(mod$y), "/")) +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPC for proportion zeros
ppc_stat(mod$y[,1], yrep[indx,], stat = function(x) mean(x == 0)) +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

#----------------------------------------------------------------
# All responses
#----------------------------------------------------------------

# Normal QQ plot of year-level random effects
ranef(mod)$fyear %>% rename(intercept = `(Intercept)`) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# Normal QQ plot of turtle-level random effects
ranef(mod)$name %>% rename(intercept = `(Intercept)`) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), plot.title = element_text(size = 11)) +
  xlab("Normal quantiles") + ylab("Quantiles of turtle-specific intercepts") +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))
ggsave(filename=here("analysis", "results", paste0(mod_name, "_turtle-intercept_qqnorm.png")),
       width=7, height=7, units="in", dpi=300, type="cairo-png")


#================================================================
# Save stanreg and loo objects 
#================================================================

save(list = ls()[sapply(mget(ls()), function(x)
  any(class(x) == "stanreg") | any(sapply(x, class) == "loo"))], 
  file = here("analysis","results","leatherback_nesting_models.RData"))


