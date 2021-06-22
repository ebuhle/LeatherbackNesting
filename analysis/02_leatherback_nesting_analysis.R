#################################################################
# LEATHERBACK SEA TURTLE NESTING ANALYSES
#################################################################

## TO-DO
##
## * In addition to analyzing encounter DOY, what about nest DOY?
## * Model daily encounter or nest counts as observations
## * Include neophyte vs. remigrant as predictor in phenology models
## * Zone random effects (which models? all?)
## * Growth rate: linear instead of SGR
## * Emergence rate w/ encounter: include clutch sequence as predictor?

## Qs about survey protocol
##
## * Surveys are not daily; is that by design? What is the frequency?
## * What is done with nests constructed in between survey dates? 
##   Are they ignored? Added to the next survey's total?
## * Brost et al. 2015 note potential bias from excluding predated or washed-out
##   nests from nest success estimates. Are those data (coded as NA) available?
## * Brost et al. 2015 also seem to imply that post-emergence estimates of clutch size 
##   may be unreliable. Any reason for concern?
## * There are a few lat measurements that appear to be outside the recorded zone
##   plot(lat ~ zone, data = nest_raw)

#================================================================
# SETUP
#================================================================

## @knitr setup
library(rstanarm)
library(brms)
library(shinystan)
library(ggplot2)
library(patchwork)
library(bayesplot)
library(ggridges)
library(here)
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

## @knitr
if(.Platform$OS.type == "windows") options(device = windows)


#================================================================
# MODELING
#
# Fit hierarchical models to different response variables
# using rstanarm
#================================================================

#----------------------------------------------------------------
# Breeding phenology: encounter date
#----------------------------------------------------------------

## Linear mixed-effects models with rstanarm ##

# turtle-level and year-level hierarchical intercepts
## @knitr lmer_doy0-fit
lmer_doy0 <- stan_lmer(doy_encounter ~ (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr lmer_doy0-print
print(lmer_doy0)
## @knitr
summary(lmer_doy0, pars = "alpha", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy0, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
## @knitr lmer_doy1-fit
lmer_doy1 <- stan_lmer(doy_encounter ~ year_ctr + (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr lmer_doy1-print
print(lmer_doy1)
## @knitr
summary(lmer_doy1, pars = c("alpha","beta"), probs = c(0.025, 0.5, 0.975))
summary(lmer_doy1, pars = "varying", regex_pars = "igma")

# turtle-level and year-level hierarchical intercepts
# linear trend
# all weather variables (average annual anomalies)
# except dew point (correlation with temperature is 0.91)
## @knitr lmer_doy2-fit
lmer_doy2 <- stan_lmer(doy_encounter ~ year_ctr + humid_std + ws_std + p_std + t_std + 
                         ppt_std + sst_std + (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr lmer_doy2-print
print(lmer_doy2)
## @knitr
summary(lmer_doy2, pars = c("alpha","beta"), probs = c(0.025, 0.5, 0.975))
summary(lmer_doy2, pars = "varying", regex_pars = "igma")

## Hierarchical skew-normal models with brms ##

# turtle-level and year-level hierarchical intercepts
## @knitr sn_doy0-fit
sn_doy0 <- brm(doy_encounter ~ (1 | name) + (1 | fyear), 
               data = turtle, family = skew_normal,
               chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr sn_doy0-summary
summary(sn_doy0)
## @knitr

# turtle-level and year-level hierarchical intercepts
# year-level dispersion (sigma) and skew (alpha)
## @knitr sn_doy1-fit
sn_doy1 <- brm(bf(doy_encounter ~ (1 | name) + (1 | fyear), 
                  alpha ~ (1 | fyear)),
               data = turtle, family = skew_normal,
               chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr sn_doy1-summary
summary(sn_doy1)
## @knitr

## Model selection using loo
loo_doy <- lapply(list(lmer_doy0 = lmer_doy0, lmer_doy1 = lmer_doy1, lmer_doy2 = lmer_doy2, 
                       sn_doy0 = sn_doy0, sn_doy1 = sn_doy1), loo)
loo_doy
## @knitr doy-loo_compare
print(loo_compare(loo_doy)[,], 3)
## @knitr


#----------------------------------------------------------------
# Somatic growth: change in max curved carapace length 
#                 for turtles encountered multiple times
# Hierarchical von Bertalanffy models
#----------------------------------------------------------------

# turtle-level and year-level random K
# global Linf
## @knitr brm_vb0-fit
brm_vb0 <- brm(bf(ccl_max ~ ccl_max0 + (exp(logLinf) - ccl_max0) * (1 - inv_logit(logitK)^dyear), 
                  logitK ~ (1 | name) + (1 | fyear), logLinf ~ 1, nl = TRUE),
               data = size, 
               prior = c(prior(logistic(0,1), nlpar = "logitK"), 
                         prior(normal(5,5), nlpar = "logLinf")),
               chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr brm_vb0-summary
summary(brm_vb0)
## @knitr


#----------------------------------------------------------------
# Proportion of neophytes
# Binomial GLMs
#----------------------------------------------------------------

# intercept only
## @knitr glm_neo0-fit
glm_neo0 <- stan_glm(cbind(neophyte, remigrant) ~ 1, 
                     data = neophyte, family = binomial, 
                     chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr glm_neo0-print
print(glm_neo0, 2)
## @knitr
summary(glm_neo0, probs = c(0.025, 0.5, 0.975), digits = 2)

# linear trend
## @knitr glm_neo1-fit
glm_neo1 <- stan_glm(cbind(neophyte, remigrant) ~ year_ctr, 
                     data = neophyte, family = binomial, 
                     chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr glm_neo1-print
print(glm_neo1, 2)
## @knitr
summary(glm_neo1, probs = c(0.025, 0.5, 0.975), digits = 2)

## Model selection using loo
loo_neo <- lapply(stanreg_list(glm_neo0, glm_neo1), loo)
loo_neo
## @knitr neo-loo_compare
print(loo_compare(loo_neo)[,], 3)
## @knitr


#---------------------------------------------------------------
# Hatchling emergence success
# emerged = hatched - live in nest - dead in nest 
# Zero-inflated binomial GLMMs
#---------------------------------------------------------------

#---------------#
# All nest data #
#---------------#

# Binomial:
# zone-level hierarchical intercept
# year-level hierarchical intercept
# nest- (observation-) level overdispersion residual
## @knitr zib_anest0-fit
zib_anest0 <- brm(emerged | trials(clutch) ~ (1 | zone) + (1 | fyear) + (1 | nestID),
                  data = nest, family = zero_inflated_binomial(), 
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr zib_anest0-summary
summary(zib_anest0)
## @knitr

# Binomial:
# year-level hierarchical intercept
# nest- (observation-) level overdispersion residual
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
## @knitr zib_anest1-fit
zib_anest1 <- brm(emerged | trials(clutch) ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                    (1 | fyear) + (1 | nestID),
                  data = nest, family = zero_inflated_binomial(), 
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr zib_anest1-summary
summary(zib_anest1)
## @knitr

# Binomial:
# year-level hierarchical intercept
# nest- (observation-) level overdispersion residual
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
# ZI:
# fixed effects of distance to HWL and distance to dune
## @knitr zib_anest2-fit
zib_anest2 <- brm(bf(emerged | trials(clutch) ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                       (1 | fyear) + (1 | nestID),
                     zi ~ dist_hwl_std + dist_dune_std),
                  data = nest, family = zero_inflated_binomial(), 
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr zib_anest2-summary
summary(zib_anest2)
## @knitr

# ## Model selection using loo
# loo_nest <- loo(zib_anest0, zib_anest1, zib_anest2, moment_match = TRUE, compare = TRUE)
# loo_nest
# print(loo_compare(loo_nest)[,], 3)

#--------------------------------------------#
# Only nests where a female was encountered  #
# (exclude TEQ b/c only 2 encounters         #
#  NOTE brms subset() not working)           #
#--------------------------------------------#

# Binomial:
# zone-level hierarchical intercept
# year-level hierarchical intercept
# turtle-level hierarchical intercept
# nest- (observation-) level overdispersion residual
## @knitr zib_enest0-fit
zib_enest0 <- brm(emerged | trials(clutch) ~ (1 | zone) + (1 | fyear) + (1 | name) + (1 | nestID),
                  data = subset(nest, beach != "TEQ"), family = zero_inflated_binomial(), 
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr zib_enest0-summary
summary(zib_enest0)
## @knitr

# Binomial:
# year-level hierarchical intercept
# turtle-level hierarchical intercept
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
## @knitr zib_enest1-fit
zib_enest1 <- brm(emerged | trials(clutch) ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                    neophyte + (1 | fyear) + (1 | name) + (1 | nestID),
                  data = subset(nest, beach != "TEQ"), family = zero_inflated_binomial(), 
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr zib_enest1-summary
summary(zib_enest1)
## @knitr

# Binomial:
# year-level hierarchical intercept
# turtle-level hierarchical intercept
# linear trend
# fixed effects of beach, distance to HWL and distance to dune
# ZI:
# fixed effects of distance to HWL and distance to dune
## @knitr zib_enest2-fit
zib_enest2 <- brm(bf(emerged | trials(clutch) ~ year_ctr + beach + dist_hwl_std + dist_dune_std + 
                    neophyte + (1 | fyear) + (1 | name) + (1 | nestID),
                    zi ~ dist_hwl_std + dist_dune_std),
                  data = subset(nest, beach != "TEQ"), family = zero_inflated_binomial(), 
                  chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
## @knitr zib_enest2-summary
summary(zib_enest2)
## @knitr


#================================================================
# Diagnostic plots 
#================================================================

mod_name <- "sn_doy1"
## @knitr ppd
mod <- get(mod_name)
yrep <- posterior_predict(mod, cores = 1)
indx <- sample(nrow(yrep), 100)
y <- switch(class(mod)[1], stanreg = as.matrix(rstanarm::get_y(mod))[,1], brmsfit = brms::get_y(mod))
form <- switch(class(mod)[1], stanreg = formula(mod), brmsfit = formula(mod)[[1]])
grp <- switch(class(mod)[1], stanreg = mod$glmod$fr$fyear, brmsfit = mod$data$fyear)
## @knitr

#----------------------------------------------------------------
# Continuous responses
#----------------------------------------------------------------

# PPD marginal density
ppc_dens_overlay(y, yrep[indx,]) + ggtitle(deparse(form, width.cutoff = 500))

ggsave(filename=here("analysis", "results", paste0(mod_name, "_ppc_dens_overlay.png")),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

# PPD marginal density grouped by year (takes a while)
ppc_dens_overlay_grouped(y, yrep[indx,], group = grp) +
  ggtitle(deparse(form, width.cutoff = 500))

# PPD scatterplots grouped by year
ppc_scatter_avg_grouped(y, yrep, group = grp) + geom_abline(intercept = 0, slope = 1) + 
  ggtitle(deparse(form, width.cutoff = 500))

#----------------------------------------------------------------
# Binomial responses
#----------------------------------------------------------------

# PPD marginal histogram for binomial responses
## @knitr ppc_rootogram
ppc_rootogram(y, yrep[1:5,]) +   
  ggtitle(deparse(form, width.cutoff = 500))
## @knitr
ggsave(filename=here("analysis", "results", paste0(mod_name, "_ppc_rootogram.png")),
       width=7, height=7, units="in", dpi=300, type="cairo-png")

# PPD marginal histogram for binomial responses
ppc_hist(y/rowSums(y), sweep(yrep[1:5,], 2, rowSums(y), "/"), binwidth = 0.05) + 
  scale_x_continuous(limits = c(0,1)) +
  ggtitle(deparse(form, width.cutoff = 500))

# PPD empirical CDF overlay for binomial responses
ppc_ecdf_overlay(y/rowSums(y), sweep(yrep[indx,], 2, rowSums(y), "/")) +
  ggtitle(deparse(form, width.cutoff = 500))

# PPC for proportion zeros
ppc_stat(y, yrep[indx,], stat = function(x) mean(x == 0)) +
  ggtitle(deparse(form, width.cutoff = 500))

#----------------------------------------------------------------
# All responses
#----------------------------------------------------------------

# Normal QQ plot of year-level random effects
as.data.frame(ranef(mod)$fyear) %>% rename(intercept = 1) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + ggtitle(deparse(form, width.cutoff = 500))

# Normal QQ plot of turtle-level random effects
as.data.frame(ranef(mod)$name) %>% rename(intercept = 1) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), plot.title = element_text(size = 11)) +
  xlab("Normal quantiles") + ylab("Quantiles of turtle-specific intercepts") +
  ggtitle(deparse(form, width.cutoff = 500))

ggsave(filename=here("analysis", "results", paste0(mod_name, "_turtle-intercept_qqnorm.png")),
       width=7, height=7, units="in", dpi=300, type="cairo-png")


#================================================================
# Save stanreg and loo objects 
#================================================================

save(list = ls()[sapply(mget(ls()), function(x)
  any(class(x) %in% c("stanreg","brmsfit")) | any(sapply(x, class) == "loo"))], 
  file = here("analysis","results","leatherback_nesting_models.RData"))


