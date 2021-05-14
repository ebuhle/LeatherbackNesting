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
summary(lmer_doy0, pars = "alpha", regex_pars = "igma", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy0, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# linear trend
lmer_doy1 <- stan_lmer(doy_encounter ~ year_ctr + (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy1)
summary(lmer_doy1, pars = c("alpha","beta"), regex_pars = "igma", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy1, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# linear trend
# all weather variables (average annual anomalies)
# except dew point (correlation with temperature is 0.91)
lmer_doy2 <- stan_lmer(doy_encounter ~ year_ctr + humid_std + ws_std + p_std + t_std + 
                         ppt_std + sst_std + (1 | name) + (1 | fyear), data = turtle, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_doy2)
summary(lmer_doy2, pars = c("alpha","beta"), regex_pars = "igma", probs = c(0.025, 0.5, 0.975))
summary(lmer_doy2, pars = "varying")

# # DOY of a female's *first* encounter in a given year
# # turtle-level and year-level hierarchical intercepts
# # linear change over time (year)
# # => very similar to lmer_doy but fewer obs/turtle so more uncertainty
# lmer_doy1st1 <- stan_lmer(doy_encounter ~ year_ctr + (1 | name) + (1 | fyear),
#                          data = turtle, subset = first_of_year,
#                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
# print(lmer_doy1st1)


#----------------------------------------------------------------
# Somatic growth: 
# change in max curved carapace length / year 
# for turtles encountered multiple times
#----------------------------------------------------------------

# turtle-level and year-level hierarchical intercepts
lmer_sgr0 <- stan_lmer(sgr ~ (1 | name) + (1 | fyear), data = size, 
                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr0, 3)
summary(lmer_sgr0, pars = "alpha", regex_pars = "igma", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr0, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# initial size effect
lmer_sgr1 <- stan_lmer(sgr ~ ccl_max0_std + (1 | name) + (1 | fyear), data = size, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr1, 3)
summary(lmer_sgr1, pars = "alpha", regex_pars = "igma", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr1, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# initial size effect with turtle-varying slopes 
# (diagonal random effects covariance matrix)
lmer_sgr2 <- stan_lmer(sgr ~ ccl_max0_std + (ccl_max0_std || name) + (1 | fyear), data = size, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr2, 3)
summary(lmer_sgr2, pars = "alpha", regex_pars = "igma", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr2, pars = "varying")

# turtle-level and year-level hierarchical intercepts
# initial size effect with turtle-varying slopes 
# linear time trend
lmer_sgr3 <- stan_lmer(sgr ~ ccl_max0_std + year_ctr + (ccl_max0_std || name) + (1 | fyear), data = size, 
                       chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(lmer_sgr3, 3)
summary(lmer_sgr3, pars = "alpha", regex_pars = "igma", probs = c(0.025, 0.5, 0.975), digits = 3)
summary(lmer_sgr3, pars = "varying")


#----------------------------------------------------------------
# Proportion of neophytes
#----------------------------------------------------------------

# intercept only
glm_neo0 <- stan_glm(cbind(neophyte, remigrant) ~ 1, 
                     data = neophyte, family = binomial, 
                     chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glm_neo0, 2)
summary(glm_neo0, probs = c(0.025, 0.5, 0.975))

# linear trend
glm_neo1 <- stan_glm(cbind(neophyte, remigrant) ~ year_ctr, 
                     data = neophyte, family = binomial, 
                     chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
print(glm_neo1, 2)
summary(glm_neo1, probs = c(0.025, 0.5, 0.975))

## Little evidence of overdispersion
## GLMMs also had lots of bad Pareto k estimates
# # intercept only
# # observation-level (year) residuals
# glmer_neo0 <- stan_glmer(cbind(neophyte, remigrant) ~ (1 | fyear), 
#                          data = neophyte, family = binomial, 
#                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
# print(glmer_neo0, 2)
# summary(glmer_neo0, probs = c(0.025, 0.5, 0.975))
# 
# # linear trend
# # observation-level (year) residuals
# glmer_neo1 <- stan_glmer(cbind(neophyte, remigrant) ~ year_ctr + (1 | fyear), 
#                          data = neophyte, family = binomial, 
#                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000)
# print(glmer_neo1, 2)
# summary(glmer_neo1, probs = c(0.025, 0.5, 0.975))

## Model selection using loo
stanreg_list_neo <- stanreg_list(glm_neo0, glm_neo1)
loo_neo <- lapply(stanreg_list_neo, loo)
loo_neo
loo_compare(loo_neo)


#----------------------------------------------------------------
# Diagnostic plots 
#----------------------------------------------------------------

mod_name <- "glm_neo1"
mod <- get(mod_name)
yrep <- posterior_predict(mod)
indx <- sample(nrow(yrep), 100)

# PPD marginal density
ppc_dens_overlay(mod$y, yrep[indx,]) + ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))
ggsave(filename=here("analysis", "results", paste0(mod_name, "_ppc_dens_overlay.png")),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

# PPD marginal histogram
ppc_hist(mod$y[,"neophyte"], yrep[1:7,], binwidth = diff(range(yrep)) / 15) + 
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPD marginal density grouped by year (takes a while)
ppc_dens_overlay_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$year) +
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

# PPD scatterplots grouped by year
ppc_scatter_avg_grouped(mod$y, yrep, group = mod$glmod$fr$fyear) +
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle(deparse(mod$glmod$formula, width.cutoff = 500))

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





