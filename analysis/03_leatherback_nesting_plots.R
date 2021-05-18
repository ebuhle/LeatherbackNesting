#################################################################
# LEATHERBACK SEA TURTLE NESTING FIGURES
#################################################################

# Possible exploratory plots to show interesting patterns
#
# * Number of nests per year (shows rebound after decline from 2015-2017
#
# * Number of beaches visited by females
# nest %>% group_by(name) %>% summarize(n_beach = length(unique(beach))) %>% 
#   group_by(n_beach) %>% summarize(count = n())
#
# * Do females have a typical nesting height as well as date?
# * Do females show site-fidelity to a specific beach? How much beach-switching
#   occurs w/in or b/w years?


#----------------------------------------------------------------
# Weather data
#----------------------------------------------------------------

# Annual climatologies of each variable
dev.new(width = 10, height = 5)

weather %>% group_by(date) %>% select(-dp_avg) %>% 
  summarize(across(c(ends_with("avg"), ppt, sst), mean, na.rm = TRUE), .groups = "drop") %>% 
  pivot_longer(-date, names_to = "variable") %>% 
  mutate(variable = recode_factor(variable, 
                                  humid_avg = "Humidity~('%')", ws_avg = "Wind~speed~(mph)",
                                  p_avg = "Pressure~'(in Hg)'", t_avg = "Temperature~(degree * C)", 
                                  dp_avg = "Dew~point~(degree * C)", ppt = "Precipitation~(cm)", 
                                  sst = "SST~(degree * C)")) %>% 
  mutate(doy = yday(date), .after = date) %>% group_by(variable, doy) %>% 
  summarize(lb = quantile(value, 0.05, na.rm = TRUE), med = median(value, na.rm = TRUE), 
            ub = quantile(value, 0.95, na.rm = TRUE)) %>% 
  ggplot(aes(x = as_date(as_date(doy), format = "%m-%d"), y = med)) +
  annotate("rect", xmin = min(as_date(turtle$doy_encounter), na.rm = TRUE),
           xmax = max(as_date(turtle$doy_encounter), na.rm = TRUE),
           ymin = -Inf, ymax = Inf, fill = "steelblue4", alpha = 0.5) +
  geom_line(lwd = 0.7, col = "gray40") +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "gray40", alpha = 0.7) +
  scale_x_date(date_minor_breaks = "1 month", date_labels = "%b") + 
  xlab("Date") + ylab("Daily average") +
  facet_wrap(vars(variable), ncol = 3, scales = "free_y", labeller = label_parsed) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = NA))

ggsave(filename=here("analysis", "results", "climatology.png"),
       width=10, height=5, units="in", dpi=300, type="cairo-png")


#----------------------------------------------------------------
# Breeding phenology: encounter date
#----------------------------------------------------------------

## Time series of encounter DOY, all data
## Split violin plots show data distribution vs. PPD for each year
mod_name <- "lmer_doy2"
mod <- get(mod_name)
yrep <- posterior_predict(mod)
indx <- sample(nrow(yrep), 100)

ppd <- as.data.frame(t(yrep[indx,])) %>% mutate(fyear = turtle$fyear) %>% 
  pivot_longer(-fyear, names_to = "iter", values_to = "doy_yrep") %>% 
  mutate(date_encounter = as_date(doy_yrep), "%m-%d")

dev.new(width = 7, height = 5)

turtle %>% 
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

## Time series of encounter DOY, female averages
dev.new(width = 7, height = 5)

turtle %>% group_by(year, name) %>% summarize(mean_date_encounter = mean(date_encounter), .groups = "drop") %>% 
  ggplot(aes(x = year, y = as_date(format(mean_date_encounter, "%m-%d"), format = "%m-%d"), group = name)) +
  geom_line(col = "steelblue4", alpha = 0.5) + scale_x_continuous(breaks = sort(unique(turtle$year))) +
  xlab("Year") + ylab("Encounter DOY") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

ggsave(filename=here("analysis", "results", "doy_female_avg.png"),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

## Joyplot of predicted DOY for each female in an average year
## Only include females encountered in >= 3 years (relatively informative for intercept)
mod_name <- "lmer_doy2"
mod <- get(mod_name)

dat <- turtle %>% group_by(name) %>% 
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
  geom_density_ridges(col = "lightgray", fill = "steelblue4", 
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_x_date(date_breaks = "2 weeks", date_minor_breaks = "1 week", date_labels = "%b %d",
               limits = function(x) c(x[1] + 7, x[2] - 5)) +
  xlab("Encounter DOY") + 
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11, margin = margin(r = -10)),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.line.x = element_line(size = 0.5), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename=here("analysis", "results", "doy_female_joyplot.png"),
       width=7*0.9, height=10*0.9, units="in", dpi=300, type="cairo-png")


#----------------------------------------------------------------
# Somatic growth: 
# change in max curved carapace length / year 
# for turtles encountered multiple times
#----------------------------------------------------------------

## Time series of body size for each female
dev.new(width = 7, height = 5)

turtle %>% group_by(year, name) %>% summarize(mean_ccl_max = mean(ccl_max), .groups = "drop") %>% 
  ggplot(aes(x = year, y = mean_ccl_max, group = name)) +
  geom_line(col = "steelblue4", alpha = 0.7) + 
  scale_x_continuous(breaks = sort(unique(turtle$year))) +
  xlab("Year") + ylab("Max curved carapace length (cm)") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

ggsave(filename=here("analysis", "results", "ccl_max_timeseries.png"),
       width=7, height=5, units="in", dpi=300, type="cairo-png")

## Specific growth rate as a function of initial length
## Overlay data on lmer hyper-mean regression line and female-specific lines
mod_name <- "lmer_sgr2"
mod <- get(mod_name)

dat <- data.frame(name = "0", fyear = "0", 
                  ccl_max0_std = seq(1.1*min(mod$glmod$fr$ccl_max0_std),
                                     1.1*max(mod$glmod$fr$ccl_max0_std),
                                     length = 200))
pred_hyper <- posterior_linpred(mod, newdata = dat, re.form = NA)
pred_group <- posterior_linpred(mod, newdata = dat)
ccl_ref <- size$ccl_max0[!is.na(size$sgr)]
ccl_mean <- mean(ccl_ref)
ccl_sd <- sd(ccl_ref)
dat <- dat %>% 
  mutate(ccl_max0 = ccl_max0_std * ccl_sd + ccl_mean, med_hyper = colMedians(pred_hyper), 
         lb_hyper = colQuantiles(pred_hyper, probs = 0.025), 
         ub_hyper = colQuantiles(pred_hyper, probs = 0.975), 
         lb_group = colQuantiles(pred_group, probs = 0.025), 
         ub_group = colQuantiles(pred_group, probs = 0.975))

dev.new(width = 7, height = 7)

mod$glmod$fr %>% cbind(ccl_max0 = ccl_ref) %>% 
  ggplot(aes(x = ccl_max0, group = name)) +
  geom_ribbon(aes(x = ccl_max0, ymin = lb_group, ymax = ub_group), data = dat, 
              fill = "lightgray", alpha = 0.5) +
  geom_ribbon(aes(x = ccl_max0, ymin = lb_hyper, ymax = ub_hyper), data = dat, 
              fill = "gray", alpha = 0.7) +
  geom_line(aes(x = ccl_max0, y = med_hyper), data = dat, col = "darkgray", lwd = 1) +
  geom_line(aes(y = sgr), col = "steelblue4", alpha = 0.5) +
  geom_point(aes(y = sgr), col = "steelblue4", alpha = 0.5) +
  scale_y_continuous(n.breaks = 7) + theme_bw(base_size = 16) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  xlab(bquote("Initial" ~ CCL[max] ~ "(cm)")) + 
  ylab(bquote("Specific growth rate (" * decade^-1 * " )"))

ggsave(filename=here("analysis", "results", "sgr_vs_init_ccl.png"),
       width=7, height=7, units="in", dpi=300, type="cairo-png")


#----------------------------------------------------------------
# Proportion of neophytes
#----------------------------------------------------------------

# Time series of proportion neophytes
# Show data + CI and posterior expectation + PPD
mod_name <- "glm_neo1"
mod <- get(mod_name)
epred <- posterior_epred(mod)
yrep <- sweep(posterior_predict(mod), 2, neophyte$count, "/")

dev.new(width = 7, height = 5)

neophyte %>% 
  ggplot(aes(x = year, y = p_neophyte)) +
  geom_ribbon(aes(ymin = colQuantiles(yrep, probs = 0.025),
                  ymax = colQuantiles(yrep, probs = 0.975)),
              fill = "lightgray", alpha = 0.5) +
  geom_ribbon(aes(ymin = colQuantiles(epred, probs = 0.025),
                  ymax = colQuantiles(epred, probs = 0.975)),
              fill = "gray", alpha = 0.7) +
  geom_line(y = colMedians(epred), lwd = 1, col = "darkgray") +
  geom_point(size = 4, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  scale_x_continuous(breaks = neophyte$year) + scale_y_continuous(n.breaks = 8) +
  xlab("Year") + ylab("Proportion neophytes") +
  theme(panel.grid = element_blank())
  
ggsave(filename=here("analysis", "results", "p_neophyte_timeseries.png"),
       width=7, height=5, units="in", dpi=300, type="cairo-png")


#----------------------------------------------------------------
# Hatchling emergence success
#----------------------------------------------------------------

# Probability of zero emergence vs distance from HWL by beach
## ADD binomial CIs to points
dev.new(width = 4, height = 8)

nest %>% mutate(dhwl = cut(dist_hwl, 15)) %>% group_by(beach, dhwl) %>% 
  summarize(dist_hwl = mean(dist_hwl), n_zero = sum(emergence_rate = 0, na.rm = TRUE), n = n()) %>% 
  ungroup() %>% cbind(with(., Hmisc::binconf(x = n_zero, n = n, alpha = 0.1))) %>%
  ggplot(aes(x = dist_hwl, y = PointEst)) + 
  geom_point(size = 2.5, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  xlab("Distance from HWL") + ylab("P(emergence = 0)") +
  facet_wrap(vars(beach), ncol = 1, scales = "free_y") +
  theme(strip.background = element_rect(fill = "white"))

# Mean of nonzero emergence rates vs distance from HWL by beach
## ADD SEs to points (will be busy; use vertical segments only?)
## Use probability y-axis
dev.new(width = 4, height = 8)

nest %>% filter(emergence_rate > 0) %>% 
  ggplot(aes(x = dist_hwl, y = qlogis(emergence_rate))) + 
  geom_point(size = 2, col = "steelblue4", alpha = 0.5) + 
  facet_wrap(vars(beach), ncol = 1)

# Probability of zero emergence vs. distance from dune by beach
## ADD binomial CIs to points
dev.new(width = 4, height = 8)

nest %>% mutate(ddune = cut(dist_dune, 15)) %>% group_by(beach, ddune) %>% 
  summarize(dist_dune = mean(dist_dune), n_zero = sum(emergence_rate == 0, na.rm = TRUE), 
            n = n(), .groups = "drop") %>% 
  cbind(with(., Hmisc::binconf(x = n_zero, n = n, alpha = 0.1))) %>% 
  ggplot(aes(x = dist_dune, y = PointEst)) + 
  geom_point(size = 2.5, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  xlab("Distance from toe of dune") + ylab("P(emergence = 0)") +
  facet_wrap(vars(beach), ncol = 1, scales = "free_y") +
  theme(strip.background = element_rect(fill = "white"))

# Mean of nonzero emergence rates vs distance from HWL by beach
## ADD SEs to points (will be busy; use vertical segments only?)
## Use probability y-axis
dev.new(width = 4, height = 8)

nest %>% filter(emergence_rate > 0) %>% 
  ggplot(aes(x = dist_dune, y = qlogis(emergence_rate))) + 
  geom_point(size = 2, col = "steelblue4", alpha = 0.5) + 
  facet_wrap(vars(beach), ncol = 1)




