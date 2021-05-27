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

## Annual climatologies of each variable
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
# Somatic growth: change in max curved carapace length 
#                 for turtles encountered multiple times
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

## von Bertalanffy growth curves
## Overlay data on hyper-mean and female-specific predicted curves
mod_name <- "brm_vb0"
mod <- get(mod_name)

ce_vb_hyper <- conditional_effects(mod, effects = "dyear")
ce_vb_ranef <- conditional_effects(mod, effects = "dyear", re_formula = NULL)
dat <- size %>% filter(!is.na(dyear)) %>% group_by(name, year) %>% 
  summarize(ccl_max0 = ccl_max0, fyear = fyear, dyear = seq(0, dyear, length = 10)) %>%
  ungroup()
mr_vb <- posterior_epred(mod, newdata = dat)
dat$ccl_max <- colMedians(mr_vb)

dev.new()

size %>% 
  ggplot(aes(x = dyear, y = ccl_max)) +
  geom_ribbon(aes(x = dyear, ymin = lower__, ymax = upper__), data = ce_vb_ranef$dyear,
              fill = "lightgray", alpha = 0.5) +
  geom_ribbon(aes(x = dyear, ymin = lower__, ymax = upper__), data = ce_vb_hyper$dyear,
              fill = "gray", alpha = 0.7) +
  geom_line(aes(x = dyear, y = estimate__), data = ce_vb_hyper$dyear, lwd = 1, col = "darkgray") +
  geom_line(aes(x = dyear, y = ccl_max, group = name), data = dat, col = "darkgray") +
  geom_jitter(width = 0.2, height = 0, col = "steelblue4", alpha = 0.5) +
  scale_x_continuous(breaks = unique(mod$data$dyear)) +
  xlab("Years since encounter") + ylab(bquote(CCL[max] ~ "(cm)")) + theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

ggsave(filename=here("analysis", "results", "vonB_fitted_observed.png"),
       width=7, height=7, units="in", dpi=300, type="cairo-png")


#----------------------------------------------------------------
# Proportion of neophytes
#----------------------------------------------------------------

#$ Time series of proportion neophytes
## Show data + CI and posterior expectation + PPD
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

## Probability of nest failure vs distance from HWL and dune
## note PPD not shown b/c credible interval always c(0,1)
mod_name <- "zib_anest2"
mod <- get(mod_name)

ce_zi_epred <- conditional_effects(mod, effects = c("dist_hwl_std","dist_dune_std"),
                                   dpar = "zi", method = "fitted")
ce_zi_epred$dist_hwl_std <- ce_zi_epred$dist_hwl_std %>% 
  mutate(dist_hwl = dist_hwl_std * attr(nest$dist_hwl_std, "scaled:scale") + 
           attr(nest$dist_hwl_std, "scaled:center"), .after = dist_hwl_std)
ce_zi_epred$dist_dune_std <- ce_zi_epred$dist_dune_std %>% 
  mutate(dist_dune = dist_dune_std * attr(nest$dist_dune_std, "scaled:scale") + 
           attr(nest$dist_dune_std, "scaled:center"), .after = dist_dune_std)

dev.new(width = 10, height = 5)

nest %>% mutate(dhwl = cut(dist_hwl, 10)) %>% group_by(dhwl) %>% 
  summarize(dist_hwl = mean(dist_hwl), nz = sum(emergence_rate == 0, na.rm = TRUE), n = n()) %>%
  ungroup() %>% cbind(with(., Hmisc::binconf(x = nz, n = n, alpha = 0.1))) %>%
  ggplot(aes(x = dist_hwl, y = PointEst)) + 
  geom_line(aes(x = dist_hwl, y = estimate__), data = ce_zi_epred$dist_hwl_std,
            inherit.aes = FALSE, lwd = 1, col = "darkgray") +
  geom_ribbon(aes(x = dist_hwl, ymin = lower__, ymax = upper__), 
              data = ce_zi_epred$dist_hwl_std, inherit.aes = FALSE, 
              fill = "gray", alpha = 0.7) +
  geom_point(size = 2.5, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  coord_cartesian(ylim = c(0,0.5)) + 
  xlab("Distance from HWL (m)") + ylab("Probability of nest failure") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white")) +
  
nest %>% mutate(ddune = cut(dist_dune, 10)) %>% group_by(ddune) %>% 
  summarize(dist_dune = mean(dist_dune), nz = sum(emergence_rate == 0, na.rm = TRUE), 
            n = n(), .groups = "drop") %>%
  cbind(with(., Hmisc::binconf(x = nz, n = n, alpha = 0.1))) %>% 
  ggplot(aes(x = dist_dune, y = PointEst)) + 
  geom_line(aes(x = dist_dune, y = estimate__), data = ce_zi_epred$dist_dune_std,
            inherit.aes = FALSE, lwd = 1, col = "darkgray") +
  geom_ribbon(aes(x = dist_dune, ymin = lower__, ymax = upper__), 
              data = ce_zi_epred$dist_dune_std, inherit.aes = FALSE, 
              fill = "gray", alpha = 0.7) +
  geom_point(size = 2.5, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  coord_cartesian(ylim = c(0,0.5)) + xlab("Distance from dune (m)") + ylab("") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"))

ggsave(filename=here("analysis", "results", "p_nest-failure_dist_hwl_dune.png"),
       width=10, height=5, units="in", dpi=300, type="cairo-png")

## Overall emergence rate vs distance from HWL and dune by beach
mod_name <- "zib_anest2"
mod <- get(mod_name)

condxns <- data.frame(beach = levels(nest$beach))
ce_bin_epred <- conditional_effects(mod, effects = c("dist_hwl_std","dist_dune_std"),
                                    conditions = condxns, method = "fitted")
ce_bin_epred$dist_hwl_std <- ce_bin_epred$dist_hwl_std %>% 
  mutate(dist_hwl = dist_hwl_std * attr(nest$dist_hwl_std, "scaled:scale") + 
           attr(nest$dist_hwl_std, "scaled:center"), .after = dist_hwl_std)
ce_bin_epred$dist_dune_std <- ce_bin_epred$dist_dune_std %>% 
  mutate(dist_dune = dist_dune_std * attr(nest$dist_dune_std, "scaled:scale") + 
           attr(nest$dist_dune_std, "scaled:center"), .after = dist_dune_std)

dev.new()

nest %>% mutate(dhwl = cut(dist_hwl, 10)) %>% group_by(beach, dhwl) %>% 
  summarize(dist_hwl = mean(dist_hwl), emerged = sum(emerged, na.rm = TRUE), 
            clutch = sum(clutch, na.rm = TRUE)) %>%
  ungroup() %>% cbind(with(., Hmisc::binconf(x = emerged, n = clutch, alpha = 0.1))) %>%
  ggplot(aes(x = dist_hwl, y = PointEst)) + 
  geom_line(aes(x = dist_hwl, y = estimate__), data = ce_bin_epred$dist_hwl_std,
            inherit.aes = FALSE, lwd = 1, col = "darkgray") +
  geom_ribbon(aes(x = dist_hwl, ymin = lower__, ymax = upper__), 
              data = ce_bin_epred$dist_hwl_std, inherit.aes = FALSE, 
              fill = "gray", alpha = 0.7) +
  geom_point(size = 2.5, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  xlab("Distance from HWL (m)") + ylab("Emergence rate") +
  facet_wrap(vars(beach), ncol = 1) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white")) +

  nest %>% mutate(ddune = cut(dist_dune, 10)) %>% group_by(beach, ddune) %>% 
  summarize(dist_dune = mean(dist_dune), emerged = sum(emerged, na.rm = TRUE), 
            clutch = sum(clutch, na.rm = TRUE)) %>%
  ungroup() %>% cbind(with(., Hmisc::binconf(x = emerged, n = clutch, alpha = 0.1))) %>%
  ggplot(aes(x = dist_dune, y = PointEst)) + 
  geom_line(aes(x = dist_dune, y = estimate__), data = ce_bin_epred$dist_dune_std,
            inherit.aes = FALSE, lwd = 1, col = "darkgray") +
  geom_ribbon(aes(x = dist_dune, ymin = lower__, ymax = upper__), 
              data = ce_bin_epred$dist_dune_std, inherit.aes = FALSE, 
              fill = "gray", alpha = 0.7) +
  geom_point(size = 2.5, col = "steelblue4") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "steelblue4") +
  xlab("Distance from dune (m)") + ylab("") +
  facet_wrap(vars(beach), ncol = 1) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"))
  
ggsave(filename=here("analysis", "results", "p_emergence_dist_hwl_dune.png"),
       width=7, height=7, units="in", dpi=300, type="cairo-png")



