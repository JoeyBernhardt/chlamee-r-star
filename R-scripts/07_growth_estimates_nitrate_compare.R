

#' ---
#' title: "compare growth estimates"
#' author: "Joey"
#' ---

#+
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(cache = TRUE)

#+ knitr::opts_chunk$set(message = FALSE)

library(tidyverse)
library(cowplot)
library(broom)
library(readxl)
library(janitor)
library(plotrix)
library(here)



treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629") 


# compare growth rate estimates -------------------------------------------

#+ knitr::opts_chunk$set(warning = FALSE)

logistic <- read_csv(here("data-processed", "nitrate_exp_growth_logistic.csv")) %>% 
	select(population, ancestor_id, treatment, nitrate_concentration, estimate) %>% 
	mutate(method  = "logistic") %>% 
	distinct(population, ancestor_id, treatment, nitrate_concentration, estimate) %>% 
	filter(!is.na(ancestor_id)) %>% 
	mutate(population = as.character(population)) %>% 
	mutate(method = "intrinsic rate of incr via logistic")

eyeball <- read_csv(here("data-processed", "nitrate_exp_growth_eyeball.csv")) %>% 
	select(population, ancestor_id, treatment, nitrate_concentration.x, estimate) %>%
	rename(nitrate_concentration = nitrate_concentration.x) %>% 
	distinct(population, ancestor_id, treatment, nitrate_concentration, estimate) %>% 
	filter(!is.na(population)) %>% 
	mutate(population = as.character(population)) %>% 
	mutate(method  = "eyeball-exponential")

AIC_exp <- read_csv(here("data-processed", "nitrate_exp_growth_w_growthtools_AIC.csv")) %>% 
	rename(estimate = mu) %>% 
	mutate(method = "exponential via pw regression and AIC") %>% 
	mutate(population = as.character(population))
AIC_exp <- left_join(AIC_exp, treatments)

AICc_exp <- read_csv(here("data-processed", "nitrate_exp_growth_w_growthtools_AICc.csv")) %>% 
	rename(estimate = mu) %>% 
	mutate(method = "exponential via pw regression and AICc") %>% 
	mutate(population = as.character(population))
AICc_exp <- left_join(AICc_exp, treatments)

all_growth_estimates <- bind_rows(logistic, eyeball, AIC_exp, AICc_exp)

#+ fig.width=17, fig.height=6
all_growth_estimates %>%
	filter(method %in% c("intrinsic rate of incr via logistic", "eyeball-exponential")) %>% 
	group_by(population, nitrate_concentration) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate, color = method)) + geom_point(shape = 1) +
	facet_grid(ancestor_id ~ treatment)
ggsave(here("figures", "comparison_nitrate_monod.png"), width = 15, height = 6)

#+ fig.width=17, fig.height=6
all_growth_estimates %>% 
	group_by(population, nitrate_concentration, ancestor_id, treatment, method) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(aes(x = nitrate_concentration, y = mean, color = method)) + geom_point(shape=1) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	facet_grid(ancestor_id ~ treatment) + ylab("Mean growth rate (/day)") + xlab("Nitrate (uM)")
ggsave(here("figures", "comparison_nitrate_monod_means.png"), width = 17, height = 6)

AICc_sel <- AICc_exp %>% 
	select(well_plate, best.model, estimate, nitrate_concentration, n_obs) %>% 
	rename(best_model_AICc = best.model,
		   growth_rate_AICc = estimate,
		   n_obs_AICc = n_obs)

AIC_sel <- AIC_exp %>% 
	select(well_plate, best.model, estimate, n_obs) %>% 
	rename(best_model_AIC = best.model,
		   growth_rate_AIC = estimate,
		   n_obs_AIC = n_obs)

all_AIC <- left_join(AICc_sel, AIC_sel)



# Compare growth estimates via different IC methods -----------------------

#' From Colin: For the purposes of comparing/contrasting estimation approaches, as a simplifying step it might be worth pooling all of your observations
#' across nitrate levels and treatments, and making pair-wise scatter plots of method 1 vs. 2, 1 vs. 3, etc., while adding a reference line of y=x. Makes it easier to see things like potential biases. For example, it looks like the piecewise approach with AIC leads to pretty consistently higher estimates (as you said before). 
#' Skimming plots you sent yesterday, I suspect some of this comes from cases where there are only 3 observations in the ‘exponential’ phase of lag-sat models, and very little difference between the lag and saturating abundances (see attached cherry-picked example). From my perspective, these are pretty border-line situations, and tweaking the algorithm will only get us so far - ideally, we’d like better experimental data, perhaps with lower starting abundances and/or higher nutrient concentrations, so populations remain in exponential phase longer. 
#' I guess this is harder to achieve with many replicates and an explicit nutrient treatment. JB -- yes!


### here we see that when we use AICc, we always get equal or smaller growth rate estimates

#+ fig.width=8, fig.height=6
all_AIC %>% 
	ggplot(aes(x = growth_rate_AIC,  y = growth_rate_AICc, color = best_model_AIC)) +
	geom_point(shape = 1) +
	geom_abline(slope = 1, intercept = 0) +
	# geom_text(aes(label = n_obs)) +
	ylab("Growth rate chosen via AICc") + xlab("Growth rate chosen via AIC")


## Here the points are numbers corresponding to the # of points in the exponential phase
#+ fig.width=8, fig.height=6
p <- ggplot(all_AIC, aes(x = growth_rate_AIC, y =  growth_rate_AICc, color = best_model_AIC, label = n_obs_AIC))
p + geom_text() +
	geom_abline(slope = 1, intercept = 0) 

#+ fig.width=8, fig.height=6
all_AIC %>% 
	ggplot(aes(x = growth_rate_AIC,  y = growth_rate_AICc, color = best_model_AICc)) + geom_point(shape = 1) +
	geom_abline(slope = 1, intercept = 0) +
	ylab("Growth rate chosen via AICc") + xlab("Growth rate chosen via AIC")

## Here the points are numbers corresponding to the # of points in the exponential phase
#+ fig.width=8, fig.height=6
p <- ggplot(all_AIC, aes(x = growth_rate_AIC, y =  growth_rate_AICc, color = best_model_AICc, label = n_obs_AICc))
p + geom_text() +
	geom_abline(slope = 1, intercept = 0) 


## Here we see the breakdown of number of points in the exponential, by IC and model
#+ fig.width=8, fig.height=6
all_growth_estimates %>% 
	filter(grepl("AIC", method)) %>% 
	group_by(IC_method, best.model, n_obs) %>% 
	tally() %>% 
	ggplot(aes(x = n_obs, y = n)) + geom_bar(stat = "identity") +
	facet_grid(IC_method ~ best.model, scales = "free_y") +xlab("Number of points in the exponential phase") +
	ylab("Number of time series")



#' From Colin: Might also be worth making a plot of growth rate estimate vs. nutrient level (pooled across your other treatment variables) 
#' with points color-coded by the identity of the best growth rate model - it would help reveal whether the lag-sat model
#' is consistently favored at low Nitrate with the AIC approach.


#+ fig.width=8, fig.height=4
## with all the data shown
all_growth_estimates %>% 
	group_by(population, nitrate_concentration, ancestor_id, treatment, method) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(aes(x = nitrate_concentration, y = mean, color = method)) + geom_point(shape=1) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	# facet_grid(ancestor_id ~ treatment) + 
	ylab("Mean growth rate (/day)") + xlab("Nitrate (uM)")

#+ fig.width=8, fig.height=4
## with just the means at each nitrate level

all_growth_estimates %>% 
	group_by(nitrate_concentration, method) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(aes(x = nitrate_concentration, y = mean, color = method)) + geom_point(shape=1) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	# facet_grid(ancestor_id ~ treatment) + 
	ylab("Mean growth rate (/day)") + xlab("Nitrate (uM)")


#+ fig.width=8, fig.height=4
## At low nitrate levels, lag-sat consistently is associated with higher growth rates

all_growth_estimates %>% 
	filter(grepl("AIC", method)) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate, color = best.model, shape = method)) + geom_jitter(width = 30) +
	ylab("Growth rate (/day)") + xlim(0, 40) +xlab("Nitrate (uM)")

#+ fig.width=8, fig.height=4
## At high nitrate levels, lag-sat consistently is associated with higher growth rates

all_growth_estimates %>% 
	filter(grepl("AIC", method)) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate, color = best.model, shape = method)) + geom_jitter(width = 30) +
	ylab("Growth rate (/day)") + xlim(40, 1000) +xlab("Nitrate (uM)")


## Here's a look at how the frequency of the different best models varies by IC and nitrate level

all_growth_estimates %>% 
	filter(grepl("AIC", method)) %>% 
	mutate(method = ifelse(grepl("AICc", method), "AICc", "AIC")) %>% 
	group_by(nitrate_concentration, method, best.model) %>% 
	distinct(estimate) %>% 
	tally() %>% 
	spread(key = best.model, value = n) %>% 
	knitr::kable(align = "c")

## Here, each row corresponds to a nitrate concentration (ranging from 5uM to 1000uM)
## you can see how the distribution of best models changes with nitrate level
#+ fig.width=8, fig.height=12
all_growth_estimates %>% 
	filter(grepl("AIC", method)) %>% 
	mutate(method = ifelse(grepl("AICc", method), "AICc", "AIC")) %>% 
	group_by(nitrate_concentration, method, best.model) %>% 
	distinct(estimate) %>% 
	tally() %>% 
	ggplot(aes(x = best.model, y = n, fill = best.model)) + geom_bar(stat = "identity", width = 1) +
	facet_grid(nitrate_concentration ~ method) +ylab("Number of cases") + xlab("Best model")



