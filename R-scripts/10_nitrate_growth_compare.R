
#' ---
#' title: "Compare nitrate experiment growth estimates"
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
library(growthTools)
library(rootSolve)
library(nlstools)
library(nls.multstart)

source(here("R-scripts", "predict_monod.R"))
load(here("data-processed", "fits_many_n0.RData"))
load(here("data-processed", "fitted_logistic_nitrate.RData"))

#' Read in data
treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629")
nitrate <- read_csv(here("data-processed", "nitrate-abundances-processed.csv")) %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	filter(population != "COMBO") 

population4 <- nitrate %>% 
	filter(population == 4)

write_csv(population4, here("data-processed", "population4-nitrate.csv"))

#' Eyeball exponential 2 day cutoff 

nitrate_exp <- nitrate %>% 
	mutate(exponential = case_when(days < 2 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes") 

#+ fig.width = 12, fig.height = 8
# nitrate_exp %>% 
# 	ggplot(aes(x = days, y = RFU, color = factor(nitrate_concentration))) + geom_point() +
# 	facet_wrap( nitrate_concentration~ population, scales = "free_y") + scale_color_viridis_d()
# ggsave(here("figures", "nitrate-exp.png"), width = 45, height = 45)

growth_rates <- nitrate_exp %>%
	group_by(nitrate_concentration, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


#' Fit the Monod model to the growth rates
monod_fits <- growth_rates %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))

bs_split <- monod_fits %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_n <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population")

all_predsn_2 <- left_join(all_preds_n, treatments, by = c("population")) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)
all_growth_n2 <- left_join(growth_rates, treatments) ## changed this to the 'cut' version. can switch back later



#+ fig.width = 12, fig.height = 8
all_growth_n2 %>% 
	# mutate(estimate = mu) %>%
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	# filter(treatment == "N", ancestor_id == "anc3") %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	# geom_point(aes(x = nitrate_concentration.x, y = estimate), data = filter(nitrate_eyeball_exp, population == 27), color = "blue") +
	geom_line(data=all_predsn_2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	geom_hline(yintercept = 0.1, linetype = "dotted") +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")
ggsave(here("figures", "nitrate-monod-curves-eyeball-2day.png"), width = 12, height =8)


monod_wide <-  monod_fits %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
#' Find R*
rstars <- monod_wide %>% 
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

rstars2 <- left_join(rstars, treatments, by = "population") %>%
	distinct(population, ks, umax, .keep_all = TRUE)

#+ fig.width = 6, fig.height = 4
rstars2 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave(here("figures", "nitrate-r-star-exponential-2day.png"), width = 6, height = 4)



# sequential exponential fits ---------------------------------------------

fitting_window <- function(x) {
	
	growth_rates <- nitrate %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(nitrate_concentration, population, well_plate) %>%
		do(tidy(lm(log(RFU) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

windows <- seq(3,7, by = 1)

multi_fits <- windows %>% 
	map_df(fitting_window, .id = "iteration")

#+ fig.width = 8, fig.height = 6
multi_fits %>% 
	filter(term == "days") %>% 
	ggplot(aes(x = number_of_points, y = estimate, group = well_plate)) + geom_point() + geom_line() +
	facet_wrap( ~ nitrate_concentration)

exp_fits_top <- multi_fits %>%
	filter(term == "days") %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate)

write_csv(exp_fits_top, "data-processed/nitrate-sequential-fits.csv")

#' Fit the Monod model to the growth rates
monod_fits_seq <- exp_fits_top %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))

bs_split_seq <- monod_fits_seq %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_seq <- bs_split_seq %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population")

all_preds_seq2 <- left_join(all_preds_seq, treatments, by = c("population")) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)
all_growth_seq2 <- left_join(exp_fits_top, treatments) ## changed this to the 'cut' version. can switch back later



#+ fig.width = 12, fig.height = 8
all_growth_seq2 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	geom_line(data=all_preds_seq2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	geom_hline(yintercept = 0.1, linetype = "dotted") +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")
ggsave(here("figures", "nitrate-monod-curves-sequential-method.png"), width = 12, height =8)


monod_wide_seq <-  monod_fits_seq %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
#' Find R*
rstars_seq <- monod_wide_seq %>% 
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

rstars2_seq <- left_join(rstars_seq, treatments, by = "population") %>%
	distinct(population, ks, umax, .keep_all = TRUE)

#+ fig.width = 6, fig.height = 4
rstars2_seq %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2_seq, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave(here("figures", "nitrate-r-star-sequential-method.png"), width = 6, height = 4)


# logistic ----------------------------------------------------------------

fits_many_n0 <- nitrate %>% 
	group_by(population, well_plate) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ K/(1 + (K/N0- 1)*exp(-r*days)),
												  data = .x,
												  iter = 500,
												  start_lower = c(K = 100, r = 0),
												  start_upper = c(K = 500, r = 2),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(K = 10, r = 0),
												  upper = c(K = 4000, r = 4),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

save(fits_many_n0, file = here("data-processed", "fits_many_n0.RData"))


fits_many <- fits_many_n0

fits_well_plate <- fits_many %>% 
	select(well_plate)

key <- nitrate %>% 
	select(well_plate, population, nitrate_concentration) %>% 
	distinct(well_plate, .keep_all = TRUE)	

pops <- left_join(fits_well_plate, key)

info <- fits_many %>%
	unnest(fit %>% map(glance))

# get params
params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 

# write_csv(params, here("data-processed", "logistic-parameters-nitrate.csv"))

bad_fits <- params %>% 
	filter(term == "K", estimate == 4000) 

bad_fits2 <- left_join(bad_fits, key)


# CI <- fits_many %>% 
# 	filter(fit != "NULL") %>% 
# 	unnest(fit %>% map(~ confint2(.x) %>%
# 					   	data.frame() %>% 
# 					   	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
# 	group_by(., well_plate) %>%
# 	mutate(., term = c('r', 'K')) %>%
# 	ungroup()
# 
# # merge parameters and CI estimates
# params <- merge(params, CI, by = intersect(names(params), names(CI)))

new_preds <- nitrate %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 15), stringsAsFactors = FALSE))


# max and min for each curve
max_min <- group_by(nitrate, well_plate) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

# create new predictions
preds2b <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds)) %>%
	merge(., max_min, by = 'well_plate') %>%
	group_by(., well_plate) %>%
	filter(., days > unique(min_days) & days < unique(max_days)) %>%
	rename(., RFU = .fitted) %>%
	ungroup()



preds <- fits_many %>%
	unnest(fit %>% map(augment))

new_preds <- nitrate %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 15), stringsAsFactors = FALSE))


# max and min for each curve
max_min <- group_by(nitrate, well_plate) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

key <- nitrate %>% 
	select(well_plate, population, nitrate_concentration) %>% 
	distinct(well_plate, .keep_all = TRUE)

preds3 <- left_join(preds, key, by = c("well_plate", "population")) %>% 
	left_join(., treatments) %>% 
	distinct(.fitted, population, nitrate_concentration, well_plate, .keep_all = TRUE)
preds3b <- left_join(preds2b, key, by = c("well_plate", "population"))

save(preds3b, file = here("data-processed", "fitted_logistic_nitrate.RData"))
nitrate2 <- left_join(nitrate, treatments)

#+ fig.width = 45, fig.height = 30
ggplot() +
	geom_point(aes(days, RFU, color = factor(nitrate_concentration)), size = 2, data = nitrate2) +
	geom_line(aes(days, .fitted, group = well_plate, color = factor(nitrate_concentration)), data = preds3) +
	facet_grid(ancestor_id ~ treatment, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave(here("figures", "logisitic-fits-augment.png"), width = 45, height = 30)



names(preds3b)

preds3c <- preds3b %>% 
	rename(nitrate_concentration = nitrate_concentration.y)

#+ fig.width = 12, fig.height = 30
ggplot() +
	geom_point(aes(days, RFU, color = factor(nitrate_concentration)), size = 2, data = filter(nitrate, population == 27)) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(nitrate_concentration)), data = filter(preds3c, population == 27)) +
	facet_grid(nitrate_concentration ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave(here("figures", "nitrate_rfus_logistic_pop27.png"), width = 12, height = 35)


#+ fig.width = 45, fig.height = 30
ggplot() +
	geom_point(aes(days, RFU, color = factor(nitrate_concentration)), size = 2, data = nitrate) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(nitrate_concentration)), data = preds3c) +
	facet_grid(nitrate_concentration ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave(here("figures", "nitrate_rfus_logistic_all.png"), width = 45, height = 35)




# params <- read_csv(here("data-processed", "logistic-parameters-nitrate.csv")) %>% 
# 	mutate(population = as.character(population))

growth_rates_logistic <- left_join(params, treatments, by = "population") %>% 
	left_join(., key) %>% 
	filter(term == "r")


growth_rates_logistic %>% 
	filter(population == 27) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate)) + geom_point()

#### fit Monods

#' Fit the Monod model to the growth rates
monod_fits_logistic <- growth_rates_logistic %>% 
	filter(!well_plate %in% c(bad_fits$well_plate)) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))

bs_split_logistic <- monod_fits_logistic %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_logistic <- bs_split_logistic %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population")

all_preds_logistic2 <- left_join(all_preds_logistic, treatments, by = c("population")) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)
all_growth_logistic2 <- left_join(growth_rates_logistic, treatments) ## changed this to the 'cut' version. can switch back later



#+ fig.width = 12, fig.height = 8
all_growth_logistic2 %>% 
	# mutate(estimate = mu) %>%
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	# filter(treatment == "N", ancestor_id == "anc3") %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	# geom_point(aes(x = nitrate_concentration.x, y = estimate), data = filter(nitrate_eyeball_exp, population == 27), color = "blue") +
	geom_line(data=all_preds_logistic2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	geom_hline(yintercept = 0.1, linetype = "dotted") +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")
ggsave(here("figures", "nitrate-monod-curves-logistic.png"), width = 12, height =8)


# write_csv(all_preds_logistic2, here("data-processed", "all_preds_logistic2.csv"))

monod_wide_logistic <-  monod_fits_logistic %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
#' Find R*
rstars_logistic <- monod_wide_logistic %>% 
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

rstars2_logistic <- left_join(rstars_logistic, treatments, by = "population") %>%
	distinct(population, ks, umax, .keep_all = TRUE)

#+ fig.width = 6, fig.height = 4
rstars2_logistic %>% 
	group_by(treatment) %>% 
	# filter(population != 27) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2_logistic, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave(here("figures", "nitrate-r-star-logistic.png"), width = 6, height = 4)




# growthTools -------------------------------------------------------------

#' this is the step that gets us the growth rate estimates
growth_rates_n_AICc <- nitrate %>%
	mutate(ln.fluor = log(RFU)) %>% 
	group_by(well_plate) %>% 
	do(grs = get.growth.rate(x = .$days, y = .$ln.fluor,id = .$well_plate, plot.best.Q = F)) 

#' Get growth rates via AIC
growth_rates_n_AIC <- nitrate %>%
	mutate(ln.fluor = log(RFU)) %>% 
	group_by(well_plate) %>% 
	do(grs = get.growth.rate(x = .$days, y = .$ln.fluor,id = .$well_plate, plot.best.Q = F, model.selection = "AIC")) 


#' Pull out the things we want
growth_sum_n_AICc <- growth_rates_n_AICc %>%
	summarise(well_plate, mu = grs$best.slope, n_obs = grs$best.model.slope.n,
			  slope_r2 = grs$best.model.slope.r2,
			  best.model_r2 = grs$best.model.rsqr,
			  best.model = grs$best.model, best.se = grs$best.se,
			  contents = grs$best.model.contents) 

growth_sum_n_AIC <- growth_rates_n_AIC %>%
	summarise(well_plate, mu = grs$best.slope, n_obs = grs$best.model.slope.n,
			  slope_r2 = grs$best.model.slope.r2,
			  best.model_r2 = grs$best.model.rsqr,
			  best.model = grs$best.model, best.se = grs$best.se,
			  contents = grs$best.model.contents) 

all_models <- left_join(growth_sum_n_AICc, growth_sum_n_AIC, by = "well_plate")

sat_models <- growth_sum_n_AIC %>% 
	unnest(contents %>% map(tidy, .id = "number")) %>% 
	# filter(well_plate %in% c(mismatches$well_plate)) %>% 
	filter(best.model == "gr.sat", term == "B1") %>% 
	rename(cutoff_point = estimate) %>% 
	select(well_plate, cutoff_point, best.model)

lag_sat_models <- growth_sum_n_AIC %>% 
	unnest(contents %>% map(tidy, .id = "number")) %>% 
	# filter(well_plate %in% c(mismatches$well_plate)) %>% 
	filter(best.model == "gr.lagsat", term == "B2") %>% 
	rename(cutoff_point = estimate) %>% 
	select(well_plate, cutoff_point, best.model)

exponential_lag_models <- growth_sum_n_AIC %>% 
	unnest(contents %>% map(tidy, .id = "number")) %>% 
	filter(!best.model %in% c("gr.sat", "gr.lagsat")) %>% 
	# filter(!well_plate %in% c(mismatches$well_plate)) %>%
	distinct(well_plate, .keep_all = TRUE) %>% 
	mutate(cutoff_point = 20) %>% 
	select(well_plate, cutoff_point, best.model)

all_cutoffs <- bind_rows(sat_models, lag_sat_models, exponential_lag_models)


#' now trim the time series of nitrate and force fit an exponential model

nitrate_with_cutoffs <- left_join(nitrate, all_cutoffs, by = "well_plate") %>%
	filter(days < cutoff_point)

nitrate_exponential_growth_rates <- nitrate_with_cutoffs %>%
	filter(population != "COMBO") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	group_by(nitrate_concentration, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 

#### fit Monods

#' Fit the Monod model to the growth rates
monod_fits_cutoff <- nitrate_exponential_growth_rates %>% 
	# filter(!well_plate %in% c(bad_fits$well_plate)) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))

bs_split_cutoff <- monod_fits_cutoff %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_cutoff <- bs_split_cutoff %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population")

all_preds_cutoff2 <- left_join(all_preds_cutoff, treatments, by = c("population")) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)
all_growth_cutoff2 <- left_join(nitrate_exponential_growth_rates, treatments) ## changed this to the 'cut' version. can switch back later



#+ fig.width = 12, fig.height = 8
all_growth_cutoff2 %>% 
	# mutate(estimate = mu) %>%
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	# filter(treatment == "N", ancestor_id == "anc3") %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	# geom_point(aes(x = nitrate_concentration.x, y = estimate), data = filter(nitrate_eyeball_exp, population == 27), color = "blue") +
	geom_line(data=all_preds_cutoff2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	geom_hline(yintercept = 0.1, linetype = "dotted") +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")
ggsave(here("figures", "nitrate-monod-curves-cutoff.png"), width = 12, height =8)



monod_wide_cutoff <-  monod_fits_cutoff %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
#' Find R*
rstars_cutoff <- monod_wide_cutoff %>% 
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

rstars2_cutoff <- left_join(rstars_cutoff, treatments, by = "population") %>%
	distinct(population, ks, umax, .keep_all = TRUE)

#+ fig.width = 6, fig.height = 4
rstars2_cutoff %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2_cutoff, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave(here("figures", "nitrate-r-star-cutoff.png"), width = 6, height = 4)


# AICc growthTools --------------------------------------------------------

#' Fit the Monod model to the growth rates
monod_fits_AICc <- growth_sum_n_AICc %>% 
	left_join(., key) %>% 
	rename(estimate = mu) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))

bs_split_AICc <- monod_fits_AICc %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_AICc <- bs_split_AICc %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(predict_monod, .id = "population")

all_preds_AICc2 <- left_join(all_preds_AICc, treatments, by = c("population")) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)
all_growth_AICc2 <- left_join(nitrate_exponential_growth_rates, treatments) ## changed this to the 'cut' version. can switch back later



#+ fig.width = 12, fig.height = 8
growth_sum_n_AICc %>% 
	mutate(estimate = mu) %>%
	left_join(., key) %>% 
	left_join(., treatments) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	# filter(treatment == "N", ancestor_id == "anc3") %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	# geom_point(aes(x = nitrate_concentration.x, y = estimate), data = filter(nitrate_eyeball_exp, population == 27), color = "blue") +
	geom_line(data=all_preds_AICc2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	geom_hline(yintercept = 0.1, linetype = "dotted") +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")
ggsave(here("figures", "nitrate-monod-curves-AICc.png"), width = 12, height =8)



monod_wide_AICc <-  monod_fits_AICc %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
#' Find R*
rstars_AICc <- monod_wide_AICc %>% 
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

rstars2_AICc <- left_join(rstars_AICc, treatments, by = "population") %>%
	distinct(population, ks, umax, .keep_all = TRUE)

#+ fig.width = 6, fig.height = 4
rstars2_AICc %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2_AICc, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave(here("figures", "nitrate-r-star-AICc.png"), width = 6, height = 4)


### compare all the growth rates

gr_AICc <- growth_sum_n_AICc %>% 
	mutate(method = "AICc") %>% 
	rename(estimate = mu) %>% 
	select(method, estimate, well_plate)

gr_exp_2day <- growth_rates %>% 
	mutate(method = "exponential-2-day")%>% 
	select(method, estimate, well_plate)

gr_cutoff <- nitrate_exponential_growth_rates %>% 
	mutate(method = "cutoff") %>% 
	select(method, estimate, well_plate)

gr_logistic <- growth_rates_logistic %>% 
	mutate(method = "logistic") %>% 
	select(method, estimate, well_plate)

gr_sequential <- exp_fits_top %>% 
	mutate(method = "sequential") %>% 
	select(method, estimate, well_plate)

gr_AIC <- growth_sum_n_AIC %>% 
	mutate(method = "AIC") %>% 
	rename(estimate = mu) %>% 
	select(method, estimate, well_plate)


all_gr <- bind_rows(gr_AIC, gr_AICc, gr_cutoff, gr_exp_2day, gr_logistic, gr_sequential)

# write_csv(all_gr, "data-processed/all_growth_rate_estimates.csv")

all_gr <- read_csv("data-processed/all_growth_rate_estimates.csv") %>% 
	left_join(., key) %>% 
	left_join(., treatments)

#+ fig.width = 14, fig.height = 8
all_gr %>% 
	left_join(., key, by = "well_plate") %>% 
	left_join(., treatments) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate, color = method)) + geom_point(shape = 1) +
	facet_grid(ancestor_id ~ treatment)
ggsave(here("figures", "growth_rate_comparison.png"), width = 12, height =8)


monod_fits_all <- all_gr %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(method, population, treatment, ancestor_id) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))


monod_wide <- monod_fits_all %>% 
	select(method, population, term, estimate, treatment, ancestor_id) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
#' Find R*
rstars_all <- monod_wide %>% 
	mutate(rstar_solve = ks*m/(umax-m))

rstars_all %>% 
	filter(grepl("AIC", method)) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = AIC, y = AICc, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AICc") + xlab("R star via AIC")
ggsave("figures/Rstar_AIC_AICc.png", width = 6, height = 4)

unique(rstars_all$method)

rstars_all %>% 
	filter(method %in% c("cutoff", "sequential")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = cutoff, y = sequential, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via sequential") + xlab("R star via cutoff")
ggsave("figures/Rstar_cutoff_sequential.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("cutoff", "logistic")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = cutoff, y = logistic, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via logistic") + xlab("R star via cutoff")
ggsave("figures/Rstar_cutoff_logistic.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("cutoff", "AIC")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = cutoff, y = AIC, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AIC") + xlab("R star via cutoff")
ggsave("figures/Rstar_cutoff_AIC.png", width = 6, height = 4)


rstars_all %>% 
	filter(method %in% c("cutoff", "AICc")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = cutoff, y = AICc, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AICc") + xlab("R star via cutoff")
ggsave("figures/Rstar_cutoff_AICc.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("cutoff", "exponential-2-day")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = cutoff, y = `exponential-2-day`, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via exponential-2-day") + xlab("R star via cutoff")
ggsave("figures/Rstar_cutoff_exponential-2-day.png", width = 6, height = 4)


### logistic comparisons
rstars_all %>% 
	filter(method %in% c("logistic", "sequential")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = logistic, y = sequential, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via sequential") + xlab("R star via logistic")
ggsave("figures/Rstar_logistic_sequential.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("logistic", "logistic")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = logistic, y = logistic, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via logistic") + xlab("R star via logistic")
ggsave("figures/Rstar_logistic_logistic.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("logistic", "AIC")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = logistic, y = AIC, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AIC") + xlab("R star via logistic")
ggsave("figures/Rstar_logistic_AIC.png", width = 6, height = 4)


rstars_all %>% 
	filter(method %in% c("logistic", "AICc")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = logistic, y = AICc, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AICc") + xlab("R star via logistic")
ggsave("figures/Rstar_logistic_AICc.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("logistic", "exponential-2-day")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = logistic, y = `exponential-2-day`, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via exponential-2-day") + xlab("R star via logistic")
ggsave("figures/Rstar_logistic_exponential-2-day.png", width = 6, height = 4)


### sequential comparisons
rstars_all %>% 
	filter(method %in% c("sequential", "cutoff")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = sequential, y = cutoff, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via cutoff") + xlab("R star via sequential")
ggsave("figures/Rstar_cutoff_sequential.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("sequential", "sequential")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = sequential, y = sequential, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via sequential") + xlab("R star via sequential")
ggsave("figures/Rstar_sequential_sequential.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("sequential", "AIC")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = sequential, y = AIC, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AIC") + xlab("R star via sequential")
ggsave("figures/Rstar_sequential_AIC.png", width = 6, height = 4)


rstars_all %>% 
	filter(method %in% c("sequential", "AICc")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = sequential, y = AICc, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via AICc") + xlab("R star via sequential")
ggsave("figures/Rstar_sequential_AICc.png", width = 6, height = 4)

rstars_all %>% 
	filter(method %in% c("sequential", "exponential-2-day")) %>%
	select(population, method, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ggplot(aes(x = sequential, y = `exponential-2-day`, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("R star via exponential-2-day") + xlab("R star via sequential")
ggsave("figures/Rstar_sequential_exponential-2-day.png", width = 6, height = 4)


sub_plot <- rstars_all %>% 
	select(method, population, rstar_solve) %>% 
	spread(key = method, value = rstar_solve) %>% 
	ungroup() 

ggplot2::plotmatrix(sub_plot[,4:9])

my_line <- function(x,y,...){
	points(x,y,...)
	abline(a = 0,b = 1,...)
}
pairs(sub_plot[1:6], pch = 21, lower.panel = my_line)

library(GGally)
ggpairs(sub_plot, columns = 4:9, upper = "blank", ggplot2::aes(colour=treatment), ggplot2::geom_abline(slope = 1, intercept = 0)) 


?ggpairs
library(ggplot2)
library(gridExtra)
library(GGally)

my_fn <- function(data, mapping, ...){
	p <- ggplot(data = data, mapping = mapping) + 
		geom_point() + 
		geom_abline(slope = 1, intercept = 0, color="black", ...)
	p
}

ggpairs(sub_plot,columns = 4:9, ggplot2::aes(colour=treatment), diag = "blank", upper = "blank", lower = list(continuous = my_fn)) +
	theme(legend.position = "bottom")
ggsave("figures/parwise_combos.png", width = 10, height = 8)


