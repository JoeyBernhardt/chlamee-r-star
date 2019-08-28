
#' ---
#' title: "Chlamee nitrate analysis"
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


#' Read in data
treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629")
nitrate <- read_csv(here("data-processed", "nitrate-abundances-processed.csv"))


#' this is the step that gets us the growth rate estimates
growth_rates_n_AICc <- nitrate %>%
	filter(population != "COMBO") %>% 
	mutate(ln.fluor = log(RFU)) %>% 
	group_by(well_plate) %>% 
	do(grs = get.growth.rate(x = .$days, y = .$ln.fluor,id = .$well_plate, plot.best.Q = F)) 

#' Get growth rates via AIC
growth_rates_n_AIC <- nitrate %>%
	filter(population != "COMBO") %>% 
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


## now try something different, and pull out all the models that were best fit by gr.sat, gr.lagsat

growth_sum_n_AIC_saturated <- growth_rates_n_AIC %>%
	summarise(well_plate, mu = grs$best.slope, n_obs = grs$best.model.slope.n,
			  slope_r2 = grs$best.model.slope.r2,
			  best.model_r2 = grs$best.model.rsqr,
			  best.model = grs$best.model, best.se = grs$best.se,
			  contents = grs$best.model.contents) %>% 
	filter(best.model %in% c("gr.sat", "gr.lagsat"))

#' 1099 out of 1480 are fit differently with AIC and AICc and picked by AIC as a lagsat or sat.
#' For these we will pull out the points before the saturated phase, and fit an exponential model.
all_models %>% 
	filter(best.model.x != best.model.y) %>% 
	filter(best.model.y %in% c("gr.lagsat","gr.sat")) %>% 
	tally() %>% 
	knitr::kable()


mismatches <- all_models %>% 
	filter(best.model.x != best.model.y) %>% 
	filter(best.model.y %in% c("gr.lagsat","gr.sat"))

key <- nitrate %>% 
	select(well_plate, population, nitrate_concentration) %>% 
	distinct(well_plate, .keep_all = TRUE)

AICc_growth_rates <- growth_sum_n_AICc %>% 
	select(-contents) %>% 
	mutate(IC_method = "AICc")

AICc_growth_rates2 <- left_join(AICc_growth_rates, key, by = "well_plate")
# write_csv(AICc_growth_rates2, here("data-processed", "nitrate_exp_growth_w_growthtools_AICc.csv"))

exp_params <- growth_sum_n_AICc %>% 
	unnest(contents %>% map(tidy, .id = "number"))  ## pull out the slopes and intercepts etc.

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

all_cutoffs %>% 
	filter(cutoff_point < 20) %>% View

nitrate_cutoffs <- left_join(nitrate, all_cutoffs, by = "well_plate")
nitrate_cutoffs %>% 
	filter(cutoff_point < 2, nitrate_concentration < 10) %>% 
	ggplot(aes(x = days, y = RFU, color = best.model)) + geom_point() + geom_point(aes(x = cutoff_point, y = 0), size = 3, color = "red") +
	facet_wrap( ~ well_plate)


#' now trim the time series of nitrate and force fit an exponential model

nitrate_with_cutoffs <- left_join(nitrate, all_cutoffs, by = "well_plate") %>%
	filter(days < cutoff_point)


### fit an exponential model

nitrate_exponential_growth_rates <- nitrate_with_cutoffs %>%
	filter(population != "COMBO") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	group_by(nitrate_concentration, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 

library(nls.multstart)




# try again exponential ---------------------------------------------------

ldata_n0 <- nitrate_with_cutoffs %>% 
	group_by(well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	ungroup()


fits_many_n0 <- ldata_n0 %>% 
	group_by(population, well_plate) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ N0 * exp(r*days),
												  data = .x,
												  iter = 500,
												  start_lower = c(r = 0.2),
												  start_upper = c(r = 1),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(r = 0),
												  upper = c(r = 5),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))


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


CI <- fits_many %>% 
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
					   	data.frame() %>% 
					   	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., well_plate) %>%
	mutate(., term = c('r')) %>%
	ungroup()

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

write_csv(params, here("data-processed", "exponential_params_cutoff_approach.csv"))

# get predictions
preds <- fits_many %>%
	unnest(fit %>% map(augment))

new_preds <- ldata_n0 %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))


# max and min for each curve
max_min <- group_by(ldata_n0, well_plate) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

# create new predictions
preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds)) %>%
	merge(., max_min, by = 'well_plate') %>%
	group_by(., well_plate) %>%
	filter(., days > unique(min_days) & days < unique(max_days)) %>%
	rename(., RFU = .fitted) %>%
	ungroup()

key <- ldata_n0 %>% 
	select(well_plate, population, nitrate_concentration) %>% 
	distinct(well_plate, .keep_all = TRUE)

preds3 <- left_join(preds2, key, by = c("well_plate", "population"))


ggplot() +
	geom_point(aes(days, RFU, color = factor(nitrate_concentration)), size = 2, data = ldata_n0) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(nitrate_concentration)), data = preds3) +
	facet_wrap( ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave("figures/nitrate_rfus_exponential_n0.pdf", width = 40, height = 35)


# Again -------------------------------------------------------------------




nitrate_with_cutoffs2 <- nitrate_with_cutoffs %>%
	filter(population != "COMBO") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) 


exponential_fits <- nitrate_with_cutoffs2 %>%
		group_by(nitrate_concentration, population, well_plate) %>%
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ N0 * exp(r*days),
												  data = .x,
												  iter = 500,
												  start_lower = c(r = 0.2),
												  start_upper = c(r = 1),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(r = 0),
												  upper = c(r = 5),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))
fits_many <- exponential_fits

params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 


info <- fits_many %>%
	unnest(fit %>% map(glance))

preds <- fits_many %>%
	unnest(fit %>% map(augment))

new_preds <-  nitrate_with_cutoffs2 %>% 
	distinct(RFU, nitrate_concentration, population, well_plate, .keep_all = TRUE) %>% 
	group_by(well_plate, population) %>% 
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 10), stringsAsFactors = FALSE))

ngrowth <- nitrate_with_cutoffs2 %>% 
	distinct(RFU, nitrate_concentration, population, days, .keep_all = TRUE) 

max_min <- dplyr::group_by(ngrowth, population, well_plate) %>%
	summarise(., min_days = min(days), max_days = max(days)) %>%
	ungroup()

# create new predictions
preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds)) %>%
	merge(., max_min, by = 'well_plate') %>%
	group_by(., well_plate) %>%
	filter(., days > unique(min_days) & days < unique(max_days)) %>%
	rename(., RFU = .fitted) %>%
	ungroup() 


ggplot() +
	geom_point(aes(days, RFU), size = 2, data = ngrowth) +
	geom_line(aes(days, RFU_fitted, group = well_plate), data = preds2) +
	facet_wrap( ~ population) +
	ylab('RFU') +
	xlab('Days')


nitrate_exponential_growth_rates_fitted <- nitrate_with_cutoffs %>%
	filter(population != "COMBO") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	ungroup() %>% 
	group_by(population, well_plate) %>%
	do(augment(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)), newdata = new_preds)) %>% 
	ungroup() 

new_preds <-  nitrate_exponential_growth_rates %>% 
	distinct(estimate, nitrate_concentration, population, .keep_all = TRUE) %>% 
	group_by(population, well_plate) %>% 
	do(., data.frame(nitrate_concentration = seq(min(.$nitrate_concentration), max(.$nitrate_concentration), length.out = 150), stringsAsFactors = FALSE))


nitrate_exponential_growth_rates_fitted %>% 
	ggplot(aes(x = days, y = .fitted, group = well_plate, color = nitrate_concentration)) + geom_line() + 
	geom_point(aes(x = days, y = RFU, color = nitrate_concentration), data = nitrate_exponential_growth_rates_fitted) +
	facet_wrap( ~ population) + scale_color_viridis_c()

nitrate_exp <- nitrate %>% 
	mutate(exponential = case_when(days < 2 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes")

growth_rates_n_AICcut <- nitrate_with_cutoffs %>% 
	# filter(cutoff_point < 20) %>% 
	# nitrate_with_cutoffs %>%
	filter(population != "COMBO") %>% 
	mutate(ln.fluor = log(RFU)) %>% 
	group_by(well_plate) %>% 
	do(grs = get.growth.rate(x = .$days, y = .$ln.fluor,id = .$well_plate, plot.best.Q = F, methods = "linear")) 


nexp <- nitrate_exponential_growth_rates %>% 
	select(estimate, population, well_plate) %>% 
	rename(exp_growth = estimate)
grtools_exp <- growth_rates_n_AICcut %>% 
	summarise(well_plate, estimate = grs$best.slope, n_obs = grs$best.model.slope.n,
			  slope_r2 = grs$best.model.slope.r2,
			  best.model_r2 = grs$best.model.rsqr,
			  best.model = grs$best.model, best.se = grs$best.se,
			  contents = grs$best.model.contents) %>% 
	select(estimate, well_plate) %>% 
	rename(gt_growth = estimate)

all_exponential_models <- left_join(nexp, grtools_exp)

all_exponential_models %>% 
	ggplot(aes(x = gt_growth, y = exp_growth)) + geom_point() +
	geom_abline(slope = 1, intercept = 0)


exponential_lag_models_original <- growth_rates_n_AICc %>%
	filter(!well_plate %in% c(mismatches$well_plate)) 

# all_growth_rates_cut <- bind_rows(growth_rates_n_AICcut, exponential_lag_models_original)

growth_sum_n_AICcut <- growth_rates_n_AICcut %>%
	summarise(well_plate, mu = grs$best.slope, n_obs = grs$best.model.slope.n,
			  slope_r2 = grs$best.model.slope.r2,
			  best.model_r2 = grs$best.model.rsqr,
			  best.model = grs$best.model, best.se = grs$best.se,
			  contents = grs$best.model.contents) 


mods <- nitrate_exponential_growth_rates
exp_params_cut <- growth_sum_n_AICcut %>% 
	unnest(contents %>% map(tidy, .id = "number"))

exp_params_aug_cut <- growth_sum_n_AICcut %>% 
	unnest(contents %>% map(augment, .id = "number")) %>% ### now add the fitted values
	rename(days = x)

exp_wide_cut <- exp_params_cut %>% 
	spread(key = term, value = estimate)

all_preds <- left_join(exp_params_aug_cut, exp_wide_cut)
all_preds2 <- left_join(all_preds, key, by = "well_plate") 
all_preds_exp <- left_join(all_preds, key, by = "well_plate") ### this is now the version with the cutoff exponentials
# %>% 
# 	group_by(well_plate) %>% 
# 	mutate(B1 = mean(B1, na.rm = TRUE)) %>% 
# 	mutate(B2 = mean(B2, na.rm = TRUE)) %>% ### here we do some wrangling to get the colour coding right for our plots
# 	mutate(exponential = case_when(best.model == "gr" ~ "yes",
# 								   best.model == "gr.sat" & days < B1 ~ "yes", 
# 								   best.model == "gr.lag" & days < B1 ~ "yes", 
# 								   best.model == "gr.lagsat" & days < B2 & days > B1 ~ "yes",
# 								   TRUE ~ "no"))


# all_preds2 %>%
# 	ggplot(aes(x = days, y = .fitted, group = well_plate, color = best.model)) + geom_line() +
# 	geom_point(aes(x = days, y = y)) +
# 	facet_grid( ~ well_plate) + ylab("Ln(RFU)") +xlab("Days")
ggsave("figures/all_exponential_nitrate_wells.png", width = 45, height = 20)


### this doesn't look like it's working, because it's still picking up time points that are clearly in the sat phase.
all_preds_exp %>%
	ggplot(aes(x = days, y = .fitted, group = well_plate, color = best.model)) + geom_line() +
	geom_point(aes(x = days, y = y)) +
	facet_grid(nitrate_concentration ~ population) + ylab("Ln(RFU)") +xlab("Days")
ggsave("figures/all_exponential_nitrate_exp_sat_cut.png", width = 45, height = 20)

#+ fig.width=45, fig.height=20
all_preds2 %>%
	ggplot(aes(x = days, y = .fitted, group = well_plate, color = best.model)) + geom_line() +
	geom_point(aes(x = days, y = y, shape = exponential)) +
	facet_grid(nitrate_concentration ~ population) + ylab("Ln(RFU)") +xlab("Days")

AICcut_growth_rates2 <- left_join(growth_sum_n_AICcut, key, by = "well_plate")

AICcut_growth_rates2 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	rename(estimate = mu) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate)) + geom_point() + facet_wrap(~population)


# Monod fits --------------------------------------------------------------


growth_rates <- AICcut_growth_rates2
growth_rates <- left_join(params, key, by = c("well_plate", "population")) %>% 
	select(population, well_plate, estimate) %>% 
	rename(exp_mult = estimate)
growth_rates_exp <- left_join(nitrate_exponential_growth_rates, key) %>% 
	# select(population, well_plate, estimate) %>% 
	rename(exp = estimate)

nitrate_exp <- nitrate %>% 
	mutate(exponential = case_when(days < 2 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) 

growth_rates <- nitrate_exp %>%
	filter(population != "COMBO") %>% 
	group_by(nitrate_concentration, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


growth_rates <- exp_fits_top
growth_rates <- exp_fits_log
growth_rates <- exp_fits_all
growth_rates <- growth_rates_exp
growth_rates <- nitrate_eyeball_exp %>% 
	rename(nitrate_concentration = nitrate_concentration.x)

growth_rates <- read_csv(here("data-processed", "nitrate_exp_growth_w_growthtools_AIC.csv")) %>% 
	filter(population != "COMBO") %>% 
	rename(estimate = mu)


nitrate_eyeball_exp <- read_csv(here("data-processed", "nitrate_exp_growth_eyeball.csv")) %>% 
	distinct(population, nitrate_concentration.x, well_plate.x, .keep_all = TRUE)# left_join(growth_rates, growth_rates_exp) %>% ok so these are the same, which is good.



#' Fit the Monod model to the growth rates
monod_fits <- growth_rates %>% 
	# AICc_growth_rates2 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	# rename(estimate = mu) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))


#' get the fitted values
prediction_function <- function(df) {
	
	monodcurve<-function(x){
		growth_rate<- (df$umax[[1]] * (x / (df$ks[[1]] +x)))
		growth_rate}
	
	pred <- function(x) {
		y <- monodcurve(x)
	}
	
	x <- seq(0, 1000, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(nitrate_concentration.x = x, 
			   growth_rate = preds)
}

bs_split <- monod_fits %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_n <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(prediction_function, .id = "population")


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


monod_wide <-  monod_fits %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve
monod_curve_mortality <- function(nitrate_concentration, umax, ks){
	res <- (umax* (nitrate_concentration / (ks+ nitrate_concentration))) - 0.1
	res
}	

#' Find R*
rstars <- monod_wide %>% 
	mutate(rstar = uniroot.all(function(x) monod_curve_mortality(x, umax, ks), c(0.0, 50))) %>% ## numerical
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

rstars2 <- left_join(rstars, treatments, by = "population") %>%
	distinct(population, ks, umax, .keep_all = TRUE)

#+ fig.width = 6, fig.height = 4
rstars2 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar), y = rstar, color = ancestor_id), size = 2, data = rstars2, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")
ggsave("figures/r-star-means-exp-cutoff.png", width = 6, height = 4)
ggsave("figures/r-star-means-exp-max-r.png", width = 6, height = 4)
ggsave("figures/r-star-means-exp-max-r-log.png", width = 6, height = 4)



# ok now try one more thing  ----------------------------------------------

### fit the exponential model to different time points, and find when the growth rate declines.


nitrate_exponential_growth_rates <- nitrate_with_cutoffs %>% 
	filter(population != "COMBO") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	group_by(nitrate_concentration, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


n2 <- nitrate %>% 
	filter(population != "COMBO") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) 


fitting_window <- function(x) {
	
	growth_rates <- n2 %>% 
	top_n(n = -x, wt = days) %>% 
		group_by(nitrate_concentration, population, well_plate) %>%
		do(tidy(nls(RFU ~ N0 * exp(r*days),
					data= .,  start=list(r=0.01),
					control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}


fitting_window_log <- function(x) {
	
	growth_rates <- n2 %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(nitrate_concentration, population, well_plate) %>%
		do(tidy(lm(log(RFU) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

windows <- seq(4,7, by = 1)

multi_fits <- windows %>% 
	map_df(fitting_window, .id = "iteration")
multi_fits_log <- windows %>% 
	map_df(fitting_window_log, .id = "iteration")


multi_fits %>% 
	ggplot(aes(x = number_of_points, y = estimate, group = well_plate)) + geom_point() + geom_line() +
	facet_wrap( ~ nitrate_concentration)

multi_fits_log %>% 
	filter(term == "days") %>% 
	ggplot(aes(x = number_of_points, y = estimate, group = well_plate)) + geom_point() + geom_line() +
	facet_wrap( ~ nitrate_concentration)

exp_fits_top <- multi_fits %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate) 

exp_fits_log <- multi_fits_log %>% 
	filter(term == "days") %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate) %>% 
	mutate(well_plate_points = paste(well_plate, number_of_points, sep = "_")) %>% 
	filter(nitrate_concentration > 40)

exp_fits_log_less_40 <- multi_fits_log %>% 
	filter(term == "days", nitrate_concentration <= 40) %>% 
	filter(number_of_points == 3) %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate) %>% 
	mutate(well_plate_points = paste(well_plate, number_of_points, sep = "_"))


exp_fits_all <- bind_rows(exp_fits_log, exp_fits_log_less_40)

fitting_window_log_augment <- function(x) {
	
	growth_rates <- n2 %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(nitrate_concentration, population, well_plate) %>%
		do(augment(lm(log(RFU) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

windows <- seq(3,7, by = 1)


multi_fits_log_augment <- windows %>% 
	map_df(fitting_window_log_augment, .id = "iteration")


multi_fits_log_augment2 <- left_join(multi_fits_log_augment, treatments) %>% 
	mutate(well_plate_points = paste(well_plate, number_of_points, sep = "_")) %>% 
	filter(well_plate_points %in% c(exp_fits_log$well_plate_points)) 


left_join(nitrate, treatments, by = "population") %>% 
	filter(treatment == "N", ancestor_id == "anc3") %>% 
	ggplot(aes(x = days, y = RFU)) + geom_point() +
	facet_wrap(~ nitrate_concentration, scales = "free")

multi_fits_log_augment2 %>% 
	filter(treatment == "N", ancestor_id == "anc3") %>% 
	ggplot(aes(x = days, y = .fitted, group = well_plate, color = factor(number_of_points))) + geom_line() +
	geom_point(aes(x = days, y = log.RFU., color = factor(number_of_points)), data = filter(multi_fits_log_augment2, treatment == "N", ancestor_id == "anc3")) + ylab("ln(RFU)") +
	facet_grid(population ~ nitrate_concentration) 
ggsave("figures/N_anc3_trajectories_exp.png", width = 20, height = 20)
ggsave("figures/log_exponential_max_r_trajectories_exp_nitrate.png", width = 30, height = 40)
