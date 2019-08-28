#' ---
#' title: "Chlamee light analysis"
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
library(nls.multstart)
library(nlstools)


#' Read in data
treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(population != "cc1629")


light_data <- read_csv(here("data-processed", "light-rstar-rfus-time.csv"))

ldata <- light_data %>% 
	mutate(percentage = ifelse(percentage == "0.5-0.7", "0.6", percentage)) %>% 
	mutate(percentage = as.numeric(percentage)) %>% 
	mutate(percentage = percentage/100) %>% 
	mutate(light = percentage*250)


#' Estimate growth rates, using piecewise regression (either exponential, lagged growth, saturated or lagged then saturated). Here we are choosing between moddels using AICc.
growth_rates_n_AICc <- ldata %>%
	filter(population != "COMBO") %>% 
	mutate(ln.fluor = log(RFU)) %>% 
	group_by(well_plate) %>% 
	do(grs = get.growth.rate(x = .$days, y = .$ln.fluor,id = .$well_plate, plot.best.Q = F)) 


#' Pull out the things we want
growth_sum_n_AICc <- growth_rates_n_AICc %>%
	summarise(well_plate, mu = grs$best.slope, n_obs = grs$best.model.slope.n,
			  slope_r2 = grs$best.model.slope.r2,
			  best.model_r2 = grs$best.model.rsqr,
			  best.model = grs$best.model, best.se = grs$best.se,
			  contents = grs$best.model.contents) 

key <- ldata %>% 
	select(well_plate, population, light) %>% 
	distinct(well_plate, .keep_all = TRUE)

AICc_growth_rates <- growth_sum_n_AICc %>% 
	select(-contents) %>% 
	mutate(IC_method = "AICc")

AICc_growth_rates2 <- left_join(AICc_growth_rates, key, by = "well_plate")
# write_csv(AICc_growth_rates2, here("data-processed", "nitrate_exp_growth_w_growthtools_AICc.csv"))

exp_params <- growth_sum_n_AICc %>% 
	unnest(contents %>% map(tidy, .id = "number"))  ## pull out the slopes and intercepts etc.

exp_params_aug <- growth_sum_n_AICc %>% 
	unnest(contents %>% map(augment, .id = "number")) %>% ### now add the fitted values
	rename(days = x)

exp_wide <- exp_params %>% 
	spread(key = term, value = estimate)

all_preds <- left_join(exp_params_aug, exp_wide)
all_preds2 <- left_join(all_preds, key, by = "well_plate") %>% 
	group_by(well_plate) %>% 
	mutate(B1 = mean(B1, na.rm = TRUE)) %>% 
	mutate(B2 = mean(B2, na.rm = TRUE)) %>% ### here we do some wrangling to get the colour coding right for our plots
	mutate(exponential = case_when(best.model == "gr" ~ "yes",
								   best.model == "gr.sat" & days < B1 ~ "yes", 
								   best.model == "gr.lag" & days < B1 ~ "yes", 
								   best.model == "gr.lagsat" & days < B2 & days > B1 ~ "yes",
								   TRUE ~ "no"))


#' Plot all the time series, here colour-coded by which model was the best fit, and the shape of the points corresponds to whether the point was in the exponential phase
#+ fig.width=45, fig.height=20
all_preds2 %>%
	ggplot(aes(x = days, y = .fitted, group = well_plate, color = best.model)) + geom_line() +
	geom_point(aes(x = days, y = y, shape = exponential)) +
	facet_grid(light ~ population) + ylab("Ln(RFU)") +xlab("Days")
ggsave(here("figures", "light_time_series_growth_tools.png"), width = 45, height = 20)



#' Fit Eilers-Peeters model

epcurve <- function(light, alpha, eopt, ps, m){
	res <-  (light/((1/(alpha * 
							eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))) + h
	res
}

AICc_growth_rates3 <- left_join(AICc_growth_rates2, treatments, by = "population")

ep_fits_i_m <- AICc_growth_rates3 %>%
	rename(r_estimate = mu) %>% 
	distinct(r_estimate, light, population, .keep_all = TRUE) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(r_estimate ~ (light/((1/(alpha * 
																		   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha)))-0.1,
												  data = .x,
												  iter = 500,
												  start_lower = c(alpha = 0.2, eopt = 100, ps = 1),
												  start_upper = c(alpha = 0.4, eopt = 200, ps = 2),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(alpha = 0, eopt = 80, ps = 0.1),
												  upper = c(alpha = 4, eopt = 400, ps = 5),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))
fits_many <- ep_fits_i_m

params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 


info <- fits_many %>%
	unnest(fit %>% map(glance))

ep_wide <- params %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate, 2:3) 


CI <- fits_many %>% 
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
					   	data.frame() %>% 
					   	dplyr::rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., population) %>%
	mutate(., term = c('alpha', 'eopt', 'ps')) %>%
	ungroup()

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))
write_csv(params, here("data-processed", "EP-with-mortality-parameters.csv"))


# get predictions
preds <- fits_many %>%
	unnest(fit %>% map(augment))

new_preds <-  AICc_growth_rates3 %>% 
	rename(estimate = mu) %>% 
	distinct(estimate, light, population, .keep_all = TRUE) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
	do(., data.frame(light = seq(min(.$light), max(.$light), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve

p_growth3 <- AICc_growth_rates3 %>% 
	rename(estimate = mu) %>% 
	distinct(estimate, light, population, .keep_all = TRUE) 

max_min <- dplyr::group_by(p_growth3, population, treatment, ancestor_id) %>%
	summarise(., min_light = min(light), max_light = max(light)) %>%
	ungroup()

# create new predictions
preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds)) %>% 
	merge(., max_min) %>% 
	group_by(., population) %>% 
	filter(., light > unique(min_light) & light < unique(max_light)) %>%
	dplyr::rename(., r_estimate = .fitted) %>%
	ungroup()  
	
#' Visualize the fits of the EP model to the growth rate data
#+ fig.width = 14, fig.height = 8
ggplot() +
	geom_point(aes(light, estimate), size = 2, data = p_growth3) +
	geom_line(aes(light, r_estimate, group = population), data = preds2) +
	facet_grid(ancestor_id ~ treatment, scales = "free") +
	ylab('Growth rate (/day)') +
	xlab('Irradiance')


#' Now find I* from root of EP, set mortality rate (respiratory losses) to 0.1 (as in m = 0.1 for the nitrate analysis) 

epcurve_m<-function(light, alpha, eopt, ps){
	res <-  (light/((1/(alpha * 
							eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha)))-0.1
	res
}

get_rstar_i <- function(df){
	roots <- uniroot.all(function(x) epcurve_m(x, alpha = df$estimate[[1]], eopt = df$estimate[[2]],
											   ps = df$estimate[[3]]), interval = c(0, 10))
	rstar <- roots[[1]]
	rstar2 <- as.data.frame(rstar)
	return(rstar)
}

ep_fits_split <- params %>% 
	split(.$population)


rstars_ep <- ep_fits_split %>% 
	map(get_rstar_i) %>% 
	unlist() %>%
	as_data_frame() %>% 
	dplyr::rename(rstar = value) %>% 
	mutate(population = rownames(.))

rstars2_ep <- left_join(rstars_ep, treatments, by = "population")

#' Plot the resulting I*
#+ fig.width = 6, fig.height = 4
rstars2_ep %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	# geom_boxplot(aes(x = treatment, y = rstar), color= "black", data = filter(rstars2_ep, rstar > 0)) +
	geom_point(aes(x = reorder(treatment, rstar), y = rstar, color = ancestor_id), data =rstars2_ep, size = 2) + ylab("I* (umols/m2/s)") +
	xlab("Selection treatment") 

