

### estimate growth rates
library(tidyverse)
library(cowplot)
library(broom)
library(readxl)
library(janitor)
library(nls.multstart)
library(rootSolve)
library(nlstools)
library(plotrix)

light_data <- read_csv("data-processed/light-rstar-rfus-time.csv")
treatments <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))

ldata <- light_data %>% 
	mutate(percentage = ifelse(percentage == "0.5-0.7", "0.6", percentage)) %>% 
	mutate(percentage = as.numeric(percentage)) %>% 
	mutate(percentage = percentage/100) %>% 
	mutate(light = percentage*250)

unique(light_data$percentage)


ldata %>% 
	ggplot(aes(x = days, y = RFU, color = light, group = well_plate)) + geom_point() + 
	scale_color_viridis_c(name = "Light") + 
	facet_wrap( ~ light, scales = "free_y") + geom_vline(xintercept = 2.5) +
	geom_line()

ldata_exp <- ldata %>% 
	mutate(exponential = case_when(light > 20 & days < 2 ~ "yes",
								   light < 20 & days < 5 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes") %>% 
	group_by(population, light, well_plate) %>% 
	mutate(N0 = RFU[[1]]) 

ldata_exp %>% 
	ggplot(aes(x = days, y = RFU, color = light, group = well_plate)) + geom_point() + 
	scale_color_viridis_c(name = "Light") + 
	geom_line() +
	facet_wrap( ~ light, scales = "free_y") + geom_vline(xintercept = 2.5)


# ldata_exp %>% 
# 	group_by(light, population, well_plate) %>%


growth_rates <- ldata_exp %>%
	filter(population != "COMBO") %>% 
	group_by(light, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


# growth_rates <- ldata_exp %>%
# 	filter(population != "COMBO") %>% 
# 	group_by(light, population, well_plate) %>%
# 	do(tidy(lm(log(RFU + 1) ~ days, data = . ), conf.int = TRUE)) %>% 
# 	filter(term == "days")

library(growthTools)

growth_rates_i <- ldata %>%
	# filter(plate == 1) %>% 
	filter(population != "COMBO") %>% 
	mutate(ln.fluor = log(RFU)) %>% 
	# filter(nitrate_level > 1) %>% 
	group_by(well_plate) %>% 
	do(grs=get.growth.rate(x=.$days, y=.$ln.fluor,id=.$well_plate,plot.best.Q=T,fpath="/Users/joeybernhardt/Documents/Narwani/ChlamEE-R-star/figures/growth_curve_fits/light/"))

growth_sum_i <- growth_rates_i %>%
	summarise(well_plate,mu=grs$best.slope,
			  best.model=grs$best.model,best.se=grs$best.se)


treatments_light <- ldata %>% 
	select(well_plate, population, light)

treatments_light2 <- left_join(treatments_light, treatments, by = "population")

all_growth_i <- left_join(growth_sum_i, treatments_light2, by = "well_plate") 


write_csv(all_growth_i, "data-processed/exponential_growth_light_growth_tools.csv")



growth2 <- left_join(growth_rates, treatments, by = "population") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))


growth2 %>% 
	ggplot(aes(x = light, y = estimate)) + geom_point() +
	facet_grid(treatment ~ ancestor_id) + geom_hline(yintercept = 0)

write_csv(growth2, "data-processed/light-r-star-growth-rates.csv")


library(phytotools)



### 

growth2 <- read_csv("data-processed/light-r-star-growth-rates.csv")

df <- growth2 %>% 
	filter(population == 4)


### define fitting function
fit_ep <- function(df) {
	PI <- fitEP(df$light, df$estimate, normalize = FALSE, lowerlim = c(0, 0, 0), upperlim = c(4, 250, 3),
				fitmethod=c("SANN"))
	alpha <- PI$alpha
	eopt <- PI$eopt
	ps <- PI$ps
	all_params <- bind_rows(alpha, eopt, ps)
	param <- data.frame(parameter = c("alpha", "iopt", "umax"))
	all_params2 <- bind_cols(param, all_params) %>% clean_names()
	return(all_params2)
}

growth_split <- growth2 %>% 
	filter(estimate > 0) %>% ## apparently the EP model is only appropriate for positive growth rates
	split(.$population)


ep_fits <- growth_split %>% 
	map_df(fit_ep, .id = "population")


epcurve<-function(light, alpha, eopt, ps){
	res <-  light/((1/(alpha * 
			   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))
	res
}


# Monod fits --------------------------------------------------------------

monod_fits_i <- growth2 %>% 
	filter(estimate > 0) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
	do(tidy(nls(r_estimate ~ umax* (light / (ks+ light)),
				data= .,  start=list(ks = 1, umax = 2), algorithm="port",
				control = nls.control(maxiter=500, minFactor=1/204800000))))

monod_fits_i_growth_tools <- all_growth_i %>% 
	# filter(estimate > 0) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= mu) %>% 
	do(tidy(nls(r_estimate ~ umax* (light / (ks+ light)),
				data= .,  start=list(ks = 1, umax = 2), algorithm="port",
				control = nls.control(maxiter=500, minFactor=1/204800000))))


monod_fits_i %>% 
	ggplot(aes(x = treatment, y = estimate, color = ancestor_id)) + geom_point(size = 3) +
	geom_boxplot(aes(x = treatment, y = estimate), data = monod_fits_i, color = "black") +
	geom_point(size = 3) +
	facet_wrap( ~ term, scales = "free_y") 
ggsave("figures/monod-params-light-rstar.pdf", width = 8, height = 4)

preds_i <- growth2 %>%
	filter(estimate > 0) %>% 
	group_by(treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
	do(augment(nls(r_estimate ~ umax* (light / (ks+ light)),
				   data= .,  start=list(ks = 1, umax = 2), algorithm="port",
				   control = nls.control(maxiter=500, minFactor=1/204800000)))) %>% 
	unite(col = unique_id, treatment, ancestor_id, remove = FALSE)



growth2 %>% 
	filter(estimate > 0) %>% 
	ggplot(aes(x= light, y= estimate, color = treatment)) +
	geom_line(data= preds_i, aes(x=light, y=.fitted, group = unique_id, color =treatment)) +
	geom_point() +
	# geom_errorbar(aes(ymin=mu-best.se, ymax=mu + best.se), width=.2) + 
	# facet_grid(rows = vars(ancestor_name), cols = vars(selection_treatment)) +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Irradiance (umols/m2/s)") + geom_hline(yintercept = 0)
ggsave("figures/monod_curves_light_r_star.pdf", width = 15, height = 15)



prediction_function <- function(df) {
		ep <-  function(x){
			res <- x/((1/(df$estimate[[1]] * df$estimate[[2]]^2)) * x^2 + (1/df$estimate[[3]] - 
		2/(df$estimate[[1]] * df$estimate[[2]])) * x + (1/df$estimate[[1]]))
	res
	}
	
	
	pred <- function(x) {
		y <- ep(x)
	}
	
	x <- seq(-2, 250, by = 1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(irradiance = x, 
			   growth = preds)
}

ep_fits_split <- ep_fits %>% 
	split(.$population)

all_preds <- ep_fits_split %>% 
	map_df(prediction_function, .id = "population")

all_preds2 <- left_join(all_preds, treatments, by = "population") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))


growth2 %>% 
	filter(estimate > 0) %>% 
	ggplot(aes(x = light, y = estimate, color = treatment)) + geom_point() +
	geom_line(aes(x = irradiance, y = growth, group = population, color = treatment), data = all_preds2) +
	facet_grid(treatment ~ ancestor_id) + ylab("Exponential growth rate (/day)") +
	xlab("Irradiance (umols/m2/s)") + geom_vline(xintercept = 0) + xlim(-1, 1) + ylim(-0.11, 1) +
	geom_hline(yintercept = 0)
ggsave("figures/ep_curves.pdf", width = 15, height = 15)


ep_fits2 <- left_join(ep_fits, treatments, by = "population") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))


ep_fits2 %>% 
	ggplot(aes(x = treatment, y = estimate, color = treatment)) + geom_point() +
	geom_boxplot(aes(x = treatment, y = estimate, color = treatment), data = ep_fits2) +
	geom_point() +
	facet_wrap( ~ parameter, scales = "free_y")
ggsave("figures/ep-params-light-r-star.pdf", width = 14, height = 4)

ep_fits2 %>% 
	ggplot(aes(x = treatment, y = estimate, color = ancestor_id)) + geom_point(size = 3) +
	geom_boxplot(aes(x = treatment, y = estimate), color= "black", data = ep_fits2) +
	geom_point(size = 3) +
	facet_wrap( ~ parameter, scales = "free_y")
ggsave("figures/ep-params-light-r-star-black-box.pdf", width = 14, height = 4)



write_csv(ep_fits2, "data-processed/ep_fits.csv")

ep_fits2 %>% 
	filter(parameter == "iopt") %>% 
	lm(estimate ~ treatment, data = .) %>% summary()

epcurve <-function(light, alpha, iopt, ps){
	res <-  light/((1/(alpha * 
					   	iopt^2)) * light^2 + (1/ps - 2/(alpha * iopt)) * light + (1/alpha))
	res
}

get_rstar_i <- function(df){
	roots <- uniroot.all(function(x) epcurve(x, alpha = df$estimate[[1]], iopt = df$estimate[[2]],
										 ps = df$estimate[[3]]), interval = c(0, 10))
	rstar <- roots[[1]]
	rstar2 <- as.data.frame(rstar)
	return(rstar)
	}

ep_fits_split <- ep_fits2 %>% 
	split(.$population)

df <- ep_fits_split[[3]]

rstars_ep <- ep_fits_split %>% 
	map(get_rstar_i) %>% 
	unlist() %>%
	as_data_frame() %>%
	rename(rstar = value) %>%
	mutate(population = rownames(.))

rstars2_ep <- left_join(rstars_ep, treatments, by = "population")



rstars2_ep %>% 
	# filter(rstar > 0) %>% 
	ggplot(aes(x = treatment, y = rstar, color = ancestor_id)) + geom_point(size = 3) +
	geom_boxplot(aes(x = treatment, y = rstar), color= "black", data = filter(rstars2_ep, rstar > 0)) +
	geom_point(size = 3) + ylab("R* (umols/ms/s)")
ggsave("figures/i-star-with-zeros.png", width = 6, height = 4)


# find R* -----------------------------------------------------------------


monod_curve <- function(light, umax, ks){
	res<- umax* (light / (ks+ light))
	res
}


monod_fits_i_split <- monod_fits_i %>% 
	split(.$population)

df <- monod_fits_i_split[[1]]

get_rstar <- function(df){
	uniroot.all(function(x) monod_curve(x, umax = df$estimate[[2]], ks = df$estimate[[1]]), c(-5, 50))
}



rstars <- monod_fits_i_split %>% 
	map(get_rstar) %>% 
	unlist() %>% 
	as_data_frame() %>%
	rename(rstar = value) %>% 
	mutate(population = rownames(.))

rstars2 <- left_join(rstars, treatments, by = "population")

rstars2 %>% 
	ggplot(aes(x = treatment, y = rstar, color = ancestor_id)) + geom_point(size = 3) +
	geom_boxplot(aes(x = treatment, y = rstar), color= "black", data =rstars2) +
	geom_point(size = 3) + ylab("R* (umols/ms/s)")
ggsave("figures/I-star-monod.pdf", width = 8, height = 4)
	

### find rstar for light from the monod curve fits
m <- 0.1

rstars_i <- monod_fits_i_growth_tools %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	# mutate(rstar = uniroot(function(x) monod_curve(x, umax, ks), c(0.0, 50))) %>% 
	mutate(rstar_solve = ks*m/(umax-m)) 


rstars2i <- left_join(rstars_i, treatments, by = "population")
rstars_i %>% 
	filter(!is.na(rstar_solve)) %>% 
	# filter(population != 3) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point(size = 3) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) +
	ylab("I* (umols/m2/s)") + xlab("Selection treatment") + scale_color_viridis_d() +
	geom_point(aes(x = treatment, y = rstar_solve, color = ancestor_id), data = rstars_i, size = 3)
ggsave("figures/I-star-monod-analytical-growthtools.pdf", width = 8, height = 4)




# fit logistic to light RFUs ---------------------------------------------------------
library(nls.multstart)
library(purrr)


ldata_n0 <- ldata %>% 
	group_by(well_plate) %>% 
	mutate(N0 = min(RFU)) %>% 
	ungroup()
	

fits_many_n0 <- ldata_n0 %>% 
	# filter(light < 27) %>% 
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
												  upper = c(K = 3000, r = 2.5),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))



fits_many_lowK <- ldata %>% 
	filter(light < 27) %>% 
	group_by(population, well_plate) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ K/(1 + (K/N0- 1)*exp(-r*days)),
												  data = .x,
												  iter = 500,
												  start_lower = c(K = 100, N0 = 10, r = 0),
												  start_upper = c(K = 500, N0 = 30, r = 3),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(K = 10, N0 = 0, r = 0),
												  upper = c(K = 1000, N0 = 150, r = 4),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

fits_many_original <- fits_many

fits_many <- fits_many_lowK
fits_many <- fits_many_lowK_n0
fits_many <- fits_many_n0

fits_well_plate <- fits_many %>% 
	select(well_plate)
pops <- left_join(fits_well_plate, key)

info <- fits_many %>%
	unnest(fit %>% map(glance))

# get params
params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 

# get confidence intervals
library(nlstools)

CI <- fits_many %>% 
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
					   	data.frame() %>% 
					   	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., well_plate) %>%
	mutate(., term = c('K', 'r')) %>%
	ungroup()

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits_many %>%
	unnest(fit %>% map(augment))

new_preds <- ldata %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 150), stringsAsFactors = FALSE))


# max and min for each curve
max_min <- group_by(ldata, well_plate) %>%
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

key <- ldata %>% 
	select(well_plate, population, light) %>% 
	distinct(well_plate, .keep_all = TRUE)

preds3 <- left_join(preds2, key, by = c("well_plate", "population"))


ggplot() +
	geom_point(aes(days, RFU, color = light), size = 2, data = ldata) +
	geom_line(aes(days, RFU, group = well_plate, color = light), data = preds3) +
	facet_wrap(light ~ population, scales = "free") +
	scale_color_viridis_c() +
	ylab('RFU') +
	xlab('Days')
ggsave("figures/light_rfus_logistic_n0_light.pdf", width = 40, height = 35)


pop35 <- ldata %>% 
	filter(population == 35)


pred35 <- preds3 %>% 
	filter(population == 35)
ggplot() +
	geom_point(aes(days, RFU, color = factor(light)), size = 2, data = filter(ldata, well_plate == "E03_28")) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(light)), data = filter(preds3, well_plate == "E03_28")) +
	facet_wrap(~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')


unique(ldata$light)
#0.25  12.50  27.50  50.00  82.50 125.00 175.00 250.00   1.50   5.00
## ok so it looks like the light levels above 27.5 inclusive are fit well by the logistic, but not below. 
ggplot() +
	geom_point(aes(days, RFU, color = factor(light)), size = 2, data = filter(ldata, population == "8")) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(light)), data = filter(preds3, population == "8")) +
	facet_wrap(light ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave("figures/logistic_pop8.pdf", width = 15, height = 10)
ggsave("figures/logistic_pop8.png", width = 15, height = 10)

ldata %>% 
	select(light) %>% View


params2 <- left_join(params, key)

params3 <- left_join(params2, treatments) ## still some issues with hitting the upper bounds for r at low lights, but better

write_csv(params3, "data-processed/r-estimates-light-rstar-logistic.csv")


unique(params$population)

params3 %>% 
	filter(population == 8, term == "r") %>% 
	ggplot(aes(x = light, y = estimate, color = factor(light))) + geom_point() +
	facet_wrap( ~ population)
ggsave("figures/pop8_light_growth.pdf", width = 8, height = 6)

params3 %>% 
	filter(term == "r") %>% 
	ggplot(aes(x = light, y = estimate, color = factor(light))) + geom_point() +
	facet_wrap( ~ population) + ylab("Intrinsic rate of increase (/day)") + xlab("Irradiance (umols/ms/s)") +
	scale_color_viridis_d(name = "Light level")
ggsave("figures/light_growth.pdf", width = 15, height = 10)
ggsave("figures/light_growth.png", width = 15, height = 10)


params3 %>% View 
	filter(population == 11, term == "K") %>% View


params3 %>% 
	ggplot(aes(x = treatment, y = estimate, color = light)) + geom_point() +
	facet_wrap(~ term, scales = "free")

params_wide <- params3 %>% 
	select(well_plate, term, estimate) %>% 
	spread(key = term, value = estimate)

params_wide2 <- left_join(params_wide, params3)

params_wide2 %>% 
	ggplot(aes(x = r, y = K, color = factor(light))) + geom_point()
ggsave("figures/r_K_light.png", width = 8, height = 6)




params3 %>% View
	filter(term == "K") %>% 
	ggplot(aes(x = r, y = K)) + geom_point()


ldata %>% 
	filter(light < 50) %>% 
	View



# r-alpha logistic --------------------------------------------------------

p_growth <- read_csv("data-processed/exponential_growth_light_growth_tools.csv")

unique(p_growth$population)
unique(key$population)


p_growth2 <- p_growth %>%
	mutate(method = "piecewise") %>% 
	rename(estimate = mu) %>% 
	select(population, ancestor_id, treatment, estimate, light, method)

unique(l_growth$population)

l_growth <- read.csv("data-processed/r-estimates-light-rstar-logistic.csv") %>%
	filter(term == "r") %>% 
	mutate(method = "logistic") %>% 
	mutate(population = as.character(population)) %>% 
	select(population, ancestor_id, treatment, estimate, light, method)


all_growth <- bind_rows(p_growth2, l_growth) %>% 
	filter(population != "COMBO")

all_growth %>% 
	ggplot(aes(x = light, y = estimate, color = method)) + geom_point() +
	facet_wrap( ~ population) + ylab("Growth rate (/day)") + xlab("Irradiance (umols/ms/s)") +
	scale_color_viridis_d(name = "Method", begin = 0.4, end = 0.9)
ggsave("figures/logistic-piecewise.pdf", width = 15, height = 10)
ggsave("figures/logistic-piecewise.png", width = 15, height = 10)
	


# fit EP from scratch with a mortality term -------------------------------

epcurve<-function(light, alpha, eopt, ps, m){
	res <-  (light/((1/(alpha * 
					   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))) + h
	res
}

ep_fits_i_m <- p_growth2 %>% 
	# filter(estimate > 0) %>%
	distinct(estimate, light, population, .keep_all = TRUE) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
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

ep_fits_i_h <- p_growth2 %>% 
	# filter(estimate > 0) %>%
	distinct(estimate, light, population, .keep_all = TRUE) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(r_estimate ~ (light/((1/(alpha * 
																		   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))) + h,
												  data = .x,
												  iter = 500,
												  start_lower = c(alpha = 0, eopt = 100, ps = 2, h = -2),
												  start_upper = c(alpha = 1, eopt = 200, ps = 4, h = -1),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(alpha = 0, eopt = 80, ps = 0.1, h = -4),
												  upper = c(alpha = 4, eopt = 600, ps = 5, h = 0),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))


ep_fits_i_h <- p_growth2 %>% 
	# filter(estimate > 0) %>%
	distinct(estimate, light, population, .keep_all = TRUE) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(r_estimate ~ (light/((1/(alpha * 
																		   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))) + h,
												  data = .x,
												  iter = 500,
												  start_lower = c(alpha = 0, eopt = 100, ps = 2, h = -2),
												  start_upper = c(alpha = 1, eopt = 200, ps = 4, h = -1),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(alpha = 0, eopt = 80, ps = 0.1, h = -4),
												  upper = c(alpha = 4, eopt = 600, ps = 5, h = 0),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))




epcurve_h<-function(light, alpha, eopt, ps, h){
	res <-  (light/((1/(alpha * 
							eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))) + h
	res
}


data <- p_growth2 %>% 
	filter(population == 11)


## try out mle to see if we get better fits. No dice.
# library(bbmle)
# fit <- mle2(data$estimate ~ dnorm(mean=epcurve_h(data$light, alpha = alpha, eopt = eopt, ps = ps, h = h), sd=s),
# 			start=list(alpha = 0.2, eopt = 200, ps = 1.5, h = -0.1, s=0.2),
# skip.hessian=TRUE, data = data)
# 
# summary(fit)	
# warnings()

# fits_many <- ep_fits_i
fits_many <- ep_fits_i_h

params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 

# fits_many <- ep_fits_i_m

info <- fits_many %>%
	unnest(fit %>% map(glance))

# get params

weird_hs <- params %>% 
	filter(estimate == 0)

weird_pops <- weird_hs$population

ep_wide <- params %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate, 2:3) 


positive_h <- ep_wide %>% 
	filter(h >= 0)


weird_pops <- positive_h$population ### these are the populations that show positive respiratory losses...we could set this at 0.

ep_wide %>% 
	ggplot(aes(x = h, y = alpha)) + geom_point()

ps <- ep_wide %>% 
	select(-population)

cor(ps)

# get confidence intervals


CI <- fits_many %>% 
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(~ confint2(.x) %>%
					   	data.frame() %>% 
					   	dplyr::rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., population) %>%
	mutate(., term = c('alpha', 'eopt', 'ps', 'h')) %>%
	ungroup()

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

write_csv(params, "data-processed/EP-with-mortality-params.csv")
write_csv(params, "data-processed/EP-with-variable-mortality-params.csv")

params <- read_csv("data-processed/EP-with-mortality-params.csv")

cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")


params <- params %>% 
	mutate(treatment = str_replace(treatment, "none", "A")) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",

							  		 							  		 "P", "B", "S", "BS"))) %>% 
	mutate(term = str_replace(term, "ps", "umax")) %>% 
	mutate(term = str_replace(term, "eopt", "Iopt"))
params %>% 
	group_by(treatment, term) %>% 
	
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(aes(x = treatment, y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	geom_point(aes(x = treatment, y = estimate, color = treatment), data = params, size = 2) +
	scale_color_manual(values = c("black", cols_no_anc), name = "Ancestor") +
	facet_wrap( ~ term, scales = "free") + ylab("Parameter value") + xlab("Selection treatment") +
	theme(legend.position = "none") + xlab("")
ggsave("figures/EP-params-ancestor-id.png", width = 8, height = 2.5)

# get predictions
preds <- fits_many %>%
	unnest(fit %>% map(augment))

new_preds <-  p_growth2 %>% 
	# filter(estimate > 0) %>%
	distinct(estimate, light, population, .keep_all = TRUE) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	mutate(r_estimate= estimate) %>% 
	do(., data.frame(light = seq(min(.$light), max(.$light), length.out = 150), stringsAsFactors = FALSE))


# max and min for each curve

p_growth3 <- p_growth2 %>% 
	# filter(estimate > 0) %>%
	distinct(estimate, light, population, .keep_all = TRUE) 

max_min <- dplyr::group_by(p_growth3, population, treatment, ancestor_id) %>%
	summarise(., min_light = min(light), max_light = max(light)) %>%
	ungroup()

# create new predictions
preds2 <- fits_many %>%
	unnest(fit %>% map(augment, newdata = new_preds)) %>% 
	merge(., max_min, by = 'population') %>% 
	group_by(., population) %>% 
	filter(., light > unique(min_light) & light < unique(max_light)) %>%
	dplyr::rename(., r_estimate = .fitted) %>%
	ungroup()



ggplot() +
	geom_point(aes(light, estimate), size = 2, data = p_growth3) +
	geom_line(aes(light, r_estimate, group = population), data = preds2) +
	facet_wrap(~ population, scales = "free") +
	ylab('Growth rate') +
	xlab('Irradiance')
ggsave("figures/fitted_ep_with_mortality.pdf", width = 15, height = 10)
ggsave("figures/fitted_ep_with_mortality.png", width = 15, height = 10)




# Now find I* from root of EP ---------------------------------------------

epcurve_m<-function(light, alpha, eopt, ps){
	res <-  (light/((1/(alpha * 
					   	eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha)))-0.1
	res
}

epcurve_h<-function(light, alpha, eopt, ps, h){
	res <-  (light/((1/(alpha * 
							eopt^2)) * light^2 + (1/ps - 2/(alpha * eopt)) * light + (1/alpha))) + h
	res
}

get_rstar_i <- function(df){
	roots <- uniroot.all(function(x) epcurve_m(x, alpha = df$estimate[[1]], eopt = df$estimate[[2]],
											 ps = df$estimate[[3]]), interval = c(0, 10))
	rstar <- roots[[1]]
	rstar2 <- as.data.frame(rstar)
	return(rstar)
}

df <- ep_fits_split[[7]]

get_rstar_ih <- function(df){
	roots <- uniroot.all(function(x) epcurve_h(x, alpha = df$estimate[[1]], eopt = df$estimate[[2]],
											   ps = df$estimate[[4]], h = df$estimate[[3]]), interval = c(0, 50))
	rstar <- roots[[1]]
	rstar2 <- as.data.frame(rstar)
	return(rstar)
}




ep_wide <- params %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate, 2:3) 


istars <- ep_wide %>% 
	mutate(istar = uniroot.all(function(x) epcurve_h(x, alpha, eopt, h, ps), c(-10, 150))) ## numerical

ep_fits_split <- params %>% 
	split(.$population)

ep_fits_split[[1]]

rstars_ep <- ep_fits_split %>% 
	map(get_rstar_ih) %>% 
	unlist() %>%
	as_data_frame() %>% 
	dplyr::rename(rstar = value) %>% 
	mutate(population = rownames(.))

rstars2_ep <- left_join(rstars_ep, treatments, by = "population")

write_csv(rstars2_ep, "data-processed/I-stars-from-EP-with-m01.csv")
write_csv(rstars2_ep, "data-processed/I-stars-from-EP-with-variable-m-capped-0.csv")

rstars_fixed_m <- read_csv("data-processed/I-stars-from-EP-with-m01.csv")


rstars2_ep %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	# geom_boxplot(aes(x = treatment, y = rstar), color= "black", data = filter(rstars2_ep, rstar > 0)) +
	geom_point(aes(x = reorder(treatment, rstar), y = rstar, color = ancestor_id), data =rstars2_ep, size = 2) + ylab("I* (umols/m2/s)") +
	xlab("Selection treatment") + scale_color_brewer(name = "Ancestor", type = "qual", palette = 4)

ggsave("figures/i_star_ep.pdf", width = 6, height = 4)
ggsave("figures/i_star_ep.png", width = 6, height = 4)



# fixed h = 0.1 figure ----------------------------------------------------


#### fixed h = 0.1
rstars_fixed_m %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	# geom_boxplot(aes(x = treatment, y = rstar), color= "black", data = filter(rstars2_ep, rstar > 0)) +
	geom_point(aes(x = reorder(treatment, rstar), y = rstar, color = ancestor_id), data = rstars_fixed_m, size = 2) + ylab("I* (umols/m2/s)") +
	xlab("Selection treatment") + scale_color_brewer(name = "Ancestor", type = "qual", palette = 4)
ggsave("figures/i_star_ep_m01.png", width = 6, height = 4)


## plot all the EP curves

ggplot() +
	# geom_point(aes(light, estimate), size = 2, data = p_growth3) +
	geom_line(aes(light, r_estimate, group = population, color = treatment.x), data = preds2) +
	# facet_wrap(~ population, scales = "free") +
	ylab('Growth rate') +
	xlab('Irradiance') +
	scale_color_brewer(name = "Selection treatment", type = "qual", palette = 3)


# function valued trait stuff ---------------------------------------------



### PCA -- so far nothing seems to make sense here
library(cowplot)
library(vegan)
library(ggbiplot)



growth_16 <- preds2 %>% 
	select(population, light, r_estimate) %>% 
	filter(light %in% light_10) %>% 
	distinct(population, light, r_estimate) %>% 
	spread(key = light, value = r_estimate, 2:3) %>% 
	select(-population)

light_levels <- unique(preds2$light)
light_10 <- sample(x = light_levels, size =  10, replace = FALSE)
# tpc_temps <- c(1, 5, 10, 15, 20, 25, 30, 35, 40)

cov_mat <- cov(growth_16)
pca_res <- prcomp(cov_mat, center = TRUE,scale. = TRUE)
summary(pca_res)
ggbiplot(pca_res, labels= light_10)

loadings16 <- scores(pca_res,choices=c(1,2))
summary(eigenvals(pca_res))

pcs16 <- as_data_frame((loadings16)) %>% 
	select(PC1)
pc1_16 <- pcs16 %>% 
	mutate(light = light_10) 
pc1_16 %>% 
	ggplot(aes(x = light, y = PC1)) + geom_point() +
	xlim(0, 250) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Light") + geom_line() +
	ylab("PC1 loadings")
ggsave("figures/light-curves-pc1.pdf", width = 6, height = 4)
