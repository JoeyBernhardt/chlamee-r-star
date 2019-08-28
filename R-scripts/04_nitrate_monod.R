library(tidyverse)
library(cowplot)
library(broom)
library(readxl)
library(janitor)
library(plotrix)
library(here)

plate_info <- read_csv(here("data-processed", "chlamee-nitrate-rstar-plate-info.csv")) %>% 
	unite(col = "well_plate", well, plate, remove = FALSE) %>% 
	mutate(population = ifelse(population == "cc11629", "COMBO", population))



treatments2 <- left_join(plate_info, treatments) %>% 
	separate(n_level, into = c("N", "nitrate_level"), sep = 1) %>% 
	mutate(nitrate_concentration = NA) %>% 
	mutate(nitrate_concentration = case_when(nitrate_level == "1" ~ "5",
											 nitrate_level == "2" ~ "10",
											 nitrate_level == "3" ~ "20",
											 nitrate_level == "4" ~ "40",
											 nitrate_level == "5" ~ "60",
											 nitrate_level == "6" ~ "80",
											 nitrate_level == "7" ~ "100",
											 nitrate_level == "8" ~ "400",
											 nitrate_level == "9" ~ "600",
											 nitrate_level == "10" ~ "1000",
											 TRUE ~ "no")) %>%
	mutate(population_corrected = case_when(plate == "19" & population == "11" ~ "anc5",
								  plate == "19" & population == "anc5" ~ "11",
								  plate == "24" & population == "5" ~ "30",
								  plate == "24" & population == "30" ~ "5",
								  TRUE ~ population)) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	filter(!is.na(ancestor_id))


treatments3 <- treatments2 %>% 
	select(population, plate, nitrate_level, nitrate_concentration, population_corrected)

treatments4 <- left_join(treatments3, treatments, by = c("population_corrected" = "population"))

treatments5 <- left_join(treatments4, plate_info) %>% 
	distinct(well_plate, .keep_all = TRUE)
treatments6 <- left_join(treatments4, plate_info, by = c("population_corrected" = "population")) %>% 
	distinct(well_plate, .keep_all = TRUE) %>% 
	distinct(population_corrected, ancestor_id, treatment) %>% 
	rename(population = population_corrected)

write_csv(treatments2, "data-processed/nitrate-treatments-processed.csv") ### old version until Feb 20, when I realized there was a mistake here
write_csv(treatments6, "data-processed/nitrate-treatments-processed.csv") ### new as of Feb 20
 
 nitrate <- read_csv("data-processed/nitrate-rstar-rfus-time.csv") %>% 
 	mutate(nitrate_concentration = NA) %>% 
 	mutate(nitrate_concentration = case_when(nitrate_level == "1" ~ "5",
 											 nitrate_level == "2" ~ "10",
 											 nitrate_level == "3" ~ "20",
 											 nitrate_level == "4" ~ "40",
 											 nitrate_level == "5" ~ "60",
 											 nitrate_level == "6" ~ "80",
 											 nitrate_level == "7" ~ "100",
 											 nitrate_level == "8" ~ "400",
 											 nitrate_level == "9" ~ "600",
 											 nitrate_level == "10" ~ "1000",
 								   TRUE ~ "no")) %>%
 	mutate(population_corrected = case_when(plate == "19" & population == "11" ~ "anc5",
 								  plate == "19" & population == "anc5" ~ "11",
 								  plate == "24" & population == "5" ~ "30",
 								  plate == "24" & population == "30" ~ "5",
 								  TRUE ~ population)) %>% 
 	select(-population) %>% 
 	rename(population = population_corrected)
 
names(nitrate)

write_csv(nitrate, "data-processed/nitrate-abundances-processed.csv")

nitrate <- read_csv("data-processed/nitrate-abundances-processed.csv")

nitrate %>% 
	filter(is.na(ancestor_id)) %>% View
 
nitrate %>% 
 	# filter(days < 1.4) %>% 
 	ggplot(aes(x = days, y = RFU, color = factor(nitrate_level), group = well_plate)) + geom_point() +
 	geom_line() + facet_wrap( ~ population, scales = "free_y") + scale_color_viridis_d(name = "Nitrate level")
 
nitrate_exp <- nitrate %>% 
	mutate(exponential = case_when(days < 2 ~ "yes",
								   TRUE ~ "no")) %>% 
	filter(exponential == "yes") %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) 


write_csv(nitrate_exp, "data-processed/nitrate-RFUs-time-exponential-phase.csv")


growth_rates <- nitrate_exp %>%
	filter(population != "COMBO") %>% 
	group_by(nitrate_concentration, population, well_plate) %>%
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


growth2 <- left_join(growth_rates, treatments, by = "population") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))


growth2 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate)) + geom_point() +
	facet_grid(treatment ~ ancestor_id) + geom_hline(yintercept = 0) + ylab("Exponential growth rate (/day)") +
	xlab("Nitrate concentration (uM)") 

write_csv(growth2, "data-processed/nitrate-r-star-growth-rates.csv")
growth2 <- read_csv("data-processed/nitrate-r-star-growth-rates.csv")

monod_fits <- growth2 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))


preds <- growth2 %>%
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>%
	group_by(treatment, ancestor_id) %>% 
	do(augment(nls(estimate ~ umax* (nitrate_concentration/ (ks+ nitrate_concentration)),
				   data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				   control = nls.control(maxiter=500, minFactor=1/204800000))))


growth2 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>%
	ggplot(aes(x= nitrate_concentration, y= estimate)) + geom_point() +
	# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
	geom_line(data=preds, aes(x=nitrate_concentration, y=.fitted), color = "purple", size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")


prediction_function <- function(df) {
	
	
	monodcurve<-function(x){
		growth_rate<- (df$umax[[1]] * (x / (df$ks[[1]] +x)))
		growth_rate}
	
	pred <- function(x) {
		y <- monodcurve(x)
	}
	
	x <- seq(0, 1000, by = 1)
	
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


all_preds_n <- bs_split %>% 
	map_df(prediction_function, .id = "population")


all_predsn_2 <- left_join(all_preds_n, treatments, by = c("population"))
growth2 %>% 
	# all_growth_n %>% 
		# mutate(estimate = mu) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + geom_point() +
	# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
	geom_line(data=all_predsn_2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)")
ggsave("figures/nitrate_monod.pdf", width = 15, height = 10)

library(plotrix)

m2 <- left_join(monod_fits, treatments, by = "population")

m2 %>% 
	group_by(treatment, term) %>% 
	summarise_each(funs(mean, std.error), estimate) %>%
	ggplot(aes(x = treatment, y = mean)) + geom_point() + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) +
	facet_wrap( ~ term, scales = "free")
ggsave("figures/nitrate_monod_params.pdf", width = 8, height = 6)



# growth rates with growthTools -------------------------------------------
library(growthTools)

b2 <- nitrate %>% 
	filter(nitrate_level == 1) %>% 
	filter(well_plate == "B06_1") %>% 
	mutate(ln.fluor = log(RFU)) 

nitrate %>%
	filter(plate == 1) %>%
	ggplot(aes(x = days, y = RFU, group = well_plate, color = population)) + geom_point() + geom_line()


res <- get.growth.rate(b2$days, b2$ln.fluor, plot.best.Q =T, id = "B02")

growth_rates_n <- nitrate %>%
	filter(population != "COMBO") %>% 
	mutate(ln.fluor = log(RFU)) %>% 
	# filter(nitrate_level > 1) %>% 
	group_by(well_plate) %>% 
	do(grs=get.growth.rate(x=.$days, y=.$ln.fluor,id=.$well_plate,plot.best.Q=T,fpath="/Users/joeybernhardt/Documents/Narwani/ChlamEE-R-star/figures/growth_curve_fits/nitrate/"))

growth_sum_n <- growth_rates_n %>%
	summarise(well_plate,mu=grs$best.slope,
			  best.model=grs$best.model,best.se=grs$best.se)

all_growth_n <- left_join(growth_sum_n, treatments5, by = c("well_plate")) %>% 
	filter(!is.na(ancestor_id))
write_csv(all_growth_n, "data-processed/exponential_growth_nitrate.csv")

all_growth_n <- read_csv("data-processed/exponential_growth_nitrate.csv")
all_growth_n %>% 
	filter(!is.na(ancestor_id)) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x = nitrate_concentration, y = mu)) + geom_point() +
	facet_grid(treatment ~ ancestor_id) + geom_hline(yintercept = 0) + ylab("Exponential growth rate (/day)") +
	xlab("Nitrate concentration (uM)") 



# Monod fits with the growthTools estimates -------------------------------

monod_fits <- all_growth_n %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	rename(estimate = mu) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))

preds <- all_growth_n %>%
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>%
	rename(estimate = mu) %>% 
	group_by(treatment, ancestor_id) %>% 
	do(augment(nls(estimate ~ umax* (nitrate_concentration/ (ks+ nitrate_concentration)),
				   data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				   control = nls.control(maxiter=500, minFactor=1/204800000))))


bs_split <- monod_fits %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_n <- bs_split %>% 
	map_df(prediction_function, .id = "population")


all_predsn_2 <- left_join(all_preds_n, treatments5, by = c("population")) %>% 
	filter(!is.na(ancestor_id))


all_predsn_2 %>% 
	# filter(treatment == "none", ancestor_id == "anc5") %>% 
	filter(population == 11) %>% View

treatments2 %>% 
	filter(population == 11) %>% View

names(all_predsn_2)

all_predsn_3 <- all_predsn_2 %>% 
	distinct(population, nitrate_concentration.x, nitrate_level, well_plate, .keep_all = TRUE)


	all_growth_n %>% 
	mutate(estimate = mu) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
		# filter(treatment == "none", ancestor_id == "anc5") %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + geom_point() +
	# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
	geom_line(data= all_predsn_3, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment, group = population), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	# facet_wrap(~population) +
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = 0) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)") + xlim(-15, 1000) + ylim(-1, 2.5)
ggsave("figures/nitrate_monod_growth_tools2.pdf", width = 15, height = 10)
	
	m2 <- left_join(monod_fits, treatments5, by = "population")
	
	m2 %>% 
		filter(term == "ks") %>% View
	
	m2 %>% 
		group_by(treatment, term) %>% 
		summarise_each(funs(mean, std.error), estimate) %>%
		ggplot(aes(x = treatment, y = mean)) + geom_point() + 
		geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) +
		facet_wrap( ~ term, scales = "free")
	ggsave("figures/nitrate_monod_params.pdf", width = 8, height = 6)
	
	

# find R* -----------------------------------------------------------------


	
library(rootSolve)	

	m <- 0.1 ## mortality rate
	

	
monod_wide <-  monod_fits %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)


### define the Monod curve, with a mortality rate of 0.1
monod_curve_mortality <- function(nitrate_concentration, umax, ks){
	res <- (umax* (nitrate_concentration / (ks+ nitrate_concentration))) - 0.1
	res
}	

m <- 0.1 ## set mortality rate, which we use in the rstar_solve

rstars <- monod_wide %>% 
	# mutate(rstar = uniroot.all(function(x) monod_curve_mortality(x, umax, ks), c(0.0, 50))) %>% ## numerical
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

write_csv(monod_wide, "data-processed/nitrate_monod_parameters.csv")

monod_wide <- read_csv("data-processed/nitrate_monod_parameters.csv")
monod_wide <- monod_fits %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)


find_rstar <- function(m) {
	rstar <- monod_wide %>% 
		mutate(rstar_solve = ks*m/(umax-m)) %>% 
		mutate(mortality_rate = m)
	return(rstar)
}

ms <- seq(0.00, 0.3, by = 0.01)

all_rstars <- ms %>% 
	map_df(find_rstar, .id = "ms")


rstars2 <- left_join(all_rstars, treatments, by = "population")
rstars2 %>% 
	filter(!is.na(rstar_solve)) %>% 
	# filter(population != 3) %>% 
	group_by(treatment, mortality_rate) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean, color = mortality_rate)) + geom_point() + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error, color = mortality_rate), width = 0.2) +
	ylab("R* (uM N)") + xlab("Selection treatment") + scale_color_viridis_c() +
	facet_wrap( ~ mortality_rate, scales = "free") +
	geom_point(shape = 1, color = "black")

ggsave("figures/nitrate-r-star-mortality-rates.pdf", width = 20, height = 10)


rstars2 %>% 
	# filter(!is.na(rstar_solve)) %>% 
	# filter(treatment == "B") %>% 
	group_by(treatment, mortality_rate) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = mortality_rate, y = mean, color = mortality_rate)) + geom_point() + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error, color = mortality_rate), width = 0.01) +
	ylab("R* (uM N)") + xlab("Mortality rate") + scale_color_viridis_c() +
	facet_wrap( ~ treatment, scales = "free") +
	geom_point(shape = 1, color = "black") + geom_smooth(method = "lm")
ggsave("figures/nitrate-r-star-mortality-rate-effect.pdf", width = 8, height = 6)




monod_wide <- read_csv("data-processed/nitrate_monod_parameters.csv")

rstars_a <- monod_wide %>% 
	# mutate(rstar = uniroot.all(function(x) monod_curve_mortality(x, umax, ks), c(0.0, 50))) %>% ## numerical
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical


treatments_old <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629") %>% 
	mutate(population_old = population) %>% 
	mutate(treatment_old = treatment) 
	# filter(population != "cc1690") %>% 

treatments_new <- read_csv(here("data-processed", "nitrate-treatments-processed.csv")) %>% 
	distinct(population, ancestor_id, treatment) %>%
	# filter(population != "cc1690") %>% 
	mutate(population_new = population) %>% 
	mutate(treatment_new = treatment)

left_join(treatments_new, treatments_old) %>% 
	select(contains("treatment"), everything()) %>% View

treatments7 <- treatments6 %>% 
	mutate(population_new = population)
left_join(treatments7, treatments_old) %>% 
	select(contains("treatment"), everything()) %>% View


rstars3 <- left_join(rstars_a, treatments_new, by = "population") %>% 
	distinct(population, ks, umax, .keep_all = TRUE)



rstars3 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars3, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor") + geom_point()
ggsave("figures/nitrate-r-star-means.pdf", width = 6, height = 4)
ggsave("figures/nitrate-r-star-means.png", width = 6, height = 4)



# nitrate logistic --------------------------------------------------------
library(nls.multstart)

nitrate <- read_csv("data-processed/nitrate-abundances-processed.csv")

ldata_n0 <- nitrate %>% 
	group_by(well_plate) %>% 
	mutate(N0 = min(RFU)) %>% 
	ungroup()


fits_many_n0 <- ldata_n0 %>% 
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


fits_many_lowN <- ldata_n0 %>% 
	filter(nitrate_concentration < 50) %>% 
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
												  upper = c(K = 300, r = 4),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))


fits_many <- fits_many_n0
fits_many <- fits_many_lowN

fits_well_plate <- fits_many %>% 
	select(well_plate)

key <- nitrate %>% 
	select(well_plate, population, nitrate_concentration) %>% 
	distinct(well_plate, .keep_all = TRUE)	

pops <- left_join(fits_well_plate, key)

info <- fits_many %>%
	unnest(fit %>% map(glance))

# get params
params_low <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 


# params %>% 
# 	filter(term == "r") %>% View

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

library(broom)
?augment

key <- ldata_n0 %>% 
	select(well_plate, population, nitrate_concentration) %>% 
	distinct(well_plate, .keep_all = TRUE)

preds3 <- left_join(preds2, key, by = c("well_plate", "population"))


ggplot() +
	geom_point(aes(days, RFU, color = factor(well_plate)), size = 2, data = filter(ldata_n0, population == 10, nitrate_concentration == 5)) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(well_plate)), data = filter(preds3, population == 10, nitrate_concentration == 5)) +
	facet_wrap(nitrate_concentration ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave("figures/nitrate_rfus_logistic_n0.pdf", width = 40, height = 35)


r_estimates_low <- params_low %>% 
	filter(term == "r")

treatments <- read_csv(here("data-processed", "nitrate-treatments-processed.csv"))


r2_low <- left_join(r_estimates_low, key)
r3_low <- left_join(r2_low, treatments) %>% 
	filter(!is.na(ancestor_id))

write_csv(r3, "data-processed/nitrate_exp_growth_logistic.csv")

r3_low %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate)) + geom_point() +
	facet_grid(ancestor_id ~ treatment) + geom_point(aes(x = nitrate_concentration, y = estimate), data = r3, color = "green", alpha = 0.5, shape = 1) +
	xlim(0, 100) +
	geom_point(aes(x = nitrate_concentration.x, y = estimate), data = growth2, color = "blue", alpha = 0.5, shape = 1)

ggsave("figures/nitrate_logistic_high_low_K_vs_exp.png", width = 10, height = 6)


ggplot() +
geom_point(aes(x = nitrate_concentration.x, y = estimate), data = growth2, color = "blue", alpha = 0.5, shape = 1)+
	facet_grid(ancestor_id ~ treatment)

growth2 %>% 
	group_by(nitrate_concentration.x, ancestor_id, treatment) %>% 
	distinct(nitrate_concentration.x, population, estimate) %>% 
	tally()

write_csv(growth2, "data-processed/nitrate_exp_growth_eyeball.csv")

monod_fits_logistic <- r3 %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	group_by(population) %>% 
	do(tidy(nls(estimate ~ umax* (nitrate_concentration / (ks+ nitrate_concentration)),
				data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
				control = nls.control(maxiter=500, minFactor=1/204800000))))


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

bs_split <- monod_fits_logistic %>% 
	select(population, term, estimate) %>% 
	dplyr::ungroup() %>% 
	spread(key = term, value = estimate) %>%
	split(.$population)


all_preds_logistic <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
	map_df(prediction_function, .id = "population")


all_preds_logistic2 <- left_join(all_preds_logistic, treatments) %>% 
	distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)

##	Plot the Monod fits


# all_growth_n2 <- left_join(r3, treatments)

r3 %>% 
	# filter(nitrate_concentration < 50) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
	geom_line(data=all_preds_logistic2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)") 
ggsave(here("figures", "nitrate_monod_logistic.png"), width = 12, height = 8)

r3 %>% 
	filter(population == 10, nitrate_concentration == 5) %>% View

monod_wide <-  monod_fits_logistic %>% 
	select(population, term, estimate) %>% 
	spread(key = term, value = estimate)

m <- 0.1 ## set mortality rate, which we use in the rstar_solve


## find the R*
rstars <- monod_wide %>% 
	# mutate(rstar = uniroot.all(function(x) monod_curve_mortality(x, umax, ks), c(0.0, 50))) %>% ## numerical
	mutate(rstar_solve = ks*m/(umax-m)) ## analytical

pop_key <- treatments %>% 
	select(population, nitrate_concentration, ancestor_id, treatment) %>% 
	distinct(population, .keep_all = TRUE)

rstars2 <- left_join(rstars, pop_key, by = "population")

write_csv(rstars2, "data-processed/rstars-nitrate-spline.csv")	

rstars2 %>%
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar_solve) %>% 
	ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
	ylab("R* (umol N)") + xlab("Selection treatment") +
	geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2, alpha = 0.5) +
	scale_color_discrete(name = "Ancestor")

ggsave(here("figures", "nitrate-r-star-means-logistic.png"), width = 6, height = 4)



# compare growth rate estimates -------------------------------------------

logistic <- read_csv("data-processed/nitrate_exp_growth_logistic.csv") %>% 
	select(population, ancestor_id, treatment, nitrate_concentration, estimate) %>% 
	mutate(method  = "logistic") %>% 
	distinct(population, ancestor_id, treatment, nitrate_concentration, estimate) %>% 
	filter(!is.na(ancestor_id)) %>% 
	mutate(population = as.character(population)) %>% 
	mutate(method = "intrinsic rate of incr via logistic")
eyeball <- read_csv("data-processed/nitrate_exp_growth_eyeball.csv") %>% 
	select(population, ancestor_id, treatment, nitrate_concentration.x, estimate) %>%
	rename(nitrate_concentration = nitrate_concentration.x) %>% 
	distinct(population, ancestor_id, treatment, nitrate_concentration, estimate) %>% 
	filter(!is.na(population))%>% 
	mutate(population = as.character(population)) %>% 
	mutate(method  = "eyeball-exponential")
AIC_exp <- read_csv("data-processed/nitrate_exp_growth_w_growthtools_AIC.csv") %>% 
	rename(estimate = mu) %>% 
	mutate(method = "exponential via pw regression and AIC") %>% 
	mutate(population = as.character(population))
AIC_exp <- left_join(AIC_exp, treatments)

AICc_exp <- read_csv("data-processed/nitrate_exp_growth_w_growthtools.csv") %>% 
	rename(estimate = mu) %>% 
	mutate(method = "exponential via pw regression and AICc") %>% 
	mutate(population = as.character(population))
AICc_exp <- left_join(AICc_exp, treatments)

all_growth_estimates <- bind_rows(logistic, eyeball, AIC_exp, AICc_exp)


all_growth_estimates %>%
	group_by(population, nitrate_concentration) %>% 
	ggplot(aes(x = nitrate_concentration, y = estimate, color = method)) + geom_point(shape = 1) +
	facet_grid(ancestor_id ~ treatment)
ggsave("figures/comparison_nitrate_monod.png", width = 15, height = 6)


all_growth_estimates %>% 
	group_by(population, nitrate_concentration, ancestor_id, treatment, method) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(aes(x = nitrate_concentration, y = mean, color = method)) + geom_point(shape=1) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	facet_grid(ancestor_id ~ treatment) + ylab("Mean growth rate (/day)") + xlab("Nitrate (uM)")
ggsave("figures/comparison_nitrate_monod_means.png", width = 17, height = 6)


all_growth_estimates %>% View
