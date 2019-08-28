

library(tidyverse)
library(purrr)
library(nls.multstart)
library(here)
library(broom)
library(cowplot)

## fit logistic with no N0


nitrate <- read_csv(here("data-processed", "nitrate-abundances-processed.csv")) %>% 
	group_by(population, nitrate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	filter(population != "COMBO") 
key <- nitrate %>% 
	select(well_plate, population, nitrate_concentration, RFU) %>% 
	group_by(well_plate, population, nitrate_concentration) %>% 
	summarise(N0 = RFU[[1]]) %>% 
	distinct(well_plate, .keep_all = TRUE)	

fits_many <- nitrate %>% 
	group_by(population, well_plate) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ K/(1 + (K/N0- 1)*exp(-r*days)),
												  data = .x,
												  iter = 500,
												  start_lower = c(K = 100, N0 = 5, r = 0),
												  start_upper = c(K = 500, N0 = 10, r = 2),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(K = 10, N0 = 5, r = 0),
												  upper = c(K = 4000, N0 = 30, r = 6),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

info <- fits_many %>%
	unnest(fit %>% map(glance))



nc <- info %>% 
	filter(isConv == "FALSE") %>% 
	select(well_plate)

params <- fits_many %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 
 

params_spread <- params %>% 
	select(well_plate, population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	left_join(., key)

params_spread %>% 
	ggplot(aes(x = nitrate_concentration, y = r)) + geom_point()

params_split <- params_spread %>% 
	split(.$well_plate)

new_preds <- nitrate %>%
	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 100), stringsAsFactors = FALSE))


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
	ungroup() %>% 
	left_join(., key)



ggplot() +
	geom_point(aes(days, RFU, color = factor(nitrate_concentration)), size = 2, data = nitrate) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(nitrate_concentration)), data = preds2b) +
	facet_grid(nitrate_concentration ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave("figures/logisitic_fits_fitted_n0.png", width = 30, height = 20)

logistic_function <- function(r, N, K) {
	RFU <-  K/(1 + (K/N0 - 1)*exp(-r*days))
}


derivative <- function(f, x, ..., order = i, delta = 0.1, sig = 6) {
	# Numerically computes the specified order derivative of f at x
	vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
	grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
	vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
	for (i in 2:(order + 1)) {
		for (j in 1:(order - i + 2)) {
			stepsize <- grid[i + j - 1] - grid[i + j - 2]
			vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
		}
	}
	return(signif(vals[order + 1, 1], sig))
}



get_derivatives <- function(data) {
	all <- data
	x <- seq(0, 4, by = 0.01)
	tpc<-function(x){
		res <- all$K[1]/(1 + (all$K[1]/all$N0[1] - 1)*exp(-all$r[1]*x))
		res
	}
	
	der <- sapply(x, derivative, f = tpc, order = 1)
	derivatives <- data.frame(x, der) 
}

all <- params_spread %>% 
	filter(well_plate == "C09_13")

(N <- all$K[1]/(1 + (all$K[1]/all$N0[1] - 1)*exp(-all$r[1]*2.71)))

dndt <- (all$r[1]*N*(1-(N/all$K[1])))/N


all_derivatives <- params_split %>% 
	map_df(get_derivatives, .id = "well_plate")

max_der <- all_derivatives %>% 
	group_by(well_plate) %>% 
	filter(der == max(der)) 

all_params_der <- left_join(max_der, params_spread) %>% 
	mutate(N = K/(1 + (K/N0- 1)*exp(-r*x))) %>% 
	mutate(dndt_N = r*(1-(N/K))) %>% 
	mutate(max_der_percapita = der/N) %>% 
	left_join(., key)


all_params_der %>% 
	ggplot(aes(x = nitrate_concentration, y = max_der_percapita)) + geom_point() +
	facet_wrap( ~ population)

length(unique(max_der$well_plate))

nitrate_snip <- nitrate %>% 
	filter(well_plate == "D03_21") %>% 
	rename(x = days)

preds_snip <- preds2b %>%
	filter(well_plate == "D03_21") 
	

all_derivatives %>% 
	filter(well_plate == "D03_21") %>% 
	ggplot(aes(x = x,  y = der)) + geom_point() + 
	geom_point(aes(x = x, y = RFU), data = nitrate_snip, color = 'purple') +
	geom_line(aes(x = days, y = RFU), data = preds_snip, color = 'purple')



# Fixed N0 ----------------------------------------------------------------

nitrate <- read_csv(here("data-processed", "nitrate-abundances-processed.csv")) %>% 
	group_by(population, nitrate_concentration, well_plate) %>%
	mutate(N0 = RFU[[1]]) %>%
	filter(population != "COMBO") 

fits_many_N0 <- nitrate %>% 
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
												  upper = c(K = 4000, r = 6),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

info_N0 <- fits_many_N0 %>%
	unnest(fit %>% map(glance))


params_n0 <- fits_many_N0 %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 

params_spread_n0 <- params_n0 %>% 
	select(well_plate, population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	left_join(., key)

write_csv(params_spread_n0, "data-processed/rk-logistic-params.csv")

params_spread_n0 <- read_csv("data-processed/rk-logistic-params.csv") %>% 
	filter(K != "4000")

predictN <- function(all) {
	x <- seq(0, 4, by = 0.01)
	tpc<-function(x){
		res <- all$K[1]/(1 + (all$K[1]/all$N0[1] - 1)*exp(-all$r[1]*x))
		res
	}
	
	N <- sapply(x, tpc)
	df <- data.frame(x, N)
	return(df)
}

splits <- params_spread_n0 %>% 
	split(.$well_plate)

predictions <- splits %>% 
	map_df(predictN, .id = "well_plate") 

predictions2 <- predictions %>% 
	left_join(., key) %>% 
	rename(RFU = N,
		   days = x)




# max_min <- group_by(nitrate, well_plate) %>%
# 	summarise(., min_days = min(days), max_days = max(days)) %>%
# 	ungroup()
# 
# new_preds <- nitrate %>%
# 	do(., data.frame(days = seq(min(.$days), max(.$days), length.out = 10), stringsAsFactors = FALSE))


# preds2b_n0 <- fits_many_N0 %>%
# 	unnest(fit %>% map(augment, newdata = new_preds)) %>%
# 	merge(., max_min, by = 'well_plate') %>%
# 	group_by(., well_plate) %>%
# 	filter(., days > unique(min_days) & days < unique(max_days)) %>%
# 	rename(., RFU = .fitted) %>%
# 	ungroup() %>% 
# 	left_join(., key)

ggplot() +
	geom_point(aes(days, RFU, color = factor(nitrate_concentration)), size = 2, data = nitrate) +
	geom_line(aes(days, RFU, group = well_plate, color = factor(nitrate_concentration)), data = predictions2) +
	facet_grid(nitrate_concentration ~ population, scales = "free") +
	scale_color_viridis_d() +
	ylab('RFU') +
	xlab('Days')
ggsave("figures/logisitic_fits_fitted_n0_fixed.png", width = 30, height = 20)



params_split_n0 <- params_spread_n0 %>% 
	split(.$well_plate)



get_derivatives_fixed <- function(data) {
	all <- data
	x <- seq(0, 4, by = 0.01)
	tpc<-function(x){
		res <- all$K[1]/(1 + (all$K[1]/all$N0[1] - 1)*exp(-all$r[1]*x))
		res
	}
	
	der <- sapply(x, derivative, f = tpc, order = 1)
	derivatives <- data.frame(x, der) 
}

all_derivatives <- params_split_n0 %>% 
	map_df(get_derivatives_fixed, .id = "well_plate")

all_derivatives %>% 
	left_join(., key) %>% 
	left_join(., params_spread_n0) %>% 
	mutate(N = K/(1 + (K/N0- 1)*exp(-r*x))) %>% 
	mutate(dndt_N = r*(1-(N/K))) %>%
	mutate(der_N = der/N) %>%
	ggplot(aes(x = x, y = der_N, color = nitrate_concentration)) + geom_point() +
	facet_wrap( ~ population)
ggsave("figures/derivative_time_nitrate.png", width = 20, height = 30)

max_der <- all_derivatives %>% 
	left_join(., key) %>% 
	left_join(., params_spread_n0) %>% 
	mutate(N = K/(1 + (K/N0- 1)*exp(-r*x))) %>% 
	mutate(dndt_N = r*(1-(N/K))) %>%
	mutate(der_N = der/N) %>% 
	group_by(well_plate, nitrate_concentration, population, r) %>% 
	summarise_each(funs(max), der_N, der)

all_params_der <- left_join(max_der, params_spread_n0) %>% 
	mutate(N = K/(1 + (K/N0- 1)*exp(-r*x))) %>% 
	mutate(dndt_N = r*(1-(N/K))) %>% 
	mutate(max_der_percapita = der/N) %>% 
	left_join(., key)

all_params_der %>% 
	ggplot(aes(x = r, y = dndt_N, color = nitrate_concentration)) + geom_point()
ggsave("figures/max_deriv_r.png", width = 6, height = 4)

max_der %>% 
	ggplot(aes(x = nitrate_concentration, y = der_N)) + geom_point() +
	facet_wrap( ~ population)
ggsave("figures/max_per_capita_deriv_r.png", width = 26, height = 24)
	


derivs <- all_derivatives %>% 
	left_join(., key) %>% 
	left_join(., params_spread_n0) %>% 
	mutate(N = K/(1 + (K/N0- 1)*exp(-r*x))) %>% 
	mutate(dndt_N = r*(1-(N/K))) %>%
	mutate(der_N = der/N)


max_der %>% 
	ggplot(aes(x = r, y = der_N, color = nitrate_concentration)) + geom_point() +
	geom_abline(slope = 1, intercept = 0)

derivs %>% 
	# filter(population == 30) %>% 
	select(population, nitrate_concentration, x, der, der_N, well_plate, N) %>% 
	gather(key = type, value = derivative, der, der_N, N) %>% 
	# filter(grepl("B", well_plate)) %>% 
	ggplot(aes(x = x, y = derivative, group = well_plate, color = factor(nitrate_concentration))) +
	geom_line() +
	facet_grid(type ~ population, scales = "free") + scale_color_viridis_d()
ggsave("figures/derivs-nitrate.png", width = 46, height = 26)


all_params_der %>% 
	filter(K != "4000") %>% 
	ggplot(aes(x = nitrate_concentration, y = max_der_percapita)) + geom_point() +
	facet_wrap( ~ population) + ylab("Maximum derivative / N") + xlab("Nitrate concentration (uM)")
ggsave("figures/max_derivative_nitrate_concentration.png", width = 30, height = 20)
