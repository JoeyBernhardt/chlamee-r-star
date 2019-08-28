library(here) ## this is a handy package to help you deal with file paths more easily
library(tidyverse)
library(cowplot)
library(broom)
library(readxl)
library(janitor) ## for cleaning up data.frames
library(plotrix)
library(growthTools) ### note this isn't yet available on CRAN, it's from Colin Kremer
library(rootSolve)


all_growth_n2 %>% 
	mutate(estimate = mu) %>%
	filter(nitrate_concentration < 50) %>% 
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x= nitrate_concentration, y= estimate)) + 
	geom_point() +
	# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
	geom_line(data=all_predsn_2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)") + xlim(0,10)

all_predsn_2 %>% 
	# distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE) %>% 
	filter(ancestor_id == "anc5", treatment == "C", nitrate_concentration.x < 10) %>% View
	mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
	ggplot(aes(x= nitrate_concentration.x, y= growth_rate)) + 
	geom_line() +
	facet_grid(treatment ~ ancestor_id) +
	ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)") + xlim(0,10)
	?get.growth.rate()	
	all_preds2 %>% 
		ggplot(aes(x = days, y = .fitted, group = well_plate, color = best.model)) + geom_line() +
		geom_point(aes(x = days, y = y, shape = exponential)) +
		facet_wrap(~ well_plate) + ylab("Ln(RFU)") +xlab("Days")
	# ggsave(here("figures", "all_growth_tools_plots_well_plate.pdf"), width = 45, height = 45)
	ggsave(here("figures", "all_growth_tools_plots_well_plate_AIC_lag.png"), width = 45, height = 45)
	
	
	
	

# try via spline approach -------------------------------------------------

	
	example_data <- read_csv(here("data-processed", "nitrate-abundances-processed.csv")) %>% 
		filter(population != "COMBO") %>% 
		# filter(grepl("C", well_plate)) %>% 
		filter(!is.na(log(RFU)), !is.na(days))
	
	
	example_split <-  example_data %>% 
		# filter(grepl("C", well_plate)) %>% 
		split(.$well_plate)
	
	
	### define Nathaniel's fitting function
	nderiv<-function(fit, x, eps=1e-5){(predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}
	
	spline.slope <- function(df, n=101, eps=1e-5, span=0.5){
		x <- df$days
		y <- df$RFU
		time_point <- seq(min(x), max(x), length=n)
		growth <- nderiv(loess(log(y) ~ x, degree=1, span=span), time_point)
		output <- top_n(data.frame(growth, time_point), n = 1, wt = growth)
		return(output)
	}
	

	
	
	## fit each well
	growth_rates_example <- example_split %>% 
		map_df(spline.slope, .id = "well_plate")
	
	## join the growth rate results back to initial df
	example_results <- left_join(example_data, growth_rates_example)
	
	## plot it, to visualize where the max growth rate was found along the time series
	
	
	fsample <- example_results %>% 
		filter(well_plate == "F02_1") 
	
	example_results %>% 
		mutate(time_point = as.numeric(time_point)) %>% 
		mutate(max_growth_time = round(time_point, digits = 2)) %>% 
		mutate(time = ifelse(days >= (time_point -0.1) & days <= (time_point + 0.1), "max_point", "time")) %>%
		ggplot(aes(x = days, y = log(RFU), color = time)) + geom_point(size = 0.5) +
		# geom_vline(xintercept = max_growth_time, data = example_results) +
		facet_wrap( ~ well_plate) 
	ggsave(here("figures", "spline_practice.png"), width = 45, height = 45)
	
	
	key <- nitrate %>% 
		select(well_plate, population, nitrate_concentration) %>% 
		distinct(well_plate, .keep_all = TRUE)	
	
	treatments <- read_csv(here("data-processed", "nitrate-treatments-processed.csv"))	
	
	
all_results <- left_join(example_results, treatments)	
	
	all_results %>% 
		distinct(nitrate_concentration, growth, time, .keep_all = TRUE) %>% 
		ggplot(aes(x = nitrate_concentration, y = growth)) + geom_point() +
		facet_grid(ancestor_id ~ treatment)
	
	
	monod_fits_spline <- all_results %>% 
		mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
		rename(estimate = growth) %>% 
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
	
	bs_split <- monod_fits_spline %>% 
		select(population, term, estimate) %>% 
		dplyr::ungroup() %>% 
		spread(key = term, value = estimate) %>%
		split(.$population)
	
	
	all_preds_spline <- bs_split %>% ### here we just use the fitted parameters from the Monod to get the predicted values 
		map_df(prediction_function, .id = "population")
	
	
	all_preds_spline2 <- left_join(all_preds_spline, treatments, by = c("population")) %>% 
		distinct(ancestor_id, treatment, nitrate_concentration.x, .keep_all = TRUE)
	
##	Plot the Monod fits
	
	
	all_growth_n2 <- left_join(all_growth_n, treatments)
	
	all_results %>% 
		mutate(estimate = growth) %>%
		# filter(nitrate_concentration < 50) %>% 
		mutate(nitrate_concentration = as.numeric(nitrate_concentration)) %>% 
		ggplot(aes(x= nitrate_concentration, y= estimate)) + 
		geom_point() +
		# geom_errorbar(aes(ymin=estimate-best.se, ymax=estimate + best.se), width=.2) + 
		geom_line(data=all_preds_spline2, aes(x=nitrate_concentration.x, y=growth_rate, color = treatment), size = 1) +
		facet_grid(treatment ~ ancestor_id) +
		ylab("Exponential growth rate (/day)") + xlab("Nitrate concentration (uM)") 
	ggsave(here("figures", "nitrate_monod_spline.png"), width = 12, height = 8)
	
	
	
	monod_wide <-  monod_fits_spline %>% 
		select(population, term, estimate) %>% 
		spread(key = term, value = estimate)
	
	m <- 0.1 ## set mortality rate, which we use in the rstar_solve
	monod_curve_mortality <- function(nitrate_concentration, umax, ks){
		res <- (umax* (nitrate_concentration / (ks+ nitrate_concentration))) - 0.1
		res
	}	
	
	
	## find the R*
	rstars <- monod_wide %>% 
		# mutate(rstar = uniroot.all(function(x) monod_curve_mortality(x, umax, ks), c(0.0, 50))) %>% ## numerical
		mutate(rstar_solve = ks*m/(umax-m)) ## analytical
	
	treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>% 
		clean_names() %>% 
		mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
		filter(population != "cc1629") 
	
	rstars2 <- left_join(rstars, treatments, by = "population")
	
write_csv(rstars2, "data-processed/rstars-nitrate-spline.csv")	
	
	rstars2 %>%
		group_by(treatment) %>% 
		summarise_each(funs(mean, std.error), rstar_solve) %>% 
		ggplot(aes(x = reorder(treatment, mean), y = mean)) + geom_point() +
		geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error),width = 0.1) +
		ylab("R* (umol N)") + xlab("Selection treatment") + geom_point(aes(x = reorder(treatment, rstar_solve), y = rstar_solve, color = ancestor_id), size = 2, data = rstars2, alpha = 0.5) +
		scale_color_discrete(name = "Ancestor")
	
	ggsave(here("figures", "nitrate-r-star-means-spline.png"), width = 6, height = 4)
	