
library(tidyverse)
library(nls.multstart)
library(cowplot)
### estimate salt tolerances from soup to nuts

library(here)
library(janitor)
library(readxl)
library(minpack.lm)
library(nlstools)

treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")


salt_data <- read.csv("data-raw/chlamee-salt-old/chlamee-salt-rfus-time.csv") %>% 
	mutate(population = str_replace(population, "Anc ", "anc")) %>% 
	mutate(well_plate = paste(well, plate, sep = "_")) %>% 
	group_by(population, treatment, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	rename(days = time) %>% 
	filter(!is.na(RFU)) %>% 
	rename(salt_level = treatment)

salt_split <- salt_data %>% 
	arrange(days) %>% 
	split(.$well_plate) 


s2 <- salt_split %>% 
	map_df(rownames_to_column, var = "time_point") 



unique(salt_data$population)

fitting_window <- function(x) {
	growth_rates <- salt_data %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(population, treatment, well_plate) %>% 
		do(tidy(lm(log(RFU + 1) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

windows <- seq(3,8, by = 1)

multi_fits <- windows %>% 
	map_df(fitting_window, .id = "iteration")

#+ fig.width = 8, fig.height = 6
multi_fits %>% 
	filter(term == "days") %>% 
	ggplot(aes(x = number_of_points, y = estimate, group = well_plate)) + geom_point() + geom_line() +
	facet_wrap( ~ treatment)

exp_fits_top <- multi_fits %>%
	filter(term == "days") %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate)

exp_fits_top %>% 
	mutate(salt_concentration = str_replace(treatment, "S", "")) %>% 
	mutate(salt_concentration = as.numeric(salt_concentration)) %>% 
	ggplot(aes(x = salt_concentration, y = estimate)) + geom_point() +
	facet_wrap( ~ population)


salt_new <- exp_fits_top %>% 
	mutate(salt_concentration = str_replace(treatment, "S", "")) %>% 
	mutate(salt_concentration = as.numeric(salt_concentration)) %>% 
	select(-treatment) %>% 
	left_join(., treatments) %>% 
	rename(mu = estimate) %>% 
	rename(concentration = salt_concentration)


fits3 <- salt_new %>% 
	# filter(treatment %in%  c("A")) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(mu ~ a/(1 + exp(-b * (concentration-c))),
												  data = .x,
												  iter = 500,
												  start_lower = c(a = 1, b = -3, c = 3),
												  start_upper = c(a = 2, b = 1, c = 5),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(a = 0.1, b = -10, c = 1),
												  upper = c(a = 10, b = 20, c = 100),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

params3 <- fits3 %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 
p3 <- params3 %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(salt_tolerance = c)

ssplit3 <- params3 %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	split(.$population)
predictions_salt3 <- ssplit3 %>% 
	map_df(predict_growth, .id = "population")

salt_new %>% 
	ggplot(aes(x = concentration, y = mu)) + geom_point() +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt3) +
	geom_vline(aes(xintercept = salt_tolerance), data = p3, color = "lightblue") + 
	facet_wrap( ~ population, scales = "free") + xlim(0, 10) + ylab("Population growth rate") + xlab("Salt concentration (g/L)")
ggsave("figures/predictions-salt-sigmoid-new.png", width = 15, height = 10)

p3 %>% 
	mutate(treatment = str_replace(treatment, "none", "A")) %>% 
	ggplot(aes(x = treatment, y = salt_tolerance)) + geom_point() +
	geom_point(aes(x = treatment, y = salt_tolerance), data = p1, color = "orange") 




fits3 <- salt_new %>% 
	# filter(treatment %in%  c("A")) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(mu ~ a/(1 + exp(-b * (concentration-c))),
												  data = .x,
												  iter = 500,
												  start_lower = c(a = 1, b = -3, c = 3),
												  start_upper = c(a = 2, b = 1, c = 5),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(a = 0.1, b = -10, c = 1),
												  upper = c(a = 10, b = 20, c = 100),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))


### now try via the direct approach

View(exp_fits_top)

time_to_include <- exp_fits_top %>% 
	select(well_plate, number_of_points) 

salt_abund <- left_join(s2, treatments, by = "population") %>% 
	left_join(., time_to_include) %>% 
	filter(time_point <= number_of_points) %>% 
	group_by(population) %>% 
	mutate(N0_mean = mean(N0)) %>% 
	mutate(salt_concentration = str_replace(salt_level, "S", "")) %>% 
	mutate(salt_concentration = as.numeric(salt_concentration)) 
	


salt_abund_split <- salt_abund %>% 
	# filter(population %in% c(11)) %>% 
	filter(!is.na(N0_mean), !is.na(salt_concentration), !is.na(days), !is.na(RFU), !is.na(population)) %>% 
	split(.$population)



fit1 <- nlsLM(RFU ~ N0_mean * exp((a/(1 + exp(-b * (salt_concentration-c))))*(days)),
	  data= salt_abund,  
	  start= c(a = 1, b = -3, c = 3),
	  lower = c(a = 1, b = -3, c = 1),
	  upper = c(a = 2, b = 1, c = 10),
	  control = nls.control(maxiter=1024, minFactor=1/204800000))
summary(fit1)


View(nsplit[[5]])
i <-5
nsplit <- salt_abund_split
results3<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(RFU ~ N0_mean * exp((a/(1 + exp(-b * (salt_concentration-c))))*(days)),
				  data= nsplit[[i]],  
				  start= c(a = 1, b = -3, c = 3),
				  lower = c(a = 1, b = -10, c = 1),
				  upper = c(a = 2, b = 1, c = 10),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results3 <- bind_rows(results3, hold)
}


salt_tol_summary <- results3 %>%
	group_by(population_index) %>% 
	summarise(salt_tol_lower=quantile(c, probs=0.025),
			  salt_tol_upper=quantile(c, probs=0.975),
			  salt_tol_mean = mean(c)) %>% 
	rename(population = population_index) %>% 
	left_join(., treatments)
	

salt_tol_summary %>% 
	ggplot(aes(x = treatment, y = salt_tol_mean)) + geom_point()


### ok forget the direct approach starting abundances are too uneven, back to the growthTools growth

salt <- read.csv("data-processed/salt_growth_rates.csv") %>% 
	separate(salt_level, into = c("S", "concentration"), sep = 1) %>% 
	mutate(concentration = as.numeric(concentration)) %>%
	mutate(population = as.character(population)) %>% 
	left_join(., treatments, by = "population") %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS")))


### try again
library(broom)
library(nls.multstart)
fits <- salt %>% 
	# filter(treatment %in%  c("A")) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(mu ~ a/(1 + exp(-b * (concentration-c))),
												  data = .x,
												  iter = 500,
												  start_lower = c(a = 1, b = -3, c = 3),
												  start_upper = c(a = 2, b = 1, c = 5),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(a = 0.1, b = -10, c = 1),
												  upper = c(a = 10, b = 20, c = 100),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))


params <- fits %>%
	filter(fit != "NULL") %>% 
	unnest(fit %>% map(tidy)) 
 
library(nlstools)
cis <- fits %>% 
unnest(fit %>% map(~ confint2(.x) %>%
				   	data.frame() %>%
				   	rename(., conf.low = X2.5.., conf.high = X97.5..))) %>% 
	group_by(., population) %>%
	mutate(., term = c('a', 'b', 'c')) %>%
	ungroup()

params3 <- merge(params, cis, by = intersect(names(params), names(cis)))


params3 %>% 
	filter(term == "c") %>% 
	ggplot(aes(x = treatment, y = estimate)) + geom_point() +
	geom_pointrange(aes(x = treatment, y = estimate, ymin = conf.low, ymax = conf.high))

p1 <- params %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(salt_tolerance = c)

write_csv(p1, "data-processed/salt-tolerances-sigmoid.csv")


ssplit <- params %>% 
	select(1:5) %>% 
	spread(key = term, value = estimate) %>% 
	split(.$population)

## now we need to plot the fits
predict_growth <- function(df){
	
	predict_salt <- function(x){
		growth.rate <- df$a[[1]]/(1 + exp(-df$b[[1]] * (x-df$c[[1]])))
		return(growth.rate)
	}
	
	pred <- function(x) {
		y <- predict_salt(x)
	}
	
	x <- seq(0, 60, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(salt_concentration = x, 
			   growth_rate = preds)
	
}


predictions_salt <- ssplit %>% 
	map_df(predict_growth, .id = "population")

View(p1)
salt %>% 
	ggplot(aes(x = concentration, y = mu)) + geom_point() +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt) +
	geom_vline(aes(xintercept = salt_tolerance), data = p1, color = "lightblue") + 
	facet_wrap( ~ population, scales = "free") + xlim(0, 10) + ylab("Population growth rate") + xlab("Salt concentration (g/L)")
ggsave("figures/predictions-salt-sigmoid-new2.png", width = 15, height = 10)



# start_lower = c(a = 1, b = -3, c = 3),
# start_upper = c(a = 2, b = 1, c = 5),
# supp_errors = 'N',
# na.action = na.omit,
# lower = c(a = 0.1, b = -10, c = 1),
# upper = c(a = 10, b = 20, c = 100),


nsplit <- salt %>% 
	split(.$population)
results3<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(mu ~ a/(1 + exp(-b * (concentration-c))),
				  data= nsplit[[i]],  
				  start= c(a = 1, b = -3, c = 3),
				  lower = c(a = 1, b = -20, c = 1),
				  upper = c(a = 2, b = 10, c = 35), ## bound the upper limit at 35
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results3 <- bind_rows(results3, hold)
}
unique(results3$population_index)


nsplit <- salt %>% 
	filter(population %in% c(15, 18, 2)) %>% 
	split(.$population)
results4<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(mu ~ a/(1 + exp(-b * (concentration-c))),
				  data= nsplit[[i]],  
				  start= c(a = 1, b = -3, c = 3),
				  lower = c(a = 1, b = -20, c = 5),
				  upper = c(a = 2, b = 10, c = 35), ## bound the upper limit at 35
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results4 <- bind_rows(results4, hold)
}

View(results4)

results4b <- results4 %>% 
	filter(!c %in% c(5, 35))

View(results4b)

results3b <- results3 %>% 
	filter(!population_index %in% c(15, 18, 2)) ## get rid of the poorly fit curves in the first round

results3c <- bind_rows(results3b, results4b)


write_csv(results3c, "data-processed/salt-tolerances-bootstrapped-custom-high-salt.csv") ## saving with the high salt custom fit with an upper limit of 35
write_csv(results3, "data-processed/salt-tolerances-bootstrapped.csv")
results3 <- read.csv("data-processed/salt-tolerances-bootstrapped.csv") %>% 
	rename(population = population_index)
salt_tol_summary <- results3c %>%
	group_by(population_index) %>% 
	summarise(salt_tol_lower=quantile(c, probs=0.025),
			  salt_tol_upper=quantile(c, probs=0.975),
			  salt_tol_mean = mean(c)) %>% 
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	mutate(treatment = str_replace(treatment, "none", "Ancestors"))

salt_tol_ancestors <- results3c %>%
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	filter(treatment == "none") %>% 
	group_by(treatment, ancestor_id, population) %>% 
	summarise(anc_stol_q2.5=quantile(c, probs=0.025),
			  anc_stol_q97.5=quantile(c, probs=0.975),
			  anc_stol_mean = mean(c)) %>% 
	ungroup() %>% 
	select(-treatment, - population)


salt_tol_CI <- results3c %>%
	filter(c < 12) %>% 
	rename(population = population_index) %>% 
	left_join(., treatments) %>% 
	left_join(., salt_tol_ancestors) %>% 
	mutate(change_salt_tol = c - anc_stol_mean) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	summarise(change_salt_tol_lower=quantile(change_salt_tol, probs=0.025),
			  change_salt_tol_upper =quantile(change_salt_tol, probs=0.975),
			  change_salt_tol_mean = mean(change_salt_tol)) %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) 


cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")
salt_tol_CI %>% 
	ggplot(aes(x = treatment, y = change_salt_tol_mean, color = treatment)) + 
	# geom_point(position=position_dodge(width=0.5)) +
	geom_pointrange(size = .2, aes(x = treatment, y = change_salt_tol_mean, ymin = change_salt_tol_lower, ymax = change_salt_tol_upper), position=position_jitter(width=0.1)) +
	# geom_point(size = 1.5, aes(x = treatment, y = change_salt_tol_mean), position=position_jitter(width=0.1)) +
	geom_hline(yintercept = 0) +
	ylab("Trait change relative to ancestor") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position = "none") + coord_cartesian() 
ggsave("figures/salt-tol-trait-change-95CI.png", width = 6, height = 4)

write_csv(salt_tol_CI, "data-processed/salt-tol-CI-boot.csv")
write_csv(salt_tol_CI, "data-processed/salt-tol-CI-boot-12.csv") ## restricting upper limit of CI to 12

salt_tol_summary %>% 
	ggplot(aes(x = treatment, y = salt_tol_mean)) + geom_point(color = "pink") + 
	geom_pointrange(aes(x = treatment, y = salt_tol_mean, ymin = salt_tol_lower, ymax = salt_tol_upper),
					color = "pink", size = 0.5) +
	ylim(0, 20) + ylab("Salt tolerance (ug/L)") + xlab("Treatment")

salt %>% 
	ggplot(aes(x = concentration, y = mu)) + geom_point() +
	geom_line(aes(x = salt_concentration, y = growth_rate), data = predictions_salt) +
	geom_vline(aes(xintercept = salt_tolerance), data = p1, color = "lightblue") + 
	facet_wrap( ~ population, scales = "free") + xlim(0, 10) + ylab("Population growth rate") + xlab("Salt concentration (g/L)")
ggsave("figures/predictions-salt-sigmoid-new2.png", width = 15, height = 10)



