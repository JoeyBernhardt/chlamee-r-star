library(readxl)
library(tidyverse)
library(minpack.lm)
library(here)
library(janitor)
library(cowplot)
library(broom)
library(tidyverse)

treatments <- read_excel(here("data-general", "ChlamEE_Treatments_JB.xlsx")) %>%
	clean_names() %>%
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>%
	filter(population != "cc1629")

phosphate_abundances <- read_csv("data-processed/phosphate-rstar-rfus-time.csv") %>% 
	filter(population != "COMBO") %>% 
	group_by(population, phosphate_concentration, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	filter(population != "COMBO") %>% 
	filter(!is.na(RFU)) %>% 
	group_by(population) %>% 
	mutate(N0_mean = mean(N0)) 

phosphate_split <- phosphate_abundances %>% 
	arrange(days) %>% 
	split(.$well_plate) 


p2 <- phosphate_split %>% 
	map_df(rownames_to_column, var = "time_point") 


# fit sequentially to identify the exponential phase ----------------------

phosphate_abundances %>% 
	ggplot(aes(x = days, y = log(RFU), color = phosphate_concentration, group = well_plate)) + geom_point() +
	# geom_smooth(method = "lm") +
	geom_line() +
	facet_wrap( ~ population)


fitting_window <- function(x) {
	growth_rates <- phosphate_abundances %>% 
		top_n(n = -x, wt = days) %>% 
		group_by(phosphate_concentration, population, well_plate) %>% 
		do(tidy(lm(log(RFU) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

fitting_window_log_linear <- function(x) {
	growth_rates <- phosphate_abundances %>% 
		group_by(phosphate_concentration, population, well_plate) %>% 
		top_n(n = -x, wt = days) %>% 
		do(tidy(lm(log(RFU) ~ days, data = .))) %>% 
		mutate(number_of_points = x) %>% 
		ungroup()
}

windows <- seq(3,7, by = 1)

multi_fits <- windows %>% 
	map_df(fitting_window_log_linear, .id = "iteration")

#+ fig.width = 8, fig.height = 6
multi_fits %>% 
	filter(term == "days") %>% 
	ggplot(aes(x = number_of_points, y = estimate, group = well_plate)) + geom_point() + geom_line() +
	facet_wrap( ~ phosphate_concentration)

exp_fits_top <- multi_fits %>%
	filter(term == "days") %>% 
	group_by(well_plate) %>% 
	top_n(n = 1, wt = estimate)

write_csv(exp_fits_top, "data-processed/phosphate-sequential-fits.csv")


exponential_fits <- read_csv("data-processed/phosphate-sequential-fits.csv")

phosphate <- left_join(phosphate_abundances, treatments, by = "population") %>% 
	left_join(., p2) %>% 
	left_join(., exponential_fits) %>% 
	filter(time_point <= number_of_points)
write_csv(phosphate, "data-processed/phosphate-exponential.csv")	

phosphate %>% 
	# filter(population == 16, phosphate_concentration == 50) %>% View
	ggplot(aes(x = days, y = RFU, group = well_plate, color = factor(phosphate_concentration))) + geom_point() +
	scale_color_viridis_d() + facet_wrap( ~ population) + geom_line()
ggsave("figures/phosphate-exponential.png", width = 20, height = 15)




phosphate_n0 <- phosphate %>% 
	group_by(population) %>% 
	mutate(N0_mean = mean(N0))

# phosphate_n0 %>% 
# 	filter(population == 10 & phosphate_concentration == 50) %>% 
# 	ggplot(aes(x = days, y = RFU, color = well_plate)) + geom_point()
# 	mutate(outlier = ifelse(population == 10 & phosphate_concentration == 50, "yes", "no")) %>% View

results_CI <- phosphate_n0 %>% 
	filter(well_plate != "D11_30") %>% ## remove weirdo outlier well
	group_by(ancestor_id, treatment, population) %>% 
	do(tidy(nlsLM(RFU ~ N0_mean * exp((umax*(phosphate_concentration/ (ks+ phosphate_concentration)))*(days)),
				  data= .,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000)), conf.int = TRUE))

write_csv(results_CI, "data-processed/phosphate-monod-parameters.csv")
results_CI <- read_csv("data-processed/phosphate-monod-parameters.csv")

results2 <- phosphate_n0 %>% 
	group_by(ancestor_id, treatment, population) %>% 
	do(fit = nlsLM(RFU ~ N0_mean * exp((umax*(phosphate_concentration/ (ks+ phosphate_concentration)))*(days)),
				  data= .,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))) %>% 
	as.data.frame(nlsBoot(.$fit))

nlsBoot(results2$fit[[1]])

m <- 0.1
res_wide <- results_CI %>% 
	select(ancestor_id, treatment, population, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	mutate(rstar = ks*m/(umax-m))

write_csv(res_wide, "data-processed/phosphate-rstars-direct.csv")

results_CI %>% 
	filter(term == "ks") %>% 
	ggplot(aes(x = treatment, y = estimate, color = ancestor_id)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, alpha = 0.5)


library(plotrix)
res_wide %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar, umax) %>% 
	ggplot(aes(x = reorder(treatment, rstar_mean), y = rstar_mean)) + geom_point() +
	geom_errorbar(aes(ymin = rstar_mean - rstar_std.error, ymax = rstar_mean + rstar_std.error), width = 0.1) +
	ylab("R* (uM P)") + xlab("Selection treatment")
ggsave("figures/phosphate-rstar-direct.png", width = 8, height = 6)


res_wide %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), rstar, umax) %>% 
	ggplot(aes(x = reorder(treatment, umax_mean), y = umax_mean)) + geom_point() +
	geom_errorbar(aes(ymin = umax_mean - umax_std.error, ymax = umax_mean + umax_std.error), width = 0.1) +
	ylab("umax (per day)") + xlab("Selection treatment")
ggsave("figures/phosphate-umax-direct.png", width = 8, height = 6)


mod1 <- lm(umax ~ treatment, data = res_wide)
summary(mod1)


res_wide %>% 
	group_by(treatment) %>% 
	ggplot(aes(x = ancestor_id, y = rstar, color = treatment)) + geom_point() +
	# geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) +
	ylab("R* (uM P)") + xlab("Selection treatment")



# plot the fits -----------------------------------------------------------

predict_abundance <- function(df){

	growth_fun <- function(x){
	abundance <-  df$N0_mean[[1]] * exp((df$umax[[1]]*(df$phosphate_concentration[[1]]/ (df$ks[[1]]+ df$phosphate_concentration[[1]])))*(x))
	abundance
}

pred <- function(x) {
	y <- growth_fun(x)
}

x <- seq(0, 3, by = 0.1)

preds <- sapply(x, pred)
preds <- data.frame(x, preds) %>% 
	rename(days.x = x, 
		   abundance = preds)
}


phosphate_setup <- phosphate_n0 %>% 
	group_by(ancestor_id, treatment, population, well_plate) %>% 
	distinct(phosphate_concentration, N0_mean) 

phosphate_setup2 <- left_join(phosphate_setup, res_wide)

phosphate_split <- phosphate_setup2 %>% 
	split(.$well_plate)

abundances <- phosphate_split %>% 
	map_df(predict_abundance, .id = "well_plate")

abundances2 <- left_join(abundances, phosphate_setup)

abundances2 %>% 
	filter(ancestor_id == "anc2", treatment == "L") %>% 
	ggplot(aes(x = days.x, y = abundance, group = well_plate, color = factor(phosphate_concentration))) + geom_line() + 
	geom_point(aes(x = days, y = RFU, color = factor(phosphate_concentration)), data = phosphate_n0) +
	facet_grid(ancestor_id ~ treatment, scales = "free_y") + scale_color_viridis_d() +
	xlab("Days") + ylab("RFU")
ggsave("figures/phosphate-direct-fits-population.png", width = 15, height = 15)


abundances2 %>% 
	filter(ancestor_id == "anc2", treatment == "L") %>%  ## maybe take this one out, major outlier
	ggplot(aes(x = days.x, y = abundance, group = well_plate, color = factor(phosphate_concentration))) + geom_line() + 
	geom_point(aes(x = days, y = RFU, color = factor(phosphate_concentration)),
			   data = filter(phosphate_n0, ancestor_id == "anc2", treatment == "L")) +
	facet_grid(ancestor_id ~ treatment, scales = "free_y") + scale_color_viridis_d() +
	xlab("Days") + ylab("RFU")


abundances2 %>% 
	ggplot(aes(x = days.x, y = abundance, group = well_plate, color = factor(phosphate_concentration))) + geom_line() + 
	geom_point(aes(x = days, y = RFU, color = factor(phosphate_concentration)), data = phosphate_n0) +
	facet_wrap( ~ well_plate, scales = "free_y") + scale_color_viridis_d() +
	xlab("Days") + ylab("RFU")
ggsave("figures/phosphate-direct-fits.png", width = 45, height = 35)


#### bootstrap



fitc <- nlsLM(RFU ~ N0_mean * exp((umax*(phosphate_concentration/ (ks+ phosphate_concentration)))*(days)),
				  data= nsplit[[30]],  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
nls_boot_c <- nlsBoot(fitc, niter = 999)
boots30 <- data.frame(nls_boot_c$coefboot)

summary(fitc)
best_fit_c <- coef(fitc)
nlsResiduals(fitc)

nls_boot_coefs_c <- as_data_frame(nls_boot_c$coefboot)
nls_boot_coefs_c %>% 
	summarise_each(funs(mean, std.error), umax, ks) %>% View
nls_boot_coefs_c %>% 
	summarise(ks_q2.5=quantile(ks, probs=0.025),
		  ks_q97.5=quantile(ks, probs=0.975),
		  ks_mean = mean(ks)) %>% View


df <- nsplit[[12]]

get_bootstrap <- function(df){
	nls_boot_c <- nlsBoot(nlsLM(RFU ~ N0_mean * exp((umax*(phosphate_concentration/ (ks+ phosphate_concentration)))*(days)),
				  data= df,  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000)), niter = 999)
	# nls_boot_coefs_c <- as_tibble(nls_boot_c$coefboot)
	# nls_boot_coefs_c <- nls_boot_coefs_c %>% 
	# 	mutate(population = unique(df$population))
	
}


df <- nsplit[[1]]

library(nlstools)
boot_fits <- nsplit %>% 
	map(get_bootstrap)

boot_fits

nsplit <- phosphate_n0 %>% 
	filter(!is.na(N0_mean), !is.na(phosphate_concentration), !is.na(days), !is.na(RFU), !is.na(population)) %>% 
	split(.$population)

thing <- phosphate_n0 %>% 
	filter(!is.na(N0_mean), !is.na(phosphate_concentration), !is.na(days), !is.na(RFU), !is.na(population)) 


View(nsplit[[5]])
library(nlstools)
i <-5
results3<-data.frame()
for(i in 1:length(nsplit)){
	hold <- nlsLM(RFU ~ N0_mean * exp((umax*(phosphate_concentration/ (ks+ phosphate_concentration)))*(days)),
				  data= nsplit[[i]],  
				  start= c(umax = 1, ks = 1),
				  lower = c(umax = 0, ks= 0),
				  upper = c(umax = 3, ks = 100),
				  control = nls.control(maxiter=1024, minFactor=1/204800000))
	hold <- nlsBoot(hold, niter = 999)
	hold <- data.frame(hold$coefboot)
	hold$replicate <- rownames(hold)
	hold$population_index <- unique(nsplit[[i]]$population)
	results3 <- bind_rows(results3, hold)
	}

unique(results3$population_index)

results4 <- results3 %>% 
	rename(population = population_index)
write_csv(results4, "data-processed/phosphate-params-monod-bootstrapped.csv")

### now let's get the 95% on the r*

results4 <- read.csv("data-processed/phosphate-params-monod-bootstrapped.csv") %>% 
	mutate(population = as.character(population))
m <- 0.56

	

unique(results4$population)
m <- 0.56
ancestors_p <- results4 %>%
	left_join(., treatments) %>% 
	filter(treatment == "none") %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	summarise(anc_rstar_q2.5=quantile(rstar, probs=0.025),
			  anc_rstar_q97.5=quantile(rstar, probs=0.975),
			  anc_rstar_mean = mean(rstar)) %>% 
	ungroup() %>% 
	select(-treatment, - population)

rstar_CI_p <- results4 %>%
	left_join(., treatments) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	left_join(., ancestors_p) %>% 
	mutate(change_pstar = rstar - anc_rstar_mean) %>% 
	group_by(treatment, ancestor_id, population) %>% 
	summarise(change_rstar_lower=quantile(change_pstar, probs=0.025),
			  change_rstar_upper=quantile(change_pstar, probs=0.975),
			  change_rstar_mean = mean(change_pstar)) %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(treatment == "none", "Ancestors", treatment)) %>% 
	mutate(treatment = factor(treatment,
							  levels=c("Ancestors", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) 
cols_no_anc  <- c("#6b6b6b", "#f9a729", "#97cfd0", "#00a2b3", "#f1788d", "#cf3e53", "#b9ca5d")
rstar_CI_p %>% 
	ggplot(aes(x = treatment, y = change_rstar_mean, color = treatment)) + 
	# geom_point(position=position_dodge(width=0.5)) +
	geom_pointrange(size = .2, aes(x = treatment, y = change_rstar_mean, ymin = change_rstar_lower, ymax = change_rstar_upper), position=position_jitter(width=0.1)) +
	geom_hline(yintercept = 0) +
	ylab("Trait change relative to ancestor") +
	xlab("") + scale_color_manual(values = c("black", cols_no_anc)) +
	theme(legend.position = "none")
ggsave("figures/p-star-trait-change-95CI.png", width = 6, height = 4)

	write_csv(rstar_CI_p, "data-processed/change-pstar-monod-boot-0.56.csv")


predict_monod_p <- function(df) {
	
	monodcurve <- function(x){
		growth_rate<- (df$umax[[1]] * (x / (df$ks[[1]] +x)))
		growth_rate}
	
	pred <- function(x) {
		y <- monodcurve(x)
	}
	
	x <- seq(0, 50, by = 0.01)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(phosphate_concentration.x = x, 
			   growth_rate = preds)
}


monod_split <- res_wide %>% 
	split(.$population)

monod_split <- results3 %>% 
	mutate(unique_id = paste(population_index, replicate, sep = "_")) %>% 
	split(.$unique_id)
write_csv(results3, "data-processed/phosphate-monod-params-bootstrap.csv")
# results3 <- read_csv("data-processed/phosphate-monod-params-bootstrap.csv")


### ok, let's get the 95% CI on the R* estimates

rstar_CI <- results3 %>%
	mutate(unique_id = paste(population_index, replicate, sep = "_")) %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	group_by(population_index) %>% 
	summarise(rstar_q2.5=quantile(rstar, probs=0.025),
			  rstar_q97.5=quantile(rstar, probs=0.975),
			  rstar_mean = mean(rstar)) 

predictions_phosphate <- monod_split %>% 
	map_df(predict_monod_p, .id = "unique_id") 

## make this for the original fits only, not the bootstraps
phosphate_monod_params <- read_csv("data-processed/phosphate-rstars-direct.csv") %>% 
	split(.$population)

predictions_phosphate <- phosphate_monod_params %>% 
	map_df(predict_monod_p, .id = "population") 

predictions_phosphate_treatments <- predictions_phosphate %>% 
	left_join(., treatments)
write_csv(predictions_phosphate_treatments, "data-processed/phosphate-monod-curves.csv")


predictions_phosphate <- predictions_phosphate %>% 
	separate(unique_id, into = c("population", "replicate"))

predictions_phosphate <- predictions_phosphate %>%  
	left_join(., treatments)

predictions_phosphate %>% 
	ggplot(aes(x = phosphate_concentration.x, y = growth_rate)) + geom_line(alpha = 0.1) +
	facet_grid(ancestor_id ~ treatment) + ylab("Growth rate (per day)") + xlab("phosphate (uM N)")
ggsave("figures/phosphate-monod-bootstrap.png", width = 12, height = 8)

names(predictions_phosphate)
predictions_summary <- predictions_phosphate %>% 
	group_by(ancestor_id, treatment, population, phosphate_concentration.x) %>% 
	summarise(growth_rate_q2.5=quantile(growth_rate, probs=0.025),
			  growth_rate_q97.5=quantile(growth_rate, probs=0.975),
			  growth_rate_mean = mean(growth_rate)) 

write_csv(predictions_summary, "data-processed/phosphate-bootstrap-monod-summary.csv")

predictions_summary_p <- read_csv("data-processed/phosphate-bootstrap-monod-summary.csv")
predictions_summary_p %>% 
	ggplot(aes(x = phosphate_concentration.x, y = growth_rate_mean, group = population)) + geom_line() +
	geom_ribbon(aes(ymin = growth_rate_q2.5, ymax = growth_rate_q97.5), fill = "grey70", alpha = 0.5) +
	facet_grid(ancestor_id ~ treatment) +ylab("Growth rate (per day)") + xlab("phosphate (uM P)")
ggsave("figures/phosphate-monod-bootstrap-CI.png", width = 13, height = 7)



predictions_phosphate_treatments %>% 
	filter(treatment == "none") %>% 
	ggplot(aes(x = phosphate_concentration.x, y = growth_rate, group = ancestor_id, color = ancestor_id)) + geom_line(size = 1) +
	# facet_grid(ancestor_id ~ treatment) +
	scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "Ancestor") +
	ylab("Growth rate (per day)") + xlab("Phosphate (uM P)") + xlim(0, 0.7) + ylim(0, 0.6)
ggsave("figures/ancestors-phosphate-monod-low-P.png", width = 6, height = 4)


predictions_summary_p %>% 
	mutate(treatment = ifelse(treatment == "none", "A", treatment)) %>%
	mutate(treatment = factor(treatment,
							  levels=c("A", "C", "L", "N",
							  		 "P", "B", "S", "BS"))) %>%
	# filter(treatment == "none") %>% 
	ggplot(aes(x = phosphate_concentration.x, y = growth_rate_mean, group = ancestor_id, color = treatment)) + geom_line() +
	geom_ribbon(aes(ymin = growth_rate_q2.5, ymax = growth_rate_q97.5, fill = treatment), alpha = 0.2, linetype = "blank") +
	facet_wrap(~ treatment, nrow = 2) +ylab("Growth rate (per day)") + xlab("Phosphate (uM P)") + 
	scale_color_manual(values = cols) +
	scale_fill_manual(values = cols)
ggsave("figures/phosphate-monod-bootstrap.png", width = 12, height = 6)
